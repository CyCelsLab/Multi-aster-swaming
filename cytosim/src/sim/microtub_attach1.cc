//RCS: $Id: microtub_attach1.cc,v 2.12 2005/04/19 18:17:41 nedelec Exp $
//======================================================================
//===========================DISPATCH===================================
//======================================================================

// Divide-and-Conquer algorithm to find all MT rods close to any given point
// We first set-up an square grid covering the space: grid of class Map
// each point on this grid contains the list of all rods closer to an
// appropriate threshold H (see below), 
// the function TryToAttach(place, hand...) finds the closest grid-point to the 
// given position, and tries all the rods in the attached list with distance().

// The threshold H is calculated in order not to miss any rod that would
// be closer to MHGL (Maximum Hand Attachment Length = max of all haattachdist[])
// Finding H relies on a simple relation on the distances of 3 points:
// (A) distance( GP, rod ) < distance( GP, motor ) + distance( place, rod )
// where GP (grid-point) is the closest point on the grid to position.
// we want to find all rods for which:  (B) distance( motor, rod ) < MHGL
// because the grid is square, if GP is the closest to motor, we have:
// (C) distance( motor, GP ) < 0.5 * mtGrid.diagonal()
// with (B) (C) in (A), a list should includes all the rods closer than:
//  H =  MP.MHGL + 0.5 * mtGrid.diagonal();

// So, we need to build for all GPs a list of all the rods closer than H!
// the rod will be linked multiple time in different lists.
// This is done in an inverse way: for each rod, we construct a region  
// around the rod by inflating the segment by the width H, then we consider 
// once every GP inside this volume, and link the current rod.

// this last operation is probably similar to algorithms used to paint
// polygons on pixel arrays, which are really fast and smart I suppose. 
// We should find and use them! Here we use the custom written 'painting'
// routines defined in 'rasterizer.h', which call 'paintGrid()' defined below.


//TODO: we could perform a dispatchRods() not at every step, but only if the 
//MT points have moved by a certain threshold. This would work if we also
//extend the painted area around the rod, by the same threshold.
//we must also redo the dispathRods() when MT points are added or removed.
//This new algortihm should be a very big CPU gain on the attachment part

#include "rasterizer.h"

//================================CELLS=====================================
//    divide space into little square regions, used to find close objects
//==========================================================================


//grid for divide and conquer strategies:
Map<DIM, Array<PointMicrotub*> > MicrotubList::mtGrid;

//the current rod being painted:
PointMicrotub * MicrotubList::paint_rod;



const int SAFETY_BORDER = 8;


int MicrotubList::setMTGrid()
{
  //TODO: tune the Grid size for speed
  //in 3D the grid is (nbcell)^3. with nbcell=1000,
  //this is about 2^30, near the limit of 32 bits memory!
  
  real minsize = (DIM==3)? 1.0 : 0.5;
  Vecteur bounds = Microtub::space->getBoundingRect();
  
  Vecteur cell_size;
  int ii, cell_nb[3];
  
  for (ii = 0; ii < DIM; ++ii) {
    //we cannot at the moment handle dimension of zero:
    if ( bounds[ii] == 0 ) {
      MSG.warning("setMTGrid", " size[%d] = 0: Check parameters boxshape and boxsize", ii);
      bounds[ii] = minsize;
      //return 1;
    }
    cell_nb[ii]    = maxT(1, int( 2 * bounds[ii] / minsize ));
    cell_size[ii]  = 2 * bounds[ii] / (real)cell_nb[ii];
  }
  
  deleteMTGrid();
  
  //add a safety border of SAFETY_BORDER cells, in the non periodic dimensions
  //without changing the cell_size
  if ( Space::isPeriodic() ) {
    
    //the last dimension is not periodic
    if ( Microtub::space->getShape() == SHAPE_STRIP ) {
      ii = DIM - 1;
      cell_nb[ii]  += SAFETY_BORDER * 2;
      bounds[ii]   += SAFETY_BORDER * cell_size[ii];
    }
    
    mtGrid.setDimensions(-bounds, bounds, cell_nb);
    mtGrid.setPeriodic();
    
  } else {
    
    for (ii = 0; ii < DIM; ++ii ) {
      cell_nb[ii]  += SAFETY_BORDER * 2;
      bounds[ii]   += SAFETY_BORDER * cell_size[ii];
    }
    mtGrid.setDimensions(-bounds, bounds, cell_nb);
    
  }
  
  mtGrid.allocateCellArray();
  
  MSG(3, "grid XYZ:");
  for( ii = 0; ii < DIM; ++ii )
    MSG(3, "  %.1fum / %i cells", 2*bounds[ii], cell_nb[ii]);
  MSG(3, "\n");
  return NO_ERROR;
}

//----------------------------------------------------------------------
void MicrotubList::deleteMTGrid()
{
  if ( mtGrid.cellsAllocated() ) {
    for (int cell=0; cell < mtGrid.nbCells(); ++cell)
      mtGrid[ cell ].clearFast();
    
    mtGrid.deleteAll();
  }
}

//----------------------------------------------------------------------
bool MicrotubList::isMTGridReady()
{
  return mtGrid.cellsAllocated();
}

//----------------------------------------------------------------------

//paintGrid(x,y,z)  links paint_rod in the list at (x,y,z)
//it is given to paintFatLine(), which calls it for all integer point
//inside the volume around the rod

void MicrotubList::paintGrid( const int x, const int y, const int z )
{
  
#if   ( DIM == 1 )
  mtGrid.cell( x ).pushFront( paint_rod );
#elif ( DIM == 2 )
  mtGrid.cell( x, y ).pushFront( paint_rod );
#elif ( DIM == 3 )
  mtGrid.cell( x, y, z ).pushFront( paint_rod );
#endif
    
  /*
  //---------------- collect statistics on the grid:
  static PointMicrotub * rod = 0;
  static int cnt = 0;

  if ( rod != paint_rod ) {
    if ( rod ) MSG("paintGrid() rod=M%lx %i, cnt %i\n", 
		     rod->mt->getName(), rod->rod, cnt);
    rod = paint_rod;
    cnt = 0;
  }
  ++cnt;
  */
  //MSG("paintGrid(%3i, %3i, %3i)\n", x, y, z);

  /*
  //--------------- visual debugging by adding Grafteds
  // >>>> Do not forget to enable the clearGrafted() below <<<<
  int coo[] = { x, y, z };
  Vecteur place;
  mtGrid.cellWhere( place, coo );
  sim.link( new Grafted( paint_rod->rod%2, place ));
  */
}




//----------------------------------------------------------------------

void MicrotubList::dispatchRods()
{
  assert( true == isMTGridReady() );
  //sim.clearGrafted(); //-------------- in case of visual debugging
  
  Vecteur p, q;

  const real * minV   = mtGrid.leftSide();
  const real * deltaV = mtGrid.cellWidth();

  real width = MP.MHGL + 0.5 * mtGrid.diagonalCellLength();
  //in 3D, we need to inflate the rectangle to contain the cylinder: sqrt(2)
  if ( DIM == 3 ) width *= 1.4143;
  
  //MSG(9, "dispatchRods: MHGL=%f, width=%f\n", MP.MHGL, width);
 
  //----------------------- clear all the lists in the grid:

  //this is equivalent to a simple for(), but a little faster
  Array<PointMicrotub*> * last = &( mtGrid.lastCell() );
  for( Array<PointMicrotub*> * list = & mtGrid[ 0 ]; list <= last ; ++list)
    list->clearFast();
  
  //----------------------- paint them with the MT rods:

  real scaling;

  for( Microtub * mt=firstMicrotub(); mt; mt=mt->next() ) {
      
      p = mt->whereP(0);
      scaling = width / mt->rodLength();

      for( int pt = 1; pt < mt->nbPoints(); ++pt ) {
          
          paint_rod = mt->getRod( pt-1 );
	
          if ( pt % 2 ) {
              q = mt->whereP( pt );
              Space::moduloNear( q, p );
          } else {
              p = mt->whereP( pt );
              Space::moduloNear( p, q );
          }
          
#if   (DIM == 1)
          Rasterizer::paintFatLine1D(paintGrid, p, q, width, minV, deltaV);
#elif (DIM == 2)
          Rasterizer::paintFatLine2D(paintGrid, p, q, width, minV, deltaV, scaling);
#elif (DIM == 3)
          Rasterizer::paintFatLine3D(paintGrid, p, q, width, minV, deltaV, scaling);
#endif
      }
  }
  
  rod_dispatch_step = sim.iterationCnt();  //for the assertion in tryToAttach()
}


//======================================================================
//===========================ATTACH=====================================
//======================================================================

//tryToAttach() returns true if it attached the motor, false if it did not
bool MicrotubList::tryToAttach(const Vecteur & place, Hand& ha )
{
  assert( rod_dispatch_step == sim.iterationCnt() );

  //get the grid node list index closest to the position in space:
  int indx = mtGrid.indexFromPosition_shifted( place.getXYZ() );

  //get the list of rods associated with this cell:
  Array<PointMicrotub*> & rodlist = mtGrid[ indx ];

  //get the number of rods in this list:
  int rodi = rodlist.size();
  if ( rodi <= 0 ) return 0;

  //randomize the list, to make attachments more fair:
  //this might not be necessary, since the MT list is already mixed
  rodlist.mixWell();

  //trod will be our 'target rod':
  PointMicrotub * trod;
  int targetCnt = -1;
  
  real attach_dist = MP.haattachdist[ ha.getType() ];
  
  do {
    --rodi;
    trod = rodlist[ rodi ];
    
    //prevent illegal attachments to be made, base on key match
    //if ( trod->keyMatch( ha.bindingKey() )) continue;
    
    //we compute the distance from the hand to the rod,
    //and compare it to the maximum grabing distance of the hand.
  
    real dist = trod->distanceToRod(place);
    
    if ( dist > attach_dist ) continue;     

#ifdef NEAR_ATTACH_FASTER

    //NEW: probability increase for closer targets, following a 1/R law
    //the random test for binding is eventually done by calling RNG
    if ( dist * RNG.preal() > MP.haattachrate_dt[ ha.getType() ] ) {
      // we can now attach the hand to the target rod:
      if ( ha.attach( *trod ))
        return true;
      else
        ha.clear();
    }
    
#else    

    //the random test for binding is eventually done by calling RNG
    if ( RNG.test( MP.haattachrate_dt[ ha.getType() ] )) {
      // we can now attach the hand to the target rod:
      if ( ha.attach( *trod )) 
        return true;
      else
        ha.clear();
    }
    
#endif
    
    //we limit the number of trial according to the parameter MP.hamaxtargets.
    //With the default value MP.hamaxtargets=1, we give two chances to an hand from
    //a free complex, and one to the remaining hand from an already bound complex.
    //Usually, there is only one MT rod close enough to bind to, in which
    //case the hand has only this one chance irrespective of the complex state.

    //initialize targetCnt according to whether the hand is in a free complex or not:
    //free complex start with targetCnt = 0
    //bound complex start with targetCnt = 1
    //with MP.hamaxtargets == 1, free complex get 2 attempts, bound complex 1 attempt
    if ( targetCnt < 0 )
      targetCnt = ha.otherBound();
    //test if the number of target considered is above the specified limit:
    if ( ++targetCnt > MP.hamaxtargets )  {             //one more target
      //by default, 
      //if above the limit, we return without binding:
      return false;
    }
  } while( rodi > 0 );

  return false;      // we did not bind
}


//============================-=-=-=-=-=-=-=-=-=-===========================
//        accessory (unused) function to find the closest mt point
//--------------------------------------------------------------------------

//find and returns the closest rod, if it is closer than MP.MHGL,
PointMicrotub * MicrotubList::closestPoint(const Vecteur & place)
{
  assert( rod_dispatch_step == sim.iterationCnt() );

  //get the cell index from the position in space:
  int indx = mtGrid.indexFromPosition_shifted( place.getXYZ() );

  //get the list of rods associated with this cell:
  Array<PointMicrotub*> & rodlist = mtGrid[ indx ];

  PointMicrotub * trod, * best_rod = 0;
  real dist, best_dist = 2 * MP.MHGL;

  for( int rodi = 0; rodi < rodlist.size(); ++rodi )
    {
      trod = rodlist[ rodi ];

      //we compute the distance from the hand to the candidate rod,
      //and compare it to the best we have so far.
      dist = trod->distanceToRod( place );
      if ( dist < best_dist ) {
          best_dist = dist;
          best_rod = trod;
      }
    }
  return best_rod;
}
