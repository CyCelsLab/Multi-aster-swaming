//RCS: $Id: object.cc,v 1.7 2005/04/22 11:50:27 nedelec Exp $
/*
 *  object.cc
 *  cytosim
 *
 *  Created by Francois Nedelec on Wed Apr 07 2004.
 *  Copyright (c) 2004 EMBL All rights reserved.
 *
 *  Modified by Chaitanya Athale on Thu Feb 09 2006.
 *  Copyright (c) 2006 EMBL Heidelberg. All rights reserved.
 */

#include "object.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "sim_param.h"


//*************************************************************************
//========================DEFAULT DISTRIBUTIONS============================
//*************************************************************************

//TODO: we should move all the initial configurations to an accessory program
//this would keep the code of the main simulation leaner and cleaner
Vecteur Object::initPosition(const Space * space, const int mode)
{
  assert( space && space->isSet() );
  unsigned long ouf = 0;
  Vecteur place;

  switch ( mode ) {

    case 0:                       //anywhere
      return space->randomPlaceInVolume();

    case 1:                       //inside a central sphere / disc
      return Vecteur::randSphere(MP.initsize[0]);

    case 2:                       //at center
      return VZERO;

    case 3:                       //outside a central sphere / disc
      do {
        place = space->randomPlaceInVolume();
        if (++ouf>MAX_TRIALS_BEFORE_STUCK)
          throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::position()");
      } while ( place.normSquare() < MP.initsize[0]*MP.initsize[0] );
      return place;

    case 4:                       //surface of a central sphere / disc
      return ( MP.initsize[0] + MP.boxedge * RNG.preal() ) * Vecteur::randNormed();
      

    case 5:                       //a vertical stripe in center
      do {
        place = space->randomPlaceInVolume();
        if (++ouf>MAX_TRIALS_BEFORE_STUCK)
          throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::position()");
      } while ( fabs(place.XX) >= MP.initsize[0] );
      return place;

    case 6:                       //outside the vertical stripe in center
      do {
        place = space->randomPlaceInVolume();
        if (++ouf>MAX_TRIALS_BEFORE_STUCK)
          throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::position()");
      } while ( fabs(place.XX) < MP.initsize[0] );
      return place;

    case 7:                       //ellipse in the middle
      do {
        place = space->randomPlaceInVolume();
        for( int ii = 0; ii < 3; ++ii)
          if ( MP.initsize[ii] == 0 ) place[ii] = 0;
        if (++ouf>MAX_TRIALS_BEFORE_STUCK)
          throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::position()");
      } while ( ! insideEllipse( place, MP.initsize ));
      return place;

    case 8:                       //outside the ellipse in the middle
      do {
        place = space->randomPlaceInVolume();
        if (++ouf>MAX_TRIALS_BEFORE_STUCK)
          throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::position()");
      } while ( insideEllipse( place, MP.initsize ));
      return place;

    case 9:           //the edge of the box
      if ( space->isConfined() )
        return space->randomPlaceOnEdge( MP.boxedge );
      else
        return space->randomPlaceInVolume();

    case 10:          //outside the ellipse in the middle & inside sphere
      do {
        place = space->randomPlaceInVolume();
        if (++ouf>MAX_TRIALS_BEFORE_STUCK)
          throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::position()");
      } while (insideEllipse(place, MP.initsize) || (place.norm() > MP.initsize[3]));
      return place;

    case 11:          //inside a central sphere / disc
      return Vecteur::randSphere(MP.initsize[3]);


    case 12:          //inside the ellipse, outside the central spot
      do {
        place = space->randomPlaceInVolume();
        for( int ii = 0; ii < 3; ++ii)
          if ( MP.initsize[ii] == 0 ) place[ii] = 0;
        if (++ouf>MAX_TRIALS_BEFORE_STUCK)
          throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::position()");
      } while (!insideEllipse(place, MP.initsize) || (place.norm() < MP.initsize[3]));
      return place;

    case 13:        //two patches at the surface of a central sphere
      do {
        place = ( MP.initsize[0] + MP.boxedge * RNG.preal() ) * Vecteur::randNormed();
      } while ( place.YY*place.YY + place.ZZ*place.ZZ > 1 );
      return place;

    case 15:                       //a horizontal stripe in center
      do {
        place = space->randomPlaceInVolume();
        if (++ouf>MAX_TRIALS_BEFORE_STUCK)
          throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::position()");
      } while ( fabs(place.YY) > MP.initsize[0] );
      return place;

    case 16:                       //CHAITANYA: a VERTICAL stripe in center- with random positions on Y-axis
	      do {
	        place = space->randomPlaceInVolume();
	        if (++ouf>MAX_TRIALS_BEFORE_STUCK)
	          throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::position()");
	      } while ( fabs( place.XX ) > MP.initsize[0] );
      return place;

    case 20:
    case 21: {                    //a gradient in the X (15) or Y (17) direction
      real offset, slope, prob;
      if ( mode == 20 )
        offset = space->getBoundingRect().XX;
      else
        offset = space->getBoundingRect().YY;

      slope  = ( MP.boxgradient[1] - MP.boxgradient[0] ) /
        ( 10 * ( MP.boxgradient[1] + MP.boxgradient[0] ) * offset );

      do {
        place = space->randomPlaceInVolume();
        if ( mode == 20 )
          prob = slope * ( place.XX + offset );
        else
          prob = slope * ( place.YY + offset );
        assert( prob < 0.2 );
      } while ( !RNG.test( prob ));
    } return place;
    
    
    
	// neha khetan: 19 july 2014  : case 23 and case 24, 25 added
	// For MTOC motility simulation wherein, some % of asters are to be initialized at the NE (case23) and the remaining in the cytoplasm (case 24).
	// case 25 is buffer wherein, under case 25 option, depending on the % set: the asters enter into case 23 or case 24.
     
    case 23:                  			 // on the surface of the nucleus
      return ( MP.nuradius + MP.boxedge * RNG.preal() ) * Vecteur::randNormed();
      
    
    case 24:                       		//outside nucleus in the ooplasm
      do {
        place = space->randomPlaceInVolume();
        if (++ouf>MAX_TRIALS_BEFORE_STUCK)
          throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::position()");
      } while ( place.normSquare() < MP.nuradius*MP.nuradius
       );
      return place;
      
      
     case 25:

	// end of nk modification on 22 July 2014 ......................................................................................................
	
	
	
	case 26:
	
		{
		/*
		int boolPosgh = 0;
		
			do{
				place = Grafted::space->randomPlaceInVolume();
				real pRgh, pNgh; 					//random and gh-positioning probabilities					
				place = sim.ghNucl();			    // Retrurns the coordinates random in volume// 
				pRgh = RNG.real_uniform_range(0, 1);
				pNgh =	sim.ghMap( place[0] , place[1] );		//Calculating the prob. of positioning a grafted  //pNgh =	sim.fieldModel( posGh[0],posGh[1]);
				if (pRgh < pNgh )
					 {
						boolPosgh = 1;
					 }
				}while ( boolPosgh == 0);
				
			return place;
	
	
 */


    case 27:                       //inside nucleus
      return Vecteur::randSphere(MP.nuradius);



    }default:
      MSG.error("Space::Position","unknowned init code");
      return VZERO;
  }
}

//-------------------------------------------------------------------------CHAITANYA: case 22: in order to initially position
// objects in a specific location given here through the data.in file. Eventually by possibly reading in the data
// image file and automatically fixing position in space.
Vecteur Object::initPosition(const Space * space, const int mode, const Vecteur initPos)
{

  assert( space && space->isSet() );
  Vecteur place;
  unsigned long ouf = 0;

  switch ( mode ) {

	  case 22://the specific position
	   printf("nuVect x: %f y: %f z: %f \n", initPos[0], initPos[1], initPos[2]);
	    place = Vecteur( initPos[0], initPos[1], initPos[2]);
	    return place;

  default:
     MSG.error("Space::Position","unknown init code");
  return VZERO;

     }

}

//-------------------------------------------------------------------------
Vecteur Object::initPosition(const Space * space, const int mode, int & type)
{
  assert( space && space->isSet() );
  Vecteur place;

  switch ( mode ) {

    case 30: {                      //some horizontal wavy lines / spirals in 3D
      type = RNG.pint_inc(2);
      real offset = MP.initsize[2] * ( RNG.sflip() * type );
      place.XX = RNG.sreal() * MP.boxsize[0];
      place.YY = ( offset + sin( PI * place.XX / MP.initsize[0] )) * MP.initsize[1];
      //place.ZZ = cos( PI * place.XX / MP.initsize[0] ) * MP.initsize[1];
    } return place;

    case 31: {                      //skating-tracks
      type = RNG.pint_inc(3);
      real offset = MP.initsize[2] * ( type - 1.5 );
      place.XX = RNG.sreal() * MP.boxsize[0];
      real slope = MP.initsize[1] / MP.boxsize[0];
      if ( place.XX > 0 )
        place.YY = offset - place.XX * slope;
      else
        place.YY = offset + (MP.boxsize[0] + place.XX) * slope;
    } return place;

    case 32: {                      //skating-tracks, 4 motors
      type = RNG.pint_inc(3);
      real offset = MP.initsize[2] * ( type - 1.5 );
      int s = RNG.sflip();
      real slope = MP.initsize[1] / MP.boxsize[0];
      place.XX = 2 + RNG.preal() * ( MP.boxsize[0] - 2 );
      place.YY = offset + place.XX * s * slope;
      if ( s < 0 ) place.XX -= MP.boxsize[0];
    } return place;

    case 33: {                      //skating-tracks, 2 motors
      type = RNG.pint_inc(3);
      real offset = MP.initsize[2] * ( type - 1.5 );
      int s = RNG.sflip();
      real slope = MP.initsize[1] / MP.boxsize[0];
      place.XX = 2 + RNG.preal() * ( MP.boxsize[0] - 2 );
      place.YY = offset + place.XX * s * slope;
      if ( s < 0 ) place.XX -= MP.boxsize[0];
      type = (type % 2 ? 1 : 0 );
    } return place;

    case 34: {                      //skating-tracks, drop-shaped
      real arm = MP.initsize[0];
      real rad = MP.initsize[1];
      type = RNG.flip();
      int side = RNG.sflip();
      real rot = side * MP.initsize[2];
      real offX = ((side>0)?0:-MP.boxsize[0]);
      real offY = ((type==0)?-0.5:0.5) * MP.boxsize[1];
      real len = sqrt( arm*arm + rad*rad ); //approximate formula

      if ( RNG.test( PI*rad / ( PI*rad + 2*len ))) {
        Vecteur center( offX + arm * cos(rot), offY + arm * sin(rot), 0 );
        rot += ( RNG.preal() -0.5 ) *PI;
        return center + Vecteur( rad * cos(rot), rad * sin(rot), 0 );
      }

      real abs = RNG.preal() * len;
      real ang = atan2( rad, arm );

      if ( RNG.flip() )
        return Vecteur( offX + abs * cos( rot-ang ), offY + abs * sin( rot-ang ), 0 );
      else
        return Vecteur( offX + abs * cos( rot+ang ), offY + abs * sin( rot+ang ), 0 );
    }

    case 35: {                      //skating-tracks, with curl at the end
      real arm = MP.initsize[0];
      real rnd = MP.initsize[1];
      type = RNG.flip();
      int side = RNG.sflip();
      real rot = side * MP.initsize[2];
      real rad = MP.initsize[3];

      real offX = ((side>0)?0:-MP.boxsize[0]);
      real offY = ((type==0)?-0.5:0.5) * MP.boxsize[1];

      if ( RNG.test( rnd / ( rnd + arm ))) {
        Vecteur center( offX + arm*cos(rot) - side*rad*sin(rot), offY + arm*sin(rot) + side*rad*cos(rot), 0 );
        rot = rot -side * PI * 0.5 + side * RNG.preal() * rnd / rad;
        return center + Vecteur( rad * cos(rot), rad * sin(rot), 0 );
      }

      real abs = RNG.preal() * arm;
      return Vecteur( offX + abs * cos(rot), offY + abs * sin(rot), 0 );
    }

  default:
    return initPosition( space, mode );
}

}

