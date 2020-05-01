//RCS: $Id: pointset.cc,v 2.18 2005/04/22 12:13:14 nedelec Exp $
//-------------------------------pointset.cc----------------------------------

#include "pointset.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "iomessages.h"

//----------------------------------------------------------------------------
void PointSet::pointsetConstructorBase()
{
  psallocated    = 0;
  pspts          = 0;
  psmaxpts       = 0;
}

//----------------------------------------------------------------------------
PointSet::PointSet()
{
  pointsetConstructorBase();
}

//----------------------------------------------------------------------------
PointSet::PointSet(const PointSet & o)
{
  pointsetConstructorBase();
  allocate( o.nbPoints() );
  for (int p = 0; p < DIM * psmaxpts; ++p )
    pspts[ p ] = o.pspts[ p ];
}

//----------------------------------------------------------------------------
/** allocate(size) ensures that the set can hold <size> points
 it returns the size if new memory was allocated 
*/
int PointSet::allocate(const int size)
{
  if ( size <= psallocated ) 
    return 0;

  int psallocated_new = size + PSCHUNCK;

  real  *  pspts_new = new real [ DIM*psallocated_new ];

  if ( 0 == pspts_new )
    MSG.error("PointSet::allocate()", "memory allocation failed");

  if ( pspts ) {   // copy the current position to the new array
    for (int p = 0; p < DIM * psallocated; ++p )
      pspts_new[ p ] = pspts[ p ];
    delete[] pspts;
  }
  pspts        = pspts_new;
  psallocated  = psallocated_new;
  
  MSG(12, "PointSet::allocate %lx ask %d done %d\n", this, size, psallocated_new);
  return psallocated_new;
}


//---------------------------------------------------------------------------
void PointSet::deallocate()
{
  if ( pspts )
    delete[] pspts;   
  pspts = 0;
  psallocated = 0;
}

//----------------------------------------------------------------------------
void PointSet::setPoint(const int p, const Vecteur & w )
{
  assert(( p >= 0 ) && ( p < psallocated ));

                   pspts[ p * DIM     ] = w.XX;
  if ( DIM > 1 )   pspts[ p * DIM + 1 ] = w.YY;
  if ( DIM > 2 )   pspts[ p * DIM + 2 ] = w.ZZ;
}

//----------------------------------------------------------------------------
int PointSet::addPoint( const Vecteur & w )
{
  allocate( psmaxpts+1 );
  setPoint( psmaxpts++, w );
  return psmaxpts-1;
}

//----------------------------------------------------------------------------
void PointSet::movePoint(const int p, const Vecteur & w )
{
  assert(( p >= 0 ) && ( p < psmaxpts ));

                   pspts[ p * DIM     ] += w.XX;
  if ( DIM > 1 )   pspts[ p * DIM + 1 ] += w.YY;
  if ( DIM > 2 )   pspts[ p * DIM + 2 ] += w.ZZ;
}

//----------------------------------------------------------------------------
void PointSet::resetPoints()
{
  assert( pspts );
  for( int p = 0; p < DIM * psallocated; ++p )
    pspts[ p ] = 0;
}

//----------------------------------------------------------------------------
void PointSet::addNoise( const real howmuch )
{
  for( int p = 0; p < DIM * psmaxpts; ++p )
    pspts[ p ] += howmuch * RNG.sreal();
}

//----------------------------------------------------------------------------
void PointSet::translatePosition( const Vecteur & T )
{
  for( int p = 0; p < DIM * psmaxpts;  ) {
                     pspts[ p++ ] += T.XX;
    if ( DIM > 1 )   pspts[ p++ ] += T.YY;
    if ( DIM > 2 )   pspts[ p++ ] += T.ZZ;
  }
}

//----------------------------------------------------------------------------
void PointSet::transformPosition( const Transformation & T )
{
  for( int p = 0; p < DIM * psmaxpts; p += DIM) 
    T.rightMultOverride( &pspts[ p ] );
}

//----------------------------------------------------------------------------
Vecteur PointSet::getPosition() const
{
  Vecteur result = VZERO;
  for(int p = 0; p < DIM * psmaxpts; ) {
                   result.XX += pspts[ p++ ];
    if ( DIM > 1 ) result.YY += pspts[ p++ ];
    if ( DIM > 2 ) result.ZZ += pspts[ p++ ];
  }
  if ( psmaxpts > 1 ) result /= real( psmaxpts );
  return result;
}

//----------------------------------------------------------------------------
int PointSet::insidePosition( const Space * s ) const
{
  int result=0;
  for(int p = 0; p < psmaxpts; ++p )
    result += s->isInside( pspts+DIM*p );
  return result;
}


//----------------------------------------------------------------------------
/** modulo, around a middle point */
void PointSet::moduloPosition()
{
  int middle = psmaxpts / 2;
  Space::modulo( pspts+DIM*middle );
  for(int p = middle+1; p < psmaxpts; ++p )
    Space::moduloNear( pspts+DIM*p, pspts+DIM*(p-1));
  for(int p = middle-1; p >= 0; --p )
    Space::moduloNear( pspts+DIM*p, pspts+DIM*(p+1));
}

//----------------------------------------------------------------------------
/**modulo around the center of gravity: use this for compact objects,
this is not a good option for linear objects (MTs)*/
void PointSet::moduloPositionG()
{
  calculatePosition();
  Space::modulo( pscenter );
  for(int p = 0; p < psmaxpts; ++p )
    Space::moduloNear( pspts+DIM*p, pscenter );
}


//****************************************************************************
//----------------------------------------------------------------------------
void PointSet::write()
{
  IO.writeUInt16( psmaxpts );
  for( int p = 0; p < psmaxpts ; ++p )
    IO.writeReal32Vect( DIM, pspts + DIM * p, true );
}



//----------------------------------------------------------------------------
void PointSet::read()
{
  try {
    int newmaxpts = IO.readUInt16();
    MSG(80, "PointSet::read( s%lx ) maxpts %i new %i\n", this, psmaxpts, newmaxpts );
    allocate( newmaxpts );
    resetPoints();        //not absolutely necessary, but nicer for debugging

    psmaxpts = newmaxpts;
    for( int p = 0; p < newmaxpts ; ++p )
      IO.readReal32Vect( pspts+DIM*p );

    //we do not call modulo, as this is dependent on the objects
    
  } catch( IOException e ) {
    //MSG("S%lx %i ", name, r);
    e.addBeforeMessage("PointSet::read : ");
    psmaxpts = 0;
    throw e;
  }
}
