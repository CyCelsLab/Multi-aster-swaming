//RCS: $Id: space_combine.cc,v 2.0 2004/08/16 16:38:07 nedelec Exp $

#include "space_combine.h"


//--------------------------------------------------------------------
bool  SpaceCombine::isInside( const real point[] ) const
{
  return ( outerSpace->isInside( point ) && ! innerSpace->isInside( point ) );
}



//--------------------------------------------------------------------
void SpaceCombine::project( const real point[], real proj[] ) const
{
  static real proj2[DIM];
  
  innerSpace->project( point, proj );

  if ( outerSpace -> isConfined() ) {
    
    outerSpace->project( point, proj2 );

    real pp = 0;
    for( int dd = 0; dd < DIM; ++dd )
      pp += ( point[dd] -  proj[dd] ) * ( point[dd] -  proj[dd] )
          - ( point[dd] - proj2[dd] ) * ( point[dd] - proj2[dd] );

    if ( pp > 0 )
      for( int dd = 0; dd < DIM; ++dd )
        proj[dd] = proj2[dd];
  }
}
