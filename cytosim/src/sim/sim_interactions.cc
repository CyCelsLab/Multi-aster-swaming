//RCS: $Id: sim_interactions.cc,v 2.17 2005/03/14 18:16:19 nedelec Exp $

//==========================================================================
//The accessory functions 'near()' below are used to check the validity of
//the arguments givent to the functions interactions().
//two PointExact should point to different point,
//two PointInterpolated should not share any common point. 
//in any other case, the matrix elements would be updated in a wrong way
//if near() is true, interaction() should exit, and not update the matrix.

bool near( const PointExact & pta, const PointExact & ptb ) {
  assert( pta.looksWrong(true) == NO_ERROR );
  assert( ptb.looksWrong(true) == NO_ERROR );
  if ( pta.near(ptb) ) {
    //MSG.warning("SIMsolve::interactionLink", "IntraLink EE");
    return true;
  } else {
    return false;
  }
}
  
bool near( const PointInterpolated & pta, const PointExact & ptb ) {
  assert( pta.looksWrong(true) == NO_ERROR );
  assert( ptb.looksWrong(true) == NO_ERROR );
  if ( pta.near(ptb) ) {
    //MSG.warning("SIMsolve::interactionLink", "IntraLink IE");
    return true;
  } else {
    return false;
  }
}

bool near( const PointInterpolated & pta, const PointInterpolated & ptb ) {
  assert( pta.looksWrong(true) == NO_ERROR );
  assert( ptb.looksWrong(true) == NO_ERROR );
  if ( pta.near(ptb) ) {
    //MSG.warning("SIMsolve::interactionLink", "IntraLink II");
    return true;
  } else {
    return false;
  }
}


//==========================================================================
//$$$$$$$$$$$$$$$$$$$$$$  FORCE MATRIX ELEMENTS  $$$$$$$$$$$$$$$$$$$$$$$$$$$
//==========================================================================


//==========================================================================
//update the matrix mB to include an interaction between points aa and bb
//the force is linear with a zero resting length: f_a = weight * ( b - a )
void interactionLink( const PointExact & pta, const PointExact & ptb, 
                      const real weight = 1.0 )
{
  if ( near( pta, ptb ) ) return;
  
  int aa = pta.matIndex();        assert(( aa >= 0 ) && ( aa < nbpts ));
  int bb = ptb.matIndex();        assert(( bb >= 0 ) && ( bb < nbpts ));
  assert( aa != bb );
  
  mB( aa, aa) -= weight;
  mB( aa, bb) += weight;
  mB( bb, bb) -= weight; 
  
  if ( Space::isPeriodic() ) {
    Vecteur off;
    for( int dd = 0; dd < DIM; ++dd )
      off[ dd ] = vPTS[ DIM*aa + dd ] - vPTS[ DIM*bb + dd ];
    Space::offset( off );
    for( int dd = 0 ; dd < DIM; ++dd ) {
      if ( off[dd] ) {
        vBAS[ DIM*aa + dd ] += weight * off[dd];
        vBAS[ DIM*bb + dd ] -= weight * off[dd];
      }
    }
  }
}

//==========================================================================
//update the matrix to include an interaction between a Microtubule position,
//and another point from a pointset with the given weight
void interactionLink( const PointInterpolated & aa, const PointExact & ptb, 
                      const real weight = 1.0 )
{
  if ( near( aa, ptb ) ) return;
  
  //coefficients on the points:
  real x[] = { 1-aa.getCoef(), aa.getCoef(), -1 };
  
  //the index of the points in the matrix mB:
  int index[] = { aa.matIndex1(), aa.matIndex2(), ptb.matIndex() };
  
  for( int kk = 0;  kk < 3; ++kk )
    for( int ll = kk; ll < 3; ++ll )
      mB( index[ kk ], index[ ll ]) -= weight * x[ kk ] * x[ ll ];
  
  if ( Space::isPeriodic() ) {
    Vecteur off;
    for( int dd = 0 ; dd < DIM; ++dd )
      off[dd] = vPTS[ DIM*index[0] + dd ] - vPTS[ DIM*index[2] + dd ];
    Space::offset( off );
    for( int dd = 0 ; dd < DIM; ++dd ) {
      if ( off[dd] )
        for( int kk = 0; kk < 3 ; ++kk )
          vBAS[ DIM*index[ kk ] + dd ] += weight * x[ kk ] * off[dd];
    }
  }
}

//==========================================================================
//update the matrix mB to include an interaction between PointInterpolated a and b
//the force is linear with a zero resting length: f_a = ( b - a )
void interactionLink( const PointInterpolated & aa, const PointInterpolated & bb,
                      const real weight = 1.0 )
{
  if ( near( aa, bb ) ) return;
  
  //coefficients on the points:
  real h[] = { 1-aa.getCoef(), aa.getCoef(),  bb.getCoef()-1, -bb.getCoef() };

  //coefficients on the points:
  real wh[] = { weight*h[0], weight*h[1],  weight*h[2], weight*h[3] };

  //the index of the points in the matrix mB:
  int index[] = { aa.matIndex1(), aa.matIndex2(), bb.matIndex1(), bb.matIndex2() };
 
  for( int kk = 0;  kk < 4; ++kk )
    for( int ll = kk; ll < 4; ++ll )
      mB( index[ kk ], index[ ll ]) -= wh[ kk ] * h[ ll ];

  if ( Space::isPeriodic() ) {
    Vecteur off;
    for( int dd = 0; dd < DIM; ++dd )
      off[ dd ] = vPTS[ DIM*index[0]+dd ] - vPTS[ DIM*index[2]+dd ];
    Space::offset( off );
    for( int dd = 0 ; dd < DIM; ++dd ) {
      if ( off[ dd ] )
        for( int kk = 0; kk < 4 ; ++kk )
          vBAS[ DIM*index[ kk ] + dd ] += wh[ kk ] * off[ dd ];
    }
  }
}

//==========================================================================
//update the matrix to include an interaction between points aa and bb,
//the force is linear with non-zero resting length: f_a = ab * ( len / |ab| - 1 )
void interactionLongLink( const PointExact & aa, const PointExact & bb, real len, const real weight =  1.0 )
{
  if ( near( aa, bb ) ) return;
  
  int indexA = DIM * aa.matIndex();  // coef on aa is +weight
  int indexB = DIM * bb.matIndex();  // coef on bb is -weight

  assert( indexA != indexB );
  
  Vecteur ab;
  for(int dd = 0; dd < DIM; ++dd )
    ab[dd] = vPTS[ indexB+dd ] - vPTS[ indexA+dd ];

  Vecteur offset;
  if ( Space::isPeriodic() )
    Space::moduloAndOffset( ab, offset );
  
  real abn = ab.norm();
  if ( abn < EPSILON ) return;
  ab /= abn;

  for(int dd = 0; dd < DIM; ++dd ) {
    vBAS[ indexA + dd ] -= weight * ab[dd] * len;
    vBAS[ indexB + dd ] += weight * ab[dd] * len;
  }
  
  real m;
  len /= abn;

  //this is needed to stabilize the long links under compression
  //I do not have a good mathematical explanation of that fact...
  //this correspond to using (len==1) in the formulae if (len > 1.0)...
  bool cooked = ( len > 1.0 );
  
  for(int ii = 0; ii < DIM; ++ii ) {

    //diagonal elements
    if ( cooked )
      m = weight * ab[ii] * ab[ii]; //formula below with len=1
    else
      m = weight * ( 1.0 + len * ( ab[ii] * ab[ii] - 1.0 ));
    
    mC( indexA + ii, indexA + ii ) -= m;
    mC( indexA + ii, indexB + ii ) += m;
    mC( indexB + ii, indexB + ii ) -= m;
    
    //off-diagonal elements 
    for(int jj = ii+1; jj < DIM; ++jj ) {
      if ( cooked )
        m = weight * ab[ii] * ab[jj];
      else
        m = weight * len * ab[ii] * ab[jj];
      mC( indexA + ii, indexA + jj ) -= m;
      mC( indexA + ii, indexB + jj ) += m;
      mC( indexA + jj, indexB + ii ) += m;
      mC( indexB + ii, indexB + jj ) -= m;
    }
  }
  
  if ( Space::isPeriodic() ) {
    real s = offset * ab;  //scalar product
    for(int ii = 0; ii < DIM; ++ii ) {
      if ( cooked )
        m = weight * s * ab[ii];
      else
        m = weight * ( len * s * ab[ii] + ( 1.0 - len ) * offset[ii] );
      vBAS[ indexA + ii ] -= m;
      vBAS[ indexB + ii ] += m;
    }
  }
}

//==========================================================================
//--------------------------------------------------------------------------
//update the matrices to include an interaction with non-zero resting length
// between the two microtubules points a and b
void interactionLongLink( const PointInterpolated & aa, const PointInterpolated & bb, real len )
{
  if ( near( aa, bb ) ) return;
  
  //force coefficients on the points:
  real x[] = { 1.0-aa.getCoef(), aa.getCoef(),  bb.getCoef()-1.0, -bb.getCoef() };
  
  //index in the matrix mC:
  int index[] = { DIM*aa.matIndex1(), DIM*aa.matIndex2(), DIM*bb.matIndex1(), DIM*bb.matIndex2() };
  
  Vecteur ab = bb.where() - aa.where();
  
  Vecteur offset;
  if ( Space::isPeriodic() )
    Space::moduloAndOffset( ab, offset );
  
  real abn = ab.norm();
  if ( abn < EPSILON ) return;
  ab /= abn;

  int dk, ii, jj, kk, ll;

  for( ii = 0; ii < 4; ++ii )
    for( dk = 0; dk < DIM; ++dk )
      vBAS[ index[ii] + dk ] -= x[ii] * ab[dk] * len;

  real m;
  len /= abn;

  //this is needed to stabilize the long links under compression
  //I do not have a good mathematical explanation of why that works...
  //this correspond to using (len==1) in the formulae if (len > 1.0)...
  bool cooked = ( len > 1.0 );
  
  for( ii = 0; ii < DIM; ++ii ) {
    if ( cooked )
      m = ab[ii] * ab[ii];
    else
      m = 1.0 + len * ( ab[ii] * ab[ii] - 1.0 );
    
    for( kk = 0; kk < 4; ++kk )
      for( ll = kk; ll < 4; ++ll )
        mC( index[kk] + ii, index[ll] + ii ) -= x[kk] * x[ll] * m;
    
    for( jj = ii+1; jj < DIM; ++jj ) {
      if ( cooked )
        m = ab[ii] * ab[jj];
      else
        m = len * ab[ii] * ab[jj];
      for( kk = 0; kk < 4; ++kk )
        for( ll = 0; ll < 4; ++ll )
          mC( index[kk] + ii, index[ll] + jj ) -= x[kk] * x[ll] * m;
    }
  }
  
  if ( Space::isPeriodic() ) {
    real s = offset * ab;  //scalar product
    for( ii = 0; ii < DIM; ++ii ) {
      if ( cooked )
        m = s * ab[ii];
      else
        m = len * s * ab[ii] + ( 1.0 - len ) * offset[ii];
      if ( m ) for( jj = 0; jj < 4; ++jj )
        vBAS[ index[jj] + ii ] -= x[jj] * m;
    }
  }
}

//==========================================================================
//--------------------------------------------------------------------------
// update the matrices to include an interaction with non-zero resting length
// between a microtubule at position a, and another pointExact 
void interactionSideLink( const PointInterpolated & pta, const PointExact & ptb, const real len )
{
#if ( DIM == 2 )
  if ( near( pta, ptb ) ) return;
  assert( Space::isPeriodic() == false );
  
  int pt = ptb.matIndex();
  //force coefficients on the points:
  real aa = pta.getCoef(), bb = 1.0 - aa;
  
  int indexP = DIM * pt;

  Vecteur ga = pta.where();
  ga[0] -= vPTS[ indexP   ];
  ga[1] -= vPTS[ indexP+1 ];
  
  real ee;
  if (( ga ^ pta.dirInter() ) > 0 ) {
    ee =  len / pta.length();
  } else {
    ee = -len / pta.length();
  }
  real eeee = ee * ee;
  
  //index in the matrix mB:
  int index = pta.matIndex1();   
  
  //we put the isotropic terms in mB
  mB( index,   index   ) -=   bb * bb + eeee;
  mB( index,   index+1 ) +=  -aa * bb + eeee;
  mB( index+1, index+1 ) -=   aa * aa + eeee;

  mB( pt,      pt )  -= 1;
  mB( index,   pt )  += bb;
  mB( index+1, pt )  += aa;
  
  //index in the matrix mC:
  index *= DIM;
    
  mC(index  , index+3  ) += ee;
  mC(index+1, index+2  ) -= ee;
  
  mC(index  , indexP+1 ) -= ee;
  mC(index+1, indexP   ) += ee;
  mC(index+2, indexP+1 ) += ee;
  mC(index+3, indexP   ) -= ee;
#else
  MSG.error("Stifflink","only implemented in 2D");
#endif
}

//==========================================================================
//--------------------------------------------------------------------------
// update the matrices to include an interaction with non-zero resting length
// The second pointinterpolated has zero extension 
// TODO: weight argument in interactionAsymmetricSideLink()
void interactionAsymmetricSideLink( const PointInterpolated & pt1, const PointInterpolated & pt2, const real len )
{
  if( Space::isPeriodic() )
    MSG.error("interactionSideLink", "periodic Space not implemented");
  
#if ( DIM == 2 )
  if ( near( pt1, pt2 ) ) return;
  
  real aa1 = pt1.getCoef(), bb1 = 1.0 - aa1;
  real aa2 = pt2.getCoef(), bb2 = 1.0 - aa2;
  
  Vecteur stretch = pt2.where() - pt1.where();
  
  real ee1;
  if ((stretch ^ pt1.dirInter()) < 0 )
    ee1 =  len / pt1.length();
  else
    ee1 = -len / pt1.length();
    
  //index in the matrix mB:
  int index1 = pt1.matIndex1();
  int index2 = pt2.matIndex1();
  
  //we put the isotropic terms in mB
  real ee1ee1 = ee1 * ee1;
  mB( index1,   index1   ) -=  bb1 * bb1 + ee1ee1;
  mB( index1,   index1+1 ) += -aa1 * bb1 + ee1ee1;
  mB( index1+1, index1+1 ) -=  aa1 * aa1 + ee1ee1;
  
  mB( index2,   index2   ) -=  bb2 * bb2;
  mB( index2,   index2+1 ) += -aa2 * bb2;
  mB( index2+1, index2+1 ) -=  aa2 * aa2;
  
  mB( index1  , index2   ) +=  bb1*bb2;
  mB( index1  , index2+1 ) +=  bb1*aa2;
  mB( index1+1, index2+1 ) +=  aa1*aa2;
  mB( index1+1, index2   ) +=  aa1*bb2;
  
  //index in the matrix mC:
  index1 *= DIM;
  index2 *= DIM;
  
  mC( index1  , index1+3 ) += ee1;
  mC( index1+1, index1+2 ) -= ee1;
  
  real ee1bb2 = ee1*bb2;
  real ee1aa2 = ee1*aa2;
  
  mC( index1  , index2+1 ) -=  ee1bb2;
  mC( index1  , index2+3 ) -=  ee1aa2;
  
  mC( index1+1, index2   ) +=  ee1bb2;
  mC( index1+1, index2+2 ) +=  ee1aa2;
  
  mC( index1+2, index2+1 ) +=  ee1bb2;
  mC( index1+2, index2+3 ) +=  ee1aa2;
  
  mC( index1+3, index2   ) -=  ee1bb2;
  mC( index1+3, index2+2 ) -=  ee1aa2;
#else
  MSG.error("AsymmetricSideLink","only implemented in 2D");
#endif
}

//==========================================================================
//update the matrix to include a force 
void interactionSideLink( const PointInterpolated & pt1, const PointInterpolated & pt2, 
                            real len, const real weight=1.0 )
{
  if( Space::isPeriodic() )
    MSG.error("interactionSideLink", "periodic Space not implemented");
  
#if ( DIM == 2 )
  if ( near( pt1, pt2 ) )  return;
  
  real aa1 = pt1.getCoef(), bb1 = 1.0 - aa1;
  real aa2 = pt2.getCoef(), bb2 = 1.0 - aa2;

  Vecteur stretch = pt2.where() - pt1.where();
  
  real ee1;
  if ((stretch ^ pt1.dirInter()) < 0 )
    ee1 =  len / ( 2 * pt1.length() );
  else
    ee1 = -len / ( 2 * pt1.length() );
  
  real ee2;
  if ((stretch ^ pt2.dirInter()) < 0 )
    ee2 = -len / ( 2 * pt2.length() );
  else
    ee2 =  len / ( 2 * pt2.length() );
  
  //index in the matrix mB:
  int index1 = pt1.matIndex1();
  int index2 = pt2.matIndex1();
  
  //we put the isotropic terms in mB
  real ee1ee1 = ee1 * ee1;
  mB( index1,   index1   ) -=  bb1 * bb1 + ee1ee1;
  mB( index1,   index1+1 ) += -aa1 * bb1 + ee1ee1;
  mB( index1+1, index1+1 ) -=  aa1 * aa1 + ee1ee1;

  real ee2ee2 = ee2 * ee2;
  mB( index2,   index2   ) -=  bb2 * bb2 + ee2ee2;
  mB( index2,   index2+1 ) += -aa2 * bb2 + ee2ee2;
  mB( index2+1, index2+1 ) -=  aa2 * aa2 + ee2ee2;

  
  real ee1ee2 = ee1 * ee2;
  mB( index1  , index2   ) +=  bb1*bb2 + ee1ee2;
  mB( index1  , index2+1 ) +=  bb1*aa2 - ee1ee2;
  mB( index1+1, index2+1 ) +=  aa1*aa2 + ee1ee2;
  mB( index1+1, index2   ) +=  aa1*bb2 - ee1ee2;

  //index in the matrix mC:
  index1 *= DIM;
  index2 *= DIM;
  
  mC( index1  , index1+3 ) += ee1;
  mC( index1+1, index1+2 ) -= ee1;
  
  mC( index2  , index2+3 ) += ee2;
  mC( index2+1, index2+2 ) -= ee2;

  real bb1ee2 = bb1*ee2;
  real ee1bb2 = ee1*bb2;
  real ee1aa2 = ee1*aa2;
  real aa1ee2 = aa1*ee2;
  
  mC( index1  , index2+1 ) +=  bb1ee2 - ee1bb2;
  mC( index1  , index2+3 ) -=  bb1ee2 + ee1aa2;

  mC( index1+1, index2   ) +=  ee1bb2 - bb1ee2;
  mC( index1+1, index2+2 ) +=  ee1aa2 + bb1ee2;

  mC( index1+2, index2+1 ) +=  aa1ee2 + ee1bb2;
  mC( index1+2, index2+3 ) -=  aa1ee2 - ee1aa2;

  mC( index1+3, index2   ) -=  ee1bb2 + aa1ee2;
  mC( index1+3, index2+2 ) += -ee1aa2 + aa1ee2;
#else
  MSG.error("Stifflink","only implemented in 2D");
#endif  
}

//==========================================================================
//update the matrix to include a force between point ii and position g
//the force is linear:  f_a = weight * ( g - a );
void interactionClamp( const PointExact & pta, const real g[], const real weight = 1.0 )
{
  int aa = pta.matIndex();
  
  mB( aa, aa ) -=  weight;

  if ( Space::isPeriodic() ) {
    
    real gm[DIM];
    Space::moduloNear( gm, g, &vPTS[ DIM*aa ] );
  
    for(int dd = 0; dd < DIM; ++dd )
      vBAS[ DIM*aa + dd ] += weight * gm[ dd ];
  
  } else {
    
    for(int dd = 0; dd < DIM; ++dd )
      vBAS[ DIM*aa + dd ] += weight * g[ dd ];
  
  }
}


//==========================================================================
//update the matrix to include a linear force of given weigth between
//a microtubule point si, and a position in space g[]
void interactionClamp( const PointInterpolated & pti, const real g[], real weight = 1.0 )
{
  assert( pti.looksWrong(true) == NO_ERROR );

  int index = pti.matIndex1();

  real ca = pti.getCoef();
  real cb = 1.0 - ca;

  mB( index,   index   ) -=  weight * cb * cb;
  mB( index,   index+1 ) -=  weight * ca * cb;
  mB( index+1, index+1 ) -=  weight * ca * ca;

  if ( Space::isPeriodic() ) {
    
    real gm[DIM];
    Space::moduloNear( gm, g, &vPTS[ DIM*index ] );
    
    for(int dd=0; dd < DIM; ++dd ) {
      vBAS[ DIM*index + dd       ] += weight * cb * gm[ dd ];
      vBAS[ DIM*index + dd + DIM ] += weight * ca * gm[ dd ];
    }
    
  } else {
    
    for(int dd=0; dd < DIM; ++dd ) {
      vBAS[ DIM*index + dd       ] += weight * cb * g[ dd ];
      vBAS[ DIM*index + dd + DIM ] += weight * ca * g[ dd ];
    }
    
  }
}


//==========================================================================
//--------------------------------------------------------------------------
// update the matrices to include an interaction with non-zero resting length
// between a microtubule at position a, and a fixed point g[] 
void interactionStiffClamp( const PointInterpolated & a, const real g[], const real len )
{
#if ( DIM == 2 )
  assert( a.looksWrong(true) == NO_ERROR );
  assert( Space::isPeriodic() == false );
  

  //force coefficients on the points:
  real aa = a.getCoef(); 
  real bb = 1 - aa;
  real ee, eeee;
  
  Vecteur ga = a.where();
  ga[0] -= g[0];
  ga[1] -= g[1];
  
  if (( ga ^ a.dirInter() ) > 0 ) {
    ee =  len / a.length();
  } else {
    ee = -len / a.length();
  }
  eeee = ee * ee;
    
  //index in the matrix mB:
  int index = a.matIndex1();   

  //we put the isotropic terms in mB
  mB( index,   index   ) -=  bb * bb + eeee;
  mB( index,   index+1 ) += -aa * bb + eeee;
  mB( index+1, index+1 ) -=  aa * aa + eeee;
  
  //index in the matrix mC:
  index *= DIM;
  
  //it seems to works also fine without the term in ee* below:
  vBAS[ index     ] += bb * g[0] - ee * g[1];
  vBAS[ index + 1 ] += bb * g[1] + ee * g[0];
  vBAS[ index + 2 ] += aa * g[0] + ee * g[1];
  vBAS[ index + 3 ] += aa * g[1] - ee * g[0];
  
  mC(index  , index+3) += ee;
  mC(index+1, index+2) -= ee;
#else
  MSG.error("Stifflink","only implemented in 2D");
#endif  
}


//==========================================================================
// update the matrix to include the force from a flat slipery boundary
// w is the current position of the point on which force is acting
// ii is the index of the associated variable in vPTS
// p = projection of w on the surface of the box

void interactionPlane(const PointExact & aa, const Vecteur w, const Vecteur p, real weight = 1.0 )
{
  Vecteur r = p - w;
  real rsq = r.normSquare();

  //if the vector is too small, the direction will be imprecise
  //so we prefer not to include the interaction:
  if ( rsq < 1e-6 ) return;

  rsq = weight / rsq;

  int ii = DIM * aa.matIndex();
  
  real pr = ( w * r ) * rsq + weight;

  for( int dd = 0; dd < DIM; ++dd )
    vBAS[ ii + dd ] += pr * r[ dd ];

  mC(ii  , ii   ) -= r.XX * r.XX * rsq;

#if ( DIM >= 2 )
  mC(ii  , ii+1 ) -= r.XX * r.YY * rsq;
  mC(ii+1, ii+1 ) -= r.YY * r.YY * rsq;
#endif

#if ( DIM >= 3 )
  mC(ii  , ii+2 ) -= r.XX * r.ZZ * rsq;
  mC(ii+1, ii+2 ) -= r.YY * r.ZZ * rsq;
  mC(ii+2, ii+2 ) -= r.ZZ * r.ZZ * rsq;
#endif
}


//--------------------------------------------------------------------------
// update the matrix to include a force for a cortical grafted
// ( experimental )
void interactionPlane( Grafted * gh )
{
  assert( gh->isAttached() );

  Vecteur q = gh->whereHand();
  Vecteur p = gh->whereGrafted();

  PointExact pta( gh->getHand().getPS(), gh->getHand().getPoint1() );
  interactionPlane( pta, q, p, 1 - gh->getHand().getCoef() );
  
  pta.setTo( gh->getHand().getPS(), gh->getHand().getPoint2() );
  interactionPlane( pta, q, p,     gh->getHand().getCoef() );
}


//==========================================================================
//==========================================================================

//update the matrix to add a force from the center k, of resting length len,
//and given weight: f_a = weight * ca * ( len / norm(ca) - 1 ); ca = a - center;
void interactionLongClamp( const PointExact & aa, const Vecteur center, real len, const real weight = 1.0 )
{
  int index = DIM * aa.matIndex(); 
  
  Vecteur ca;  
  for(int dd = 0; dd < DIM; ++dd )
    ca[ dd ] = vPTS[ index+dd ] - center[ dd ];
  
  real can = ca.norm();
  if ( can < 1e-6 ) return;
  ca /= can;
    
  if ( len < can ) {
  
    len /= can;
    
    real facX = weight * len * ( can + ca * center );
    real facC = weight * ( 1.0 - len );
  
    for(int ii = 0; ii < DIM; ++ii ) {

      mC( index + ii, index + ii ) += weight * ( len * ( 1.0 - ca[ii] * ca[ii] ) - 1.0 ); 
        
      for(int jj = ii+1; jj < DIM; ++jj )
        mC( index + ii, index + jj ) -= weight * len * ca[ii] * ca[jj];
      
      vBAS[ index + ii ] += facX * ca[ ii ] + facC * center[ ii ];
    }

  } else {
    
    real facX = weight * ( len + ca * center );

    for(int ii = 0; ii < DIM; ++ii ) {
      for(int jj = ii; jj < DIM; ++jj )
        mC( index + ii, index + jj ) -= weight * ca[ii] * ca[jj];
      
      vBAS[ index + ii ] += facX * ca[ ii ];
    }

  }
}

//==========================================================================
//                             TEST-FORCES
//==========================================================================

void interactionCoulomb( Vecteur ab, int aa, int bb, real c = 0.01 )
{
  int ii, jj;
  real m, abnsq = ab.normSquare(), abn=sqrt(abnsq);

  if ( abn < EPSILON ) return;
  ab /= abn;

  real abn3 = c / abnsq;
  real abn5 = c / ( abnsq * abn );

  for( ii = 0; ii < DIM; ++ii )
  {
    vBAS[ DIM*aa + ii ] += 3 * abn3 * ab[ii];
    vBAS[ DIM*bb + ii ] -= 3 * abn3 * ab[ii];
  }
  
  for( ii = 0; ii < DIM; ++ii )
    for( jj = ii; jj < DIM; ++jj )
    {
      m = abn5 * ((ii==jj) - 3 * ab[ii] * ab[jj] );
      
      mC( DIM*aa + ii, DIM*aa + jj ) += m;
      mC( DIM*bb + ii, DIM*bb + jj ) += m;

      mC( DIM*aa + ii, DIM*bb + jj ) -= m;
      if ( ii != jj )
        mC( DIM*aa + jj, DIM*bb + ii ) -= m;
    }
}


