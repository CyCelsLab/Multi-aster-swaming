//RCS: $Id: nucleus.cc,v 2.38 2005/04/22 12:28:29 nedelec Exp $
//RCS: $Id: nucleus.cc, v2.4c 2006/02/09 20:03:47 athale Exp $
//--------------------------------nucleus.cc-----------------------------------


#include "assert_macro.h"
#include "grafted.h"
#include "nucleus.h"
#include "cblas.h"
#include "clapack.h"
#include "iowrapper.h"
#include "exceptions.h"
#include "transformation.h"
#include "sim.h"


Space * Nucleus::space = 0;

#if  ( DIM == 1 )
const int Nucleus::nbrefpts = 0;
#elif( DIM == 2 )
const int Nucleus::nbrefpts = 1;
#elif( DIM == 3 )
const int Nucleus::nbrefpts = 3;
#endif


//============================================================================
//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
//============================================================================

//------------------- construction and destruction ---------------------------


void Nucleus::constructorBase()
  //this resets the pointers to arrays,
  //it should be called only by the constructors
{
  nuRadius     = MP.nuradius;

  // set the reference points
  #if( DIM > 1 )
  ref1.XX      = nuRadius;
  ref1.YY      = 0.;
  #endif
  #if( DIM > 2 )
  ref1.ZZ      = 0.;
  ref2.XX      = 0.;
  ref2.YY      = nuRadius;
  ref2.ZZ      = 0.;
  ref3.XX      = 0.;
  ref3.YY      = 0.;
  ref3.ZZ      = nuRadius;
  #endif

  psNoise      = 0;
  transNoise   = new real[DIM];
  #if   (DIM == 1)
  rotNoise     = 0;
  #elif (DIM == 2)
  rotNoise     = new real[DIM - 1];
  #elif (DIM == 3)
  rotNoise     = new real[DIM];
  #endif

  //nuT = MZERO;
  //nuR = MZERO;

  // for the projection
  nucJComponents = 0;
  nucPPrimeForce = 0;
  allocatedP     = 0;

  // for the microtubules:
  nuFocus        = MT_MINUS_END;
}


Nucleus::Nucleus() {

  constructorBase();

  //The default constructor is mainly used when a nucleus is read from a file.
  //It is called from nucleus_list.
  //The memory for the points is allocated in Pointset::read().
  //The memory for the clamps and the mt is allocated in Nucleus::read().

  //set the mobilities of the nucleus
  setMobility();

}


Nucleus::Nucleus(int dummy) {

  Vecteur w, x, y, nuV;
  int     ghtotal;    // total number of grafted hands

  constructorBase();

  //allocate memory:
  //Memory for the points is allocated in pointset, when the points are added
  //with addPoint() further down. Memory for the clamps and the mt is
  //allocated in Nucleus::addClamp...() and Nucleus::addMicrotub() functions.

  //fix the initial location from data.in
  if ( MP.nuinit == 22  )
    w = 	( space, MP.nuinit, MP.nuVect );
  else
    w = initPosition( space, MP.nuinit );
  //set it based on Object::initPosition code

  //initialise points
  ghtotal = 0;
  for(int ii = 0; ii < MP.MAX; ii++)
    ghtotal += MP.nughmax[ii];

  assert( MP.nuptsmax  >= 0 );
  assert( MP.numtmax   >= 0 );
  assert( ghtotal      >= 0 );

  //===========================
  //----- add free points -----
  //===========================

  //free points were mainly used for testing
  switch( MP.nuptsmax ) {
  case 0:
    addPoint( w );
    break;
  case 1:
    addPoint( w );

                x[0] = -nuRadius;
    if(DIM > 1) x[1] = 0.;
    if(DIM > 2) x[2] = 0.;
    addPoint( w + x );
    break;
  case 2:
    addPoint( w );

                x[0] = nuRadius;
    if(DIM > 1) x[1] = 0.;
    if(DIM > 2) x[2] = 0.;
    addPoint( w + x );

                x[0] = -nuRadius;
    if(DIM > 1) x[1] = 0.;
    if(DIM > 2) x[2] = 0.;
    addPoint( w + x );
    break;
  case 3:
    addPoint( w );

                x[0] = nuRadius;
    if(DIM > 1) x[1] = 0.;
    if(DIM > 2) x[2] = 0.;
    addPoint( w + x );

                x[0] = -nuRadius;
    if(DIM > 1) x[1] = 0.;
    if(DIM > 2) x[2] = 0.;
    addPoint( w + x );

                x[0] = 0.;
    if(DIM > 1) x[1] = -nuRadius;
    if(DIM > 2) x[2] = 0.;
    addPoint( w + x );
    break;
  default:
    addPoint( w );

                x[0] = nuRadius;
    if(DIM > 1) x[1] = 0.;
    if(DIM > 2) x[2] = 0.;
    addPoint( w + x );

    for(int i = 2; i <= MP.nuptsmax; i++) {
      x = nuRadius * Vecteur::randNormed();
      addPoint( w + x );
    }
    break;
  }

  //=============================
  //----- add grafted hands -----
  //=============================

  int nughpint[ParamSim::MAX];
  for(int ii = 0; ii < ParamSim::MAX; ii++)
    if(MP.nughmax[ii] == 0)
      nughpint[ii] = 0;
    else {
      // we use floor(x+0.5) instead of round since icc doesn't know round
      nughpint[ii] = (int)floor( (real)MP.nughmax[ii]/(real)ghtotal + 0.5 )*100;
    }

  if ( ghtotal > 0 )
    for(int i = 0; i < ghtotal; i++) {
      x = nuRadius * Vecteur::randNormed();
      addPointWithGrafted( w+x, RNG.pint_ratio(MP.MAX, nughpint) );
    }

  //=======================================
  //----- add MTOCs with microtubules -----
  //=======================================

  if ( MP.numtmax > 0 ) {
    // create points and clamps and add tubes to it
    for(int ii = 0; ii < MP.numtmax; ii++) {
      x = nuRadius * Vecteur::randNormed();
      addMicrotub( addPointAndClamp( w + x ));
    }
  }
  
  
/////CHAITANYA-- HERE WE CAN ADD OUR distance dependent nucleation
  //===================================
  //----- add microtubule bundles -----
  //===================================

  if ( MP.nubundles > 0 ) {
    if ( MP.numtperbundle >=2 ) {

      int            indx1  = -1;
      int            indx2  = -1;
      Microtub*      mt1    =  0;
      Microtub*      mt2    =  0;
      Vecteur        place;
      Vecteur        dir;
      Vecteur        axis;
      Vecteur        storage;
      Transformation rotMat;

      //impose antiparallel overlap of microtubules
      for(int ii=0; ii < MP.nubundles; ii++ ) {
        assert( MP.numtperbundle >= 2 );

        //connect two tubes per bundle to the nucleus
        //the position of the first tube is random on the sphere
        x     = nuRadius * Vecteur::randNormed();
        indx1 = addPointAndClamp( w + x );

        //try to choose the position of the second tube wisely
        //to avoid high stress on the tubes

        //the angle theta between x and y is chosen such that the two points
        //x and y are (MP.initsize[0] + MP.initsize[1]) apart from each other
        #if ( DIM != 1 )
        real theta = 2.*asin( 0.5*(MP.initsize[1] + MP.initsize[0])/nuRadius );
        #endif
        #if ( DIM == 2 )
          rotMat.setRotationFromEulerAngles( theta, 0 ,0 );
          y = rotMat*x;
        #elif ( DIM == 3 )
          //first rotate "theta" away from x

          //find a rotation axis in the plane perpendicular to x
          storage.set( 1, 0, 0 );
          axis = x^storage;
          if ( axis.norm() == 0 ) {
            //in case x was parallel to (1,0,0), try again with (0,1,0)
            storage.set( 0, 1, 0 );
            axis = x^storage;
          }
          axis.normalize();
          rotMat.setRotationAroundAxis( axis, theta );
          storage = rotMat*x;

          //with the second rotation choose a random point on the circle
          //with distance (MP.initsize[0] + MP.initsize[1]) around x
          real phi = RNG.real_uniform_range( -PI, PI );
          axis = x;
          axis.normalize();
          rotMat.setRotationAroundAxis( axis, phi );
          y    = rotMat*storage;
        #endif
        dir   = y - x;
        indx2 = addPointAndClamp( w + y );
        mt1   = addMicrotub( indx1, w + x,  dir );
        mt2   = addMicrotub( indx2, w + y, -dir );

        //now link the two mt
        nuClamp[ indx2 ].setMT( indx1 );
        nuClamp[ indx1 ].setBundle( 1 );
        nuClamp[ indx2 ].setBundle( 1 );

        //then link free mts to the bundle
        for(int kk=1; kk < MP.numtperbundle-1; kk++ ) {
          if ( kk % 2 ) {
            addMicrotub( addClampToMT( indx2 ), w+x,  dir );
          }
          else {
            addMicrotub( addClampToMT( indx1 ), w+y, -dir );
          }
        }
      }
      
    }
    else
      MSG.warning("Nucleus, nucleus.cc: ", "bundles of need at least 2 mt, no bundles created");
  }

  //always add reference points (0 in 1D, 1 in 2D, 3 in 3D)
  //these points are needed to display the rotation of the sphere correctly
  #if( DIM > 1 )
              x[0] = ref1.XX;
              x[1] = ref1.YY;
  if(DIM > 2) x[2] = ref1.ZZ;
  addPoint( w+x );
  #endif
  #if( DIM > 2 )
  x[0] = ref2.XX;
  x[1] = ref2.YY;
  x[2] = ref2.ZZ;
  addPoint( w+x );
  x[0] = ref3.XX;
  x[1] = ref3.YY;
  x[2] = ref3.ZZ;
  addPoint( w+x );
  #endif

  // set the mobility of the nucleus
  setMobility();

}


Nucleus::~Nucleus() {
  //free memory
  if ( transNoise     ) delete[] transNoise;
  if ( rotNoise       ) delete[] rotNoise;
  if ( psNoise        ) delete[] psNoise;
  if ( nucJComponents ) delete[] nucJComponents;
  if ( nucPPrimeForce ) delete[] nucPPrimeForce;
}


int Nucleus::addClampToMT( const int mtIndx )
{
  assert(( mtIndx >= 0 ) && ( mtIndx < maxNbMicrotub() ));

  (nuClamp.pushFront()).set( -1, mtIndx, 1 );
  return nuClamp.size()-1;
}


int Nucleus::addPointAndClamp( const Vecteur where )
{
  (nuClamp.pushFront()).set( addPoint( where ), -1, 0 );
  return nuClamp.size()-1;
}


void Nucleus::addPointWithGrafted( Vecteur w, int ty )
{
  linkGrafted( new Grafted( ty, this, addPoint( w ) ) );
}


Microtub* Nucleus::addMicrotub( const int clampIndex )
{
  Vecteur dir, place;
  Microtub* addedTube;

  place = whereP( nuClamp[clampIndex].clamp );
  dir   = whereP( nuClamp[clampIndex].clamp ) - whereP( 0 );

  if ( hasMicrotub(clampIndex) ) {
    addedTube =  getMicrotub(clampIndex);
    addedTube -> setStraight( place, dir, nuFocus );
  }
  else {
    addedTube = new Microtub( place, dir, Microtub::initLength(), nuFocus );
    registerMicrotub( addedTube, clampIndex );
  }

  return addedTube;
}


Microtub* Nucleus::addMicrotub( const int clampIndex, const Vecteur& place, const Vecteur& dir )
{
  Microtub* addedTube;

  if ( hasMicrotub(clampIndex) ) {
    addedTube =  getMicrotub(clampIndex);
    addedTube -> setStraight( place, dir, nuFocus );
  }
  else {
    addedTube = new Microtub( place, dir, Microtub::initLength(), nuFocus );
    registerMicrotub( addedTube, clampIndex );
  }

  return addedTube;
}


int  Nucleus::removePoint( const int ptIndex )
//remove one point from the nucleus
//we have to clean up the c-array pspts afterwards
{
  //the center point is at the beginning of the array
  //the reference points are at the end
  assert( nbPoints() >  nbrefpts + 1 );
  assert( ptIndex    <  nbPoints()-nbrefpts );
  assert( ptIndex    >= 1 );

  //move all points above ptindex one point down
  for( int ii=ptIndex; ii<nbPoints()-1; ii++ ) {
    pspts[ii*DIM]   = pspts[(ii+1)*DIM];
    #if ( DIM > 1 )
    pspts[ii*DIM+1] = pspts[(ii+1)*DIM+1];
    #endif
    #if ( DIM > 2 )
    pspts[ii*DIM+2] = pspts[(ii+1)*DIM+2];
    #endif
  }

  setNbPoints( nbPoints()-1 );
  return nbPoints();
}


//int Nucleus::removePointAndClamp( )
//{
//}


int Nucleus::deregisterMicrotub(Microtub * old_mt)
{
  MSG( 10," Nucleus::deregisterMicrotub (name: %lx)\n",old_mt->getName() );
  //we remove the Microtub which has disapeared, getting its index in moMT
  int indx = MicrotubOrganizer::deregisterMicrotub( old_mt );
  //printf("Index of deleted tube: %d\n",indx);
  assert(( indx >= 0 ) && ( indx < maxNbMicrotub() ));
  assert( hasMicrotub(indx) == false );

  //TODO: we could have a finite nucleation rate on nuclei.

  //here we immediately re-create a new Microtub on the same spot:
  if ( nuClamp[indx].clamp == -1 ) {
    //if the tube wasn't connected to a point on the nucleus, just delete it
    //nuClamp[indx].clear();
  }
  else {
    //if the tube was connected to an MTOC on the nucleus, create a new one
    Vecteur place    = whereP( nuClamp[indx].clamp );
    Vecteur dir      = old_mt->dirEnd( nuFocus );
    Microtub* new_mt = new Microtub( place, dir, MP.mtminlength, nuFocus );
    //register the new Tube in the spot where the old one was attached:
    registerMicrotub( new_mt, indx );
    //also link the MT into the simulation
    if ( isLinked() ) sim.link( new_mt );
  }
  return indx;
}


bool Nucleus::isDeletable( Microtub* mt ) const
{
  int indx = findIndex( mt );
  assert(( indx >=0 ) && ( indx < maxNbMicrotub() ));

  return ( nuClamp[indx].bundled == 0 );
}


void Nucleus::setMobility()
{
  /// \todo we have to find a Stokes-formula for the translational and the point mobility in 2D!
  nutransmobil = MP.numobil * 1 / ( 6 * PI * MP.visc * nuRadius );
  #if (DIM == 1)
  nurotmobil = 0.;
  #elif (DIM == 2)
  nurotmobil   = MP.numobil * 3 / ( 32 * MP.visc * nuRadius * nuRadius * nuRadius );
  #elif (DIM == 3)
  nurotmobil   = MP.numobil * 1 / ( 8 * PI * MP.visc * nuRadius * nuRadius * nuRadius );
  #endif
  psMobility      = 1 / ( 6 * PI * MP.nuvisc * nuRadius );
  //printf("Mobilities: nutransmobil: %f nurotmobil: %f psmobi: %f\n", nutransmobil, nurotmobil, psMobility);
}


//============================================================================
//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
//============================================================================

//------------------- space related functions --------------------------------


void Nucleus::moduloPosition()
{
  moduloPositionG();
  MicrotubOrganizer::moduloPosition();
}


void Nucleus::moduloPositionG()
// modulo around the center of mass, i.e. the first point
{
  Space::modulo( pspts );
  for(int p = 1; p < nbPoints(); p++ )
    Space::moduloNear( pspts+DIM*p, pspts );
}


int Nucleus::insidePosition( const Space * s ) const
{
  return MicrotubOrganizer::insidePosition(s) + PointSet::insidePosition(s);
}


void Nucleus::translatePosition( const Vecteur & T )
{
  PointSet::translatePosition( T );
  MicrotubOrganizer::translatePosition( T );
}


//============================================================================
//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
//============================================================================

//------------------- technical functions and mathematics --------------------


bool Nucleus::setPointClamp( PointExact * pta, PointExact * ptb, const int indx ) const
{
  assert(( indx >= 0 ) && ( indx < maxNbMicrotub() ));

  if ( hasMicrotub(indx) && ( nuClamp[indx].clamp != -1 ) ) {
    //this is a clamp from a point on the nucleus, and the nuFocus-end of the Microtub
    pta->setTo( this, nuClamp[indx].clamp );
    ptb->setTo( getMicrotub(indx), ( nuFocus == MT_MINUS_END )? 0 : getMicrotub(indx)->lastPoint() );
    return true;
  } else {
    return false;
  }
}


bool Nucleus::setMTClamp( PointInterpolated * pta, PointInterpolated * ptb, PointInterpolated * ptc, PointInterpolated * ptd, const int indx ) const
{
  assert(( indx >= 0 ) && ( indx < maxNbMicrotub() ));

  if ( hasMicrotub(indx) && ( nuClamp[indx].mtIndex != -1 ) ) {
    Microtub* clampTube = getMicrotub( nuClamp[indx].mtIndex );
    //create the clamps only, if the other tube does really exist
    if ( clampTube != 0 ) {
      //these are two clamps between two microtubules
      clampTube           -> setInterpolation( pta, MP.initsize[0], MT_MINUS_END );
      getMicrotub( indx ) -> setInterpolation( ptb, MP.initsize[1], MT_MINUS_END );
      clampTube           -> setInterpolation( ptc, MP.initsize[1], MT_MINUS_END );
      getMicrotub( indx ) -> setInterpolation( ptd, MP.initsize[0], MT_MINUS_END );
      return true;
    }
  }
  return false;
}


void Nucleus::multMobilityVect( real* V, real s ) const
  // V <- s*V*mu
{
  MSG.error("multMobilityVect", "Reached multMobilityVect in nucleus.cc! This should not be used anymore!");
  exit(0);
}


real Nucleus::addBrownianForces( real * force ) const {
  return -1;
}


real Nucleus::addBrownianMotion( real * rhs ) const
{
  // This function is used as an alternative to addBrownianForces.
  // For the nucleus this is necessary because we have different mobilities
  // and different types of Brownian motion for the points in the pointset.
  // Also the projections are different for each part of the motion.
  // If this is used, addBrownianForces should do nothing but simply return -1!

  real transbrownian;  // translational displacement of the rigid body
  real rotbrownian;    // rotational displacement of the rigid body
  real psbrownian;     // brownian force on the MTOCs

  transbrownian = sqrt( 2. * MP.kT * MP.dt * nutransmobil);
  rotbrownian   = sqrt( 2. * MP.kT * MP.dt * nurotmobil);
  psbrownian    = sqrt( 2. * MP.kT /(MP.dt * psMobility));

  // calculate the displacements due to the Brownian motion of the rigid body
  for(int jj = 0; jj < DIM; jj++ ) {
    transNoise[jj] = transbrownian * RNG.gauss();
  }

  #if   (DIM == 2)
  *rotNoise      = rotbrownian * RNG.gauss();
  #elif (DIM == 3)
  for(int jj = 0; jj < DIM; jj++ ) {
    rotNoise[jj] = rotbrownian * RNG.gauss();
  }
  #endif

  // calculate the Brownian Forces on the MTOCs
  // we don't want MTOC-diffusion on the center point and on the reference points!
  for(int jj = 0; jj < DIM; jj++ )
    psNoise[jj] = 0.;
  for(int jj = DIM; jj < DIM*(nbPoints()-nbrefpts); jj++ )
    psNoise[jj] = psbrownian * RNG.gauss();


  // add the displacements to vRHS

  bool noiseon = true;
  if ( noiseon ) {

  if(rigidBodyTrans) {
    // translational noise on the rigid body:
    for(int i = 0; i < nbPoints(); i++) {
      blas_xaxpy(DIM,1.0,transNoise,1,rhs+i*DIM,1);
    }
  }
  if(rigidBodyRot) {
    // rotational noise on the rigid body:
    #if (DIM == 2)
    real cosphi = cos(rotNoise[0]);
    real sinphi = sin(rotNoise[0]);
    real xx,yy;
    for(int i = 1; i < nbPoints(); i++) {
      xx = pspts[i*DIM]   - pspts[0];
      yy = pspts[i*DIM+1] - pspts[1];
      rhs[i*DIM]   += (cosphi*xx - sinphi*yy) - xx;
      rhs[i*DIM+1] += (sinphi*xx + cosphi*yy) - yy;
    }
    #endif
    #if (DIM == 3)
    real cosphi   = cos(rotNoise[0]);
    real sinphi   = sin(rotNoise[0]);
    real costheta = cos(rotNoise[1]);
    real sintheta = sin(rotNoise[1]);
    real cosgamma = cos(rotNoise[2]);
    real singamma = sin(rotNoise[2]);
    real xx,yy,zz;
    for(int i = 1; i < nbPoints(); i++) {
      xx = pspts[i*DIM]   - pspts[0];
      yy = pspts[i*DIM+1] - pspts[1];
      zz = pspts[i*DIM+2] - pspts[2];
      rhs[i*DIM]   += (cosphi*costheta                      *xx - sinphi         *yy - cosphi*sintheta                               *zz) - xx;
      rhs[i*DIM+1] += ((sinphi*costheta - sintheta*singamma)*xx + cosphi*cosgamma*yy - (sinphi*sintheta*cosgamma + costheta*singamma)*zz) - yy;
      rhs[i*DIM+2] += ((sinphi*costheta + sintheta*cosgamma)*xx + cosphi*singamma*yy - (sinphi*sintheta*singamma - costheta*cosgamma)*zz) - zz;
    }
    #endif
  }
  // the noise on the MTOCs
  if(MTOCMotion)
    setMTOCSpeedsFromForces(psNoise, rhs, MP.dt);

  //we estimate the contribution of Brownian motions to vRHS.
  //assuming the projection is ~1, we have vRHS ~ MP.dt * psMobility * FBrownian:
  /** \todo the rotational part is not considered in the estimation of the
            contribution of thermal flutctuations to vRHS! */
  }

  //printf("rot %f trans %f mtoc %f\n", rotbrownian, transbrownian, psbrownian*psMobility);
  return (transbrownian < MP.dt*psMobility*psbrownian) ? transbrownian : MP.dt*psMobility*psbrownian;
}


void Nucleus::orthogonalizeRef()
{
  #if( DIM < 3 )
  return;
  #else
  real sc;

  // load the points from the array into vectors
  Vecteur tmpref1 = pspts + DIM*(nbPoints()-3);
  Vecteur tmpref2 = pspts + DIM*(nbPoints()-2);
  Vecteur tmpref3 = pspts + DIM*(nbPoints()-1);
  pscenter = pspts;

  // reduce to the center of mass an normalize
  tmpref1 -= pscenter;
  tmpref2 -= pscenter;
  tmpref3 -= pscenter;
  tmpref1.normalize();
  tmpref2.normalize();
  tmpref3.normalize();

  // make ref2 orthogonal to ref1
  sc       = tmpref1*tmpref2;
  //printf("ref1*ref2: %f\n", sc);
  tmpref2 -= sc*tmpref1;
  tmpref2.normalize();

  // make ref3 orthogonal to ref1
  sc       = tmpref1*tmpref3;
  //printf("ref1*ref3: %f\n", sc);
  tmpref3 -= sc*tmpref1;
  tmpref3.normalize();

  // make ref3 orthogonal to ref2
  sc       = tmpref2*tmpref3;
  //printf("ref2*ref3: %f\n", sc);
  tmpref3 -= sc*tmpref2;
  tmpref3.normalize();

  tmpref1 *= nuRadius;
  tmpref2 *= nuRadius;
  tmpref3 *= nuRadius;
  tmpref1 += pscenter;
  tmpref2 += pscenter;
  tmpref3 += pscenter;

  // put the corrected vectors back into the array
  for(int ii=0; ii < 3; ii++) {
    pspts[DIM*( nbPoints()-3 ) + ii] = tmpref1[ii];
    pspts[DIM*( nbPoints()-2 ) + ii] = tmpref2[ii];
    pspts[DIM*( nbPoints()-1 ) + ii] = tmpref3[ii];
  }
  #endif
}


void Nucleus::reshapeWithConservedCenterOfGravity()
  // get rid of finite-step errors but conserve the shape
  // and the center of gravity for the nucleus
{
  int nbMTOCs            = nbSegments();
  static int   allocated = 0;
  static real* E         = 0;
  static real* axis      = 0;

  if(allocated < DIM) {
    if( E )    delete[] E;
    if( axis ) delete[] axis;
    allocated = DIM;
    E    = new real[allocated];
    axis = new real[allocated];
  }

  for(int i = 0; i < DIM; i++) {
    E[i] = 0.0;
  }
  for(int j = 1; j <= nbMTOCs; j++) {
    for(int i = 0; i < DIM; i++) {
      axis[i] = pspts[i+j*DIM] - pspts[i];
    }
    real norm = blas_xnrm2(DIM,axis,1);
    for(int i = 0; i < DIM; i++) {
      axis[i]        /= norm;
      E[i]           -= (norm - nuRadius)/(real)(nbMTOCs + 1)*axis[i];
      pspts[i+j*DIM]  = pspts[i] + axis[i]*nuRadius;
    }
  }
  for(int j = 0; j <= nbMTOCs; j++) {
    for(int i = 0; i < DIM; i++) {
      pspts[i+j*DIM] += E[i];
    }
  }

  orthogonalizeRef();
}


void Nucleus::reshapeWithConservedCenter()
  // this is used only for testing purposes!
  // we get rid of finite-step errors but conserve the shape
  // by projecting back onto the sphere,
  // without changing the postition of the center point
{
  int nbMTOCs            = nbSegments();
  static int   allocated = 0;
  static real* axis      = 0;

  if(allocated < DIM) {
    if( axis ) delete[] axis;
    allocated = DIM;
    axis = new real[allocated];
  }

  for(int j = 1; j <= nbMTOCs; j++) {
    for(int i = 0; i < DIM; i++) {
      axis[i] = pspts[i+j*DIM] - pspts[i];
    }
    real norm = blas_xnrm2(DIM,axis,1);
    for(int i = 0; i < DIM; i++) {
      axis[i]        /= norm;
      pspts[i+j*DIM]  = pspts[i] + axis[i]*nuRadius;
    }
  }

  orthogonalizeRef();
}


//=============================================================================
//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
//=============================================================================

//------------------- methods for the projection ------------------------------


void Nucleus::allocateProjection()
{
  int nbMTOCs = nbSegments();
  //int vecDim  = DIM*nbPoints();

  if (allocatedP < nbMTOCs) {
    //printf("Allocating more projection memory!\n");

    if( nucJComponents ) delete[] nucJComponents;
    if( psNoise )        delete[] psNoise;
    if( nucPPrimeForce ) delete[] nucPPrimeForce;

    allocatedP     = (nbMTOCs+NUCCHUNK);
    // we don't calculate MTOC motion for the reference points
    nucJComponents = new real[(allocatedP-nbrefpts)*DIM];
    // we don't want MTOC noise on the reference points:
    psNoise        = new real[(allocatedP+nbrefpts-1)*DIM];
    nucPPrimeForce = new real[(allocatedP+1)*DIM*(allocatedP+1)*DIM];
  }
}


void Nucleus::setProjectedForces(const real * X, real * Y) const
{
  //Projection for the nucleus consists of three parts:
  //rigid body translation, rigid body rotation and the motion of the points
  //on the sphere. Furthermore we have to multiply each part with different
  //mobilities. setSpeedsFromForces was modified to do this.

  MSG.error("Nucleus::setProjectedForces", "setProjectedForces can not be used in nucleus!\nUse setSpeedsFromForces instead!");
  exit(0);
}


#if (DIM == 1)

//this is unsafe, don't use the nucleus in 1D!
void Nucleus::prepareProjection() {}
void Nucleus::setRotSpeedsFromForces(real const*, real*, const real) const {}
void Nucleus::setMTOCSpeedsFromForces(real const*, real*, const real) const {}
void Nucleus::setTransSpeedsFromForces(real const*, real*, const real) const {}
void Nucleus::setSpeedsFromForces(const real* X, real* Y, const real sc) const
{
  for( int ii = 0; ii < DIM*nbPoints(); ii++ )
    Y[ii] = 0.0;
}

#elif ( DIM == 2 || DIM == 3 )

void Nucleus::prepareProjection()
//prepare variables for the projection
{
  //allocate more memory if needed
  allocateProjection();

  //preparations for the rigid body motion:
  if(rigidBodyTrans || rigidBodyRot) {
    // set center of mass:
    pscenter=pspts;
  }

  //preparations for the motion of the MTOCs:
  //the reference points will be omitted!
  if(MTOCMotion) {
    int nbMTOCs = nbSegments() - nbrefpts;

    //calculate the nonzero components of j, omitting factors of 2
    for(int i = 0; i < nbMTOCs; i++) {
      for(int k = 0; k < DIM; k++) {
        nucJComponents[i*DIM + k] = (pspts[(i+1)*DIM + k] - pspts[k]);
      }
    }
  }
}


void Nucleus::setTransSpeedsFromForces(const real* X, real* Y, const real sc) const
{
  Vecteur      F = VZERO;    // resulting total force

  //forces are not allowed on the reference points
  for( int p = 0; p < nbPoints()-nbrefpts; p++ )
    F += X+p*DIM;

  F *= sc*nutransmobil;

  for( int p = 0; p < nbPoints(); p++ ) {
    #if   ( DIM == 2 )
    Y[ p*DIM ]   = F.XX;
    Y[ p*DIM+1 ] = F.YY;
    #elif ( DIM == 3 )
    Y[ p*DIM   ] = F.XX;
    Y[ p*DIM+1 ] = F.YY;
    Y[ p*DIM+2 ] = F.ZZ;
    #endif
  }
}


void Nucleus::setRotSpeedsFromForces(const real* X, real* Y, const real sc) const
{
  Vecteur      F = VZERO;    // resulting total force
  CrossProduct T = MZERO;

  // forces are not allowed on the reference points
  for( int p = 0; p < nbPoints()-nbrefpts; p++ ) {
    F += X+p*DIM;
    #if   ( DIM == 2 )
    T    += pspts[ p*DIM ]   * X[ p*DIM+1 ] - pspts[ p*DIM+1 ] * X[ p*DIM ];
    #elif ( DIM == 3 )
    T.XX += pspts[ p*DIM+1 ] * X[ p*DIM+2 ] - pspts[ p*DIM+2 ] * X[ p*DIM+1 ];
    T.YY += pspts[ p*DIM+2 ] * X[ p*DIM   ] - pspts[ p*DIM   ] * X[ p*DIM+2 ];
    T.ZZ += pspts[ p*DIM   ] * X[ p*DIM+1 ] - pspts[ p*DIM+1 ] * X[ p*DIM   ];
    #endif
  }

  T -= pscenter ^ F;         // reduce the torque to the center of mass
  T *= sc*nurotmobil;        // multiply by the mobility and maybe dt
  F = -T^pscenter;

  for( int p = 0; p < nbPoints(); p++ ) {
    #if   ( DIM == 2 )
    Y[ p*DIM ]   += F.XX - T*pspts[ p*DIM + 1 ];
    Y[ p*DIM+1 ] += F.YY + T*pspts[ p*DIM ];
    #elif ( DIM == 3 )
    Y[ p*DIM   ] += F.XX + T.YY * pspts[ p*DIM+2 ] - T.ZZ * pspts[ p*DIM+1 ];
    Y[ p*DIM+1 ] += F.YY + T.ZZ * pspts[ p*DIM   ] - T.XX * pspts[ p*DIM+2 ];
    Y[ p*DIM+2 ] += F.ZZ + T.XX * pspts[ p*DIM+1 ] - T.YY * pspts[ p*DIM   ];
    #endif
  }
}


void Nucleus::setMTOCSpeedsFromForces(const real* X, real* Y, const real sc) const
{
  // no MTOC-motion on the reference points
  int  nbMTOCs = nbSegments() - nbrefpts;
  real JtJX;

  // calculate JtJ*X
  for(int i = 0; i < nbMTOCs; i++) {
    JtJX = 0.;
    for(int k = 0; k < DIM; k ++ ) {
      JtJX += nucJComponents[i*DIM + k] * X[(i+1)*DIM + k];
    }
    JtJX /= nuRadius*nuRadius;
    for(int k = 0; k < DIM; k ++ ) {
      Y[(i+1)*DIM + k] += sc*psMobility*(X[(i+1)*DIM + k] - nucJComponents[i*DIM + k]*JtJX);
    }
  }
}


void Nucleus::setSpeedsFromForces(const real* X, real* Y, const real sc) const
{
  // set contribution from the translation of the rigid body
  // or set Y to zero
  if(rigidBodyTrans)
    setTransSpeedsFromForces( X, Y, sc );
  else
    blas_xzero(DIM*nbPoints(),Y,1);

  // add contribution from the rotation of the rigid body
  if(rigidBodyRot)
    setRotSpeedsFromForces( X, Y, sc );

  // add contribution from the MTOC motion
  if(MTOCMotion)
    setMTOCSpeedsFromForces( X, Y, sc );
}
#endif


//============================================================================
//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
//============================================================================

//------------------- methods for projection correction ----------------------


void Nucleus::prepareProjectionDiff( const real *force )
{
  //int nbMTOCs = nbSegments();
  int  vecDim  = nbPoints()*DIM;

  if(rigidBodyTrans) {
    // nothing to do here.
    // translations seem to be pretty stable
  }

  if(rigidBodyRot) {
    // initialize the matrix nucPPrimeForce to zero
    for(int ii = 0; ii < vecDim; ii++)
      for(int kk = 0; kk < vecDim; kk++)
        nucPPrimeForce[ ii*vecDim+kk ] = 0.;

    #if ( DIM == 2)
    // calculate PPrime * Force from scratch
    for(int kk = 1; kk < nbPoints(); kk++) {
      for(int ii = 1; ii < nbPoints(); ii++) {
        //nucPPrimeForce[  kk*DIM   *vecDim + ii*DIM   ] +=  0.*force[ 0 ];
        nucPPrimeForce[  kk*DIM   *vecDim + ii*DIM   ] += -((pspts[ ii*DIM+1 ] - pspts[ 1 ])+(pspts[ kk*DIM+1 ] - pspts[ 1 ]))*force[ 1 ];

        nucPPrimeForce[ (kk*DIM+1)*vecDim + ii*DIM   ] +=  (pspts[ ii*DIM+1 ] - pspts[ 1 ])*force[ 0 ];
        nucPPrimeForce[ (kk*DIM+1)*vecDim + ii*DIM   ] +=  (pspts[ kk*DIM   ] - pspts[ 0 ])*force[ 1 ];

        nucPPrimeForce[  kk*DIM   *vecDim + ii*DIM+1 ] +=  (pspts[ kk*DIM+1 ] - pspts[ 1 ])*force[ 0 ];
        nucPPrimeForce[  kk*DIM   *vecDim + ii*DIM+1 ] +=  (pspts[ ii*DIM   ] - pspts[ 0 ])*force[ 1 ];

        nucPPrimeForce[ (kk*DIM+1)*vecDim + ii*DIM+1 ] += -((pspts[ ii*DIM ] - pspts[ 0 ])+(pspts[ kk*DIM ] - pspts[ 0 ]))*force[ 0 ];
        //nucPPrimeForce[ (kk*DIM+1)*vecDim + ii*DIM+1 ] +=  0.*force[ 1 ];
       for(int jj = 1; jj < nbPoints(); jj++) {
          //nucPPrimeForce[  kk*DIM   *vecDim + ii*DIM   ] += 0.*force[ jj*DIM ];
          if( ii == jj )
            nucPPrimeForce[  kk*DIM   *vecDim + ii*DIM   ] +=  (pspts[ kk*DIM+1 ] - pspts[ 1 ])*force[ jj*DIM+1 ];
          if( kk == jj )
            nucPPrimeForce[  kk*DIM   *vecDim + ii*DIM   ] +=  (pspts[ ii*DIM+1 ] - pspts[ 1 ])*force[ jj*DIM+1 ];

          if( kk == jj )
            nucPPrimeForce[ (kk*DIM+1)*vecDim + ii*DIM   ] += -(pspts[ ii*DIM+1 ] - pspts[ 1 ])*force[ jj*DIM ];
          if( ii == jj )
            nucPPrimeForce[ (kk*DIM+1)*vecDim + ii*DIM   ] += -(pspts[ kk*DIM   ] - pspts[ 0 ])*force[ jj*DIM+1 ];

          if( ii == jj )
            nucPPrimeForce[  kk*DIM   *vecDim + ii*DIM+1 ] += -(pspts[ kk*DIM+1 ] - pspts[ 1 ])*force[ jj*DIM ];
          if( kk == jj )
            nucPPrimeForce[  kk*DIM   *vecDim + ii*DIM+1 ] += -(pspts[ ii*DIM   ] - pspts[ 0 ])*force[ jj*DIM+1 ];

          if( ii == jj )
            nucPPrimeForce[ (kk*DIM+1)*vecDim + ii*DIM+1 ] +=  (pspts[ kk*DIM   ] - pspts[ 0 ])*force[ jj*DIM ];
          if( kk == jj )
            nucPPrimeForce[ (kk*DIM+1)*vecDim + ii*DIM+1 ] +=  (pspts[ ii*DIM   ] - pspts[ 0 ])*force[ jj*DIM ];
          //nucPPrimeForce[ (kk*DIM+1)*vecDim + ii*DIM+1 ] += 0.*force[ jj*DIM+1 ];
        }
      }
    }
    #elif ( DIM == 3 )
    //MSG.warning("prepareProjectionDiff", "the correction has not yet been implemented for rotation in 3D!");
    #endif
  }

  if(MTOCMotion) {
    /// \todo prepareProjectionDiff has to be rewritten for the MTOC motion
  }
}


void Nucleus::addProjectionDiff( const real *X, real *Y ) const
{
  //Because of the different mobilities we can't use this function.
  //We have to multiply by the different mobilities before we add the
  //correction. This is done in setSpeedsWithCorrection.

  MSG.error("addProjectionDiff", "addProjectionDiff can not be used in nucleus!\nUse setSpeedsWithCorrection instead!");
  exit(0);
}


void Nucleus::addRotProjectionDiff( const real *X, real *Y )
{
  int vecDim = DIM*nbPoints();

  blas_xgemv('n',vecDim,vecDim,-1.,nucPPrimeForce,vecDim,X,1,1.,Y,1);
}


void Nucleus::addMTOCProjectionDiff( const real *X, real *Y )
{
  // \todo addProjectionDiff has to be rewritten for the MTOC motion
  //blas_xcopy(DIM*nbPoints(), X, 1, Y, 1);
}


void Nucleus::setSpeedsWithCorrection(const real * X, real * Z, real * Y, const real sc) const
{
  //Adding the correction terms is not working up to now, so everything has
  //been disabled below. But anyway, we don't experience instabilities in
  //the nucleus, unless we torture it artificially with the mouse.

  assert( X != Y );
  assert( X != Z );
  assert( Y != Z );

  //----- we allocate an array backup of size DIM*nbPoints at least:
  static real *backup1   = 0;
  static real *backup2   = 0;
  static int   allocated = 0;
  if ( nbPoints() > allocated ) {
    if ( backup1 ) delete[] backup1;
    if ( backup2 ) delete[] backup2;
    allocated = nbPoints() + 2;
    backup1   = new real[ DIM*allocated ];
    backup2   = new real[ DIM*allocated ];
  }

  blas_xcopy(DIM*nbPoints(),Z,1,backup1,1);
  blas_xcopy(DIM*nbPoints(),Z,1,backup2,1);

  // set contribution from the translation of the rigid body
  // or set Y to zero
  if(rigidBodyTrans)
    // there's nothing to add here since translations seem to be pretty stable
    setTransSpeedsFromForces( Z, Y, sc );
  else
    blas_xzero(DIM*nbPoints(),Y,1);

  // add contribution from the rotation of the rigid body
  if(rigidBodyRot) {
    //addRotProjectionDiff( X, backup1 );
    setRotSpeedsFromForces( backup1, Y, sc );
  }

  // add contribution from the MTOC motion
  if(MTOCMotion) {
    //addMTOCProjectionDiff( X, backup2 );
    setMTOCSpeedsFromForces( backup2, Y, sc );
  }
}


void Nucleus::addSpeedDiffFromForces( const real* X, real* Y, real* FORCE, const real sc ) const
{
  //

  if(rigidBodyRot) {
    Vecteur      F = VZERO;    // resulting total force
    CrossProduct T = MZERO;    // resulting total torque

    for( int p = 1; p < nbPoints()-1; p++ ) {
      F += FORCE+p*DIM;
      #if   ( DIM == 2 )
      T    += pspts[ p*DIM ]   * FORCE[ p*DIM+1 ] - pspts[ p*DIM+1 ] * FORCE[ p*DIM ];
      #elif ( DIM == 3 )
      T.XX += pspts[ p*DIM+1 ] * FORCE[ p*DIM+2 ] - pspts[ p*DIM+2 ] * FORCE[ p*DIM+1 ];
      T.YY += pspts[ p*DIM+2 ] * FORCE[ p*DIM   ] - pspts[ p*DIM   ] * FORCE[ p*DIM+2 ];
      T.ZZ += pspts[ p*DIM   ] * FORCE[ p*DIM+1 ] - pspts[ p*DIM+1 ] * FORCE[ p*DIM   ];
      #endif
    }

    T -= pscenter ^ F;         // reduce the torque to the center of mass
//    T *= sc*nurotmobil;        // multiply by the mobility and maybe dt
//    F *= sc*nutransmobil;
//    F -= T^pscenter;

    for( int p = 0; p < nbPoints(); p++ ) {
      #if   ( DIM == 2 )
      real sccheat = -sc*fabs(pspts[ p*DIM ])/pspts[ p*DIM ]*fabs(pspts[ p*DIM + 1 ])/pspts[ p*DIM + 1 ];
      sccheat      = sc;
      Y[ p*DIM ]   += sccheat*(   +F[ 1 ]*(pspts[ p*DIM + 1 ] - pspts[ 1 ])) * X[ 0 ];
      Y[ p*DIM ]   += sccheat*(+T -F[ 0 ]*(pspts[ p*DIM + 1 ] - pspts[ 1 ])) * X[ 1 ];
      //printf("dcenter_x: %e\n", sc*(+T -F[ 0 ]*(pspts[ p*DIM + 1 ] - pspts[ 1 ])) * X[ 1 ]);

      Y[ p*DIM+1 ] += sccheat*(-T -F[ 1 ]*(pspts[ p*DIM ]     - pspts[ 0 ])) * X[ 0 ];
      Y[ p*DIM+1 ] += sccheat*(   +F[ 0 ]*(pspts[ p*DIM ]     - pspts[ 0 ])) * X[ 1 ];
      #endif
    }
    for( int p = 1; p < nbPoints(); p++ ) {
      #if   ( DIM == 2 )
      real deltax, deltay;
      for( int ii = 1; ii < nbPoints(); ii++ ) {
        real sccheat = -sc*fabs(pspts[ p*DIM ])/pspts[ p*DIM ]*fabs(pspts[ p*DIM + 1 ])/pspts[ p*DIM + 1 ];
        sccheat      = sc;
        deltax = sccheat*(-FORCE[ ii+1 ]*(pspts[ p*DIM + 1 ] - pspts[ 1 ])) * X[ DIM*ii ]
               + sccheat*(+FORCE[ ii   ]*(pspts[ p*DIM + 1 ] - pspts[ 1 ])) * X[ DIM*ii+1 ];
        if (p == ii)
          deltax += sccheat*(-T)*X[ DIM*ii+1 ];
        printf("deltax: %f\n", deltax);
	Y[ p*DIM+1 ] += deltax;

        //Y[ p*DIM ]     += sc*(-FORCE[ ii+1 ]*(pspts[ p*DIM + 1 ] - pspts[ 1 ])) * X[ DIM*ii ];
        //Y[ p*DIM ]     += sc*(+FORCE[ ii   ]*(pspts[ p*DIM + 1 ] - pspts[ 1 ])) * X[ DIM*ii+1 ];
        //if (p == ii)
        //  Y[ p*DIM ]   += sc*(-T)*X[ DIM*ii+1 ];

        deltay = sccheat*(+FORCE[ ii+1 ]*(pspts[ p*DIM ]     - pspts[ 0 ])) * X[ DIM*ii ]
               + sccheat*(-FORCE[ ii ]  *(pspts[ p*DIM ]     - pspts[ 0 ])) * X[ DIM*ii+1 ];
        if (p == ii)
          deltay += sccheat*(+T)*X[ DIM*ii ];
        printf("deltay: %f\n", deltay);
	Y[ p*DIM+1 ] += deltay;

        //Y[ p*DIM+1 ]   += sc*(+FORCE[ ii+1 ]*(pspts[ p*DIM ]     - pspts[ 0 ])) * X[ DIM*ii ];
        //Y[ p*DIM+1 ]   += sc*(-FORCE[ ii ]  *(pspts[ p*DIM ]     - pspts[ 0 ])) * X[ DIM*ii+1 ];
        //if (p == ii)
        //  Y[ p*DIM+1 ] += sc*(+T)*X[ DIM*ii ];
      }
      #elif ( DIM == 3 )
      //Y[ p*DIM   ] = F.XX + T.YY * pspts[ p*DIM+2 ] - T.ZZ * pspts[ p*DIM+1 ];
      //Y[ p*DIM+1 ] = F.YY + T.ZZ * pspts[ p*DIM   ] - T.XX * pspts[ p*DIM+2 ];
      //Y[ p*DIM+2 ] = F.ZZ + T.XX * pspts[ p*DIM+1 ] - T.YY * pspts[ p*DIM   ];
      #endif
    }
  }
  //else {
  //  for( int p = 0; p < DIM*nbPoints(); p++ ) {
  //    Y[p] = 0.;
  //  }
  //}
}


//============================================================================
//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
//============================================================================

//------------------- reading and writing the nucleus ------------------------


void Nucleus::write()
  // write the nucleus to a file:
{
  assert( name > 0 );
  IO.writeRecordTag("nu");
  IO.writeUInt16( name );
  PointSet::write();

  IO.writeFormattingNewLine();
  IO.writeUInt16( maxNbMicrotub() );

  for( int ii = 0; ii < maxNbMicrotub(); ii++ ) {
    Microtub * mt = getMicrotub(ii);
    IO.writeFormattingNewLine();
    IO.writeUInt16( mt ? mt -> getName() : 0 );
    IO.writeInt16( nuClamp[ ii ].clamp, nuClamp[ ii ].mtIndex, nuClamp[ ii ].bundled );
  }
}


void Nucleus::read()
  // read the nucleus from a file:
{
  setFlag( sim.frameInBuffer() );
  deregisterAllMicrotubs();

  try {
    PointSet::read();

    //read the number of attached microtubs:
    int nuMTmax = IO.readUInt16();
    nuClamp.setSize( nuMTmax );

    if ( IO.getFileFormat() > 19 )
      for( int ii = 0; ii < nuMTmax; ii++ ) {
        registerMicrotub( sim.readMicrotubName(), ii);
        nuClamp[ ii ].setClamp(  IO.readInt16() );
        nuClamp[ ii ].setMT(     IO.readInt16() );
        nuClamp[ ii ].setBundle( IO.readInt16() );
      }
    else
      for( int ii = 0; ii < nuMTmax; ii++ ) {
        registerMicrotub( sim.readMicrotubName(), ii);
        nuClamp[ ii ].setClamp( IO.readUInt16() );
        nuClamp[ ii ].setMT( -1 );
        nuClamp[ ii ].setBundle( 0 );
      }

  } catch( IOException e ) {
    //MSG("S%lx %i ", name, r);
    e.addBeforeMessage("Nucleus::read : ");
    clearPoints();
    throw e;
  }
}
