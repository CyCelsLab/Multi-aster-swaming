//RCS: $Id: nucleus.h,v 2.27 2005/04/10 14:51:55 nedelec Exp $
//-----------------------------------nucleus.h---------------------------------

#ifndef NUCLEUS_H
#define NUCLEUS_H

#include "space.h"
#include "simobject.h"
#include "microtub.h"
#include "microtub_organizer.h"


///Nucleus provides a spherical object with points that can move on it. Microtubules and Grafteds can be clamped to the points.
class Nucleus : public SimObject, public MicrotubOrganizer
{
 
private:
   

  //---------------------------------------------------------------------------
  //========================== private attributes =============================
  //---------------------------------------------------------------------------


  ///accessory class that stores information about the links between Nucleus and Microtubs
  struct              Clamp;

  ///the radius of the nucleus
  real                nuRadius;

  ///size of allocated projection memory
  int                 allocatedP;
  
  ///we always allocate a little bit more projection memory
  static const int    NUCCHUNK = 2;
  
  ///the matrix J used for projection (only the needed components)
  real*               nucJComponents;

  ///the matrix PPrime*Force used for projection
  real*               nucPPrimeForce;
    
  ///list of pointers to the attached grafteds
  Array<Grafted*>     nugrafteds;

  ///which end of the microtubules is attached to the nucleus: MT_PLUS_END or MT_MINUS_END
  MTEnd               nuFocus;

  ///list of clamps that link a microtubule to the nucleus or to other microtubules
  Array<Clamp>        nuClamp;

  ///torque on the nucleus (F: not used??)
  //CrossProduct        nuT;

  ///total rotation of the nucleus (F: not used??)
  //CrossProduct        nuR;
  
  ///the translational mobility of the nucleus
  real                nutransmobil;

  ///the rotational mobility of the nucleus
  real                nurotmobil;
  /* Note: (the mobility of the points on the nucleus is
     stored in psMobility in pointset) */

  ///the translational Noise on the nucleus
  real*               transNoise;

  ///the rotational Noise on the nucleus
  real*               rotNoise;

  ///the Noise of the points on the nucleus
  real*               psNoise;

  ///switch for testing: turn the rigid body translation on or off
  static const bool   rigidBodyTrans = 1;
  
  ///switch for testing: turn the rigid body rotation on or off
  static const bool   rigidBodyRot   = 1;
    
  ///switch for testing: turn the MTOC motion on or off
  static const bool   MTOCMotion     = 1;


  //---------------------------------------------------------------------------
  //=========================== private methods ==============================
  //--------------------------------------------------------------------------

  //------------------- construction and destruction -------------------------
  
  ///basic constructor, initializes variables etc.
  void      constructorBase();

  ///add a clamp that attaches one mt to another mt with index mtIndex
  int       addClampToMT( const int mtIndex );

  ///add a point with a clamp at place where
  int       addPointAndClamp( const Vecteur where );

  ///add a microtubule at clamp clampIndex radial to the sphere
  Microtub* addMicrotub( const int clampIndex );

  ///add a microtubule at clamp clampIndex with nuFocus located at place and direction dir
  Microtub* addMicrotub( const int clampIndex, const Vecteur& place, const Vecteur& dir );

  ///check if a tube belongs to a bundle: if yes, don't delete it (called by Microtub::step)
  bool      isDeletable( Microtub* mt ) const;

  ///remove a point from the nucleus and from the array pspts
  int       removePoint( const int ptIndex );

  //------------------- technical functions and mathematics ------------------
  
  ///orthogonalize the reference Vectors, done together with reshape
  void      orthogonalizeRef();        

  //------------------- methods for the projection ---------------------------
  
  ///allocate memory for the projection matrices
  void      allocateProjection();

  ///calculates the translational speed of points in Y (rigid body motion), for the forces given in X, scaled by sc
  void      setTransSpeedsFromForces(const real* X, real* Y, const real sc=1.0) const;
  
  ///calculates the rotational speed of points in Y (rigid body motion), for the forces given in X, scaled by sc
  void      setRotSpeedsFromForces(const real* X, real* Y, const real sc=1.0) const;
  
  ///calculates the speed of MTOCs on the sphere in Y, for the forces given in X, scaled by sc
  void      setMTOCSpeedsFromForces(const real* X, real* Y, const real sc=1.0) const;
  
  //------------------- methods for projection correction --------------------
  
  ///add the correction terms for the rigid body rotation, perform Y <- Y + P' * X;
  void      addRotProjectionDiff(const real* X, real* Y);
  
  ///add the correction terms MTOC motion, perform Y <- Y + P' * X;
  void      addMTOCProjectionDiff(const real* X, real* Y);
 

public:


  //--------------------------------------------------------------------------
  //=========================== public attributes ============================
  //--------------------------------------------------------------------------


  ///Space where nuclei are happy to live
  static Space*       space;
  
  ///number and coordinates of the reference points
  static const int    nbrefpts;
  #if( DIM > 1 )
  Vecteur   ref1;
  #endif
  #if( DIM > 2 )
  Vecteur   ref2;
  Vecteur   ref3;
  #endif
      

  //--------------------------------------------------------------------------
  //============================= public methods =============================
  //--------------------------------------------------------------------------

  //------------------- construction and destruction -------------------------

  ///default constructor
  Nucleus();

  ///constructor that uses parameters to create a nucleus
  Nucleus(int dummy);

  ///destructor
  virtual ~Nucleus();
    
  ///add a new point on the nucleus with a grafted attached to it
  void        addPointWithGrafted(Vecteur w, int ty);
  
  ///link a new grafted into the nucleus list of grafted
  void        linkGrafted(Grafted * gh) { nugrafteds.pushFront( gh ); }
  
  ///decide what to do, if a microtubule was deleted
  int         deregisterMicrotub(Microtub * mt);
  
  ///calculate the three mobilities: translational, rotational and point mobility
  void        setMobility();
  
  //------------------- querying the nucleus ---------------------------------
  
  ///position of center of gravity (returns the center of the sphere)
  Vecteur     getPosition()                     const { return whereP(0); }
  
  ///return the radius of the nucleus
  real        radius()                          const { return nuRadius; }

  ///return the number of grafted attached to this nucleus
  int         nbGrafteds()                      const { return nugrafteds.size(); }
  
  ///returns a pointer to the ii-th grafted
  Grafted*    getGrafted(const int ii)          const { return nugrafteds[ ii ]; }
  
  ///should not be used, since we have three different mobilities in the nucleus! This only returns psMobility.
  real        getMobility(const int d)          const
    { MSG.warning("Nucleus::getMobility","This function should not be used. There are three different mobilities in the nucleus!");
      return psMobility; }
  
  //------------------- space related functions ------------------------------
  
  ///returns true if the space is confined
  bool        isConfined()                      const { return space->isConfined(); }
  
  ///true if point pp is inside the space
//  bool        insidePoint(const int pp)         const;
  bool        insidePoint(const int pp)         const { return (pp > 0) || space->isInside(whereP(pp)); }
  
  ///project on the space's edge
//  void        projectPoint(int pp, real p[DIM]) const;
  void        projectPoint(int pp, real p[DIM]) const { space->project(whereP(pp), p ); }
  
  ///inside function for space: returns how many MT and Nucleus points are inside
  int         insidePosition( const Space * s ) const;
  
  ///modulo function for periodic space (including nucleus and attached microtubules)
  void        moduloPosition();
  
  ///modulo around the center of gravity 
  void        moduloPositionG();
  
  ///move the nucleus with its associated microtubs
  void        translatePosition( const Vecteur & T );
  
  ///set the space for the nucleus
  static void setSpace()
  {
    real nucbox[ParamSim::MAX];
    real nucboxinflate;
  
    //set defaults
    for( int ii = 0; ii < ParamSim::MAX; ii++ ) {
      nucbox[ii] = MP.boxsize[ii];
    }
    nucboxinflate = MP.boxinflate;
  
    switch( MP.boxshape ) {
      
      case SHAPE_OVAL:
      case SHAPE_BANANA:
        if ( MP.nuradius > MP.boxsize[1] ) {
          if ( MP.numax > 0 )
            MSG.error("Nucleus::setSpace","Box too small for the nucleus!");
        } else
          nucbox[1] = MP.boxsize[1] - MP.nuradius;
        nucbox[0] = MP.boxsize[0];
        break;
      
      case SHAPE_TEE:
        nucboxinflate -= MP.nuradius;
        break;
      
      default:
        break;
    }
    
    space=newSpace(MP.boxshape, nucbox, nucboxinflate);
  }
  
  ///a pointer to the space where the nucleus is living
  const Space* getSpace() const { return space; }
  
  //------------------- technical functions and mathematics ------------------

  ///we override the next() derived from Node to fix the type
  Nucleus*    next()                            const { return static_cast<Nucleus*>( son ); }

  ///we override the prev() derived from Node to fix the type
  Nucleus*    prev()                            const { return static_cast<Nucleus*>( dad ); }
  
  ///set two PointExact that clamp the Microtub at indx to a point on the nucleus, return false if clamp should not be done
  bool        setPointClamp( PointExact *, PointExact *, const int indx ) const;
  
  ///set four PointInterpolated that clamp the Microtub at indx to another microtub, return false if clamp should not be done
  bool        setMTClamp( PointInterpolated *, PointInterpolated *, PointInterpolated *, PointInterpolated *, const int indx ) const;
  
  ///this can't be used in nucleus since we have multiple mobilities, multiplication is done in set[]SpeedsFromForces
  void        multMobilityVect(real* V, real sc=1.0) const;
  
  ///this can't be used in nucleus since we have multiple mobilities, addBrownianMotion is used instead  
  real        addBrownianForces(real * force) const;
  
  ///after projection, add the Brownian Motion to argument (this is used instead of addBrownianForces)
  real        addBrownianMotion(real * rhs) const;
  
  ///precondition()=0 tells solve() not to use preconditionning on this block
  bool        precondition()                    const { return true; }

  ///overrides reshape in pointset to call one of the two reshape functions below
  void        reshape() { reshapeWithConservedCenterOfGravity(); }
  
  ///restore the reference shape and conserve the center of gravity
  void        reshapeWithConservedCenterOfGravity();
  
  ///restore the reference shape and conserve the center point, i.e. the position of the sphere
  void        reshapeWithConservedCenter();
  
  //------------------- methods for the projection ---------------------------
  
  ///prepare for constrained projection
  void        prepareProjection();
    
  ///this is not used in nucleus, modified setSpeedsFromForces instead
  void        setProjectedForces(const real * X, real * Y) const;
  
  ///calculates the speed of points in Y, for the forces given in X, scaled by sc
  void        setSpeedsFromForces(const real* X, real* Y, const real sc=1.0) const;
  
  //------------------- methods for projection correction --------------------
  
  ///set the projection correction term, P', from the given forces in X
  void        prepareProjectionDiff(const real* X);
  
  ///this can't be used in nucleus since we have multiple mobilities, setSpeedsWithCorrection is modified instead
  void        addProjectionDiff(const real* X, real* Y) const;

  ///calculates the speed of points in Y, for the forces given in X, adding the correction term P'  
  void        setSpeedsWithCorrection(const real* X, real* Z, real* Y, const real sc=1.0) const;
  
  ///test function for adding the projection correction - not used at the moment
  void        addSpeedDiffFromForces( const real* X, real* Y, real* FORCE, const real sc=1.0 ) const;
  
  //------------------- reading and writing the nucleus ----------------------
  
  ///write the nucleus with it's clamps to IO
  void        write();
  
  ///read the nucleus with it's clamps from IO
  void        read();
};


//=============================================================================
//HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
//=============================================================================


///Microtubules on the nucleus can be clamped to points on the sphere or to other Microtubules.
struct Nucleus::Clamp
{
  int clamp;   ///<index of the Nucleus-point in pspts where mt can be attached
  int mtIndex; ///<index of another mt in moMT where this mt can be attached
  int bundled; ///<true if this microtub belongs to a bundle

  ///creator calls clear(), just for beauty
  Clamp()
  {
    clear();
  }


  ///clear members
  void clear()
  {
    clamp   = -1;
    mtIndex = -1;
    bundled =  0;
  }
  

  ///set all members
  void set( int cl, int indx, int bundl )
  {
    clamp   = cl;
    mtIndex = indx;
    bundled = bundl;
  }


  ///set members clamp
  void setClamp( int cl )
  {
    clamp   = cl;
  }
  
  
  ///set member mt
  void setMT( int indx )
  {
    mtIndex = indx;
  }


  ///set member bundled
  void setBundle( int bundl )
  {
    bundled = bundl;
  }
};

#endif     // NUCLEUS_H
