//RCS: $Id: fiber.h,v 1.4 2005/04/25 13:28:37 nedelec Exp $
//---------------------------- fiber.h ----------------------------------

#ifndef FIBER_H
#define FIBER_H

#include "array.h"
#include "vecteur.h"
#include "simobject.h"
#include "point_microtub.h"

///option to build the projection matrix explicitely:
//comparing the two methods is useful for debugging.
//( the matrix version is much SLOWER : do not enable )
//#define PROJECT_WITH_MATRIX

//--------------------- some incomplete declarations:
class Matrix;


/// Fiber is a variable length linear PointSet with rigidity. Base for Microtub
class Fiber : public SimObject
{
      
 protected:
  
  ///type of filament. currently unused
  int            mtType;
  
  ///actual section length
  real           mtcut;            
  
  ///ideal section length,(=MP.mtrodlength)
  real           mtcutbest;            

#ifdef CUT_WITH_CURVATURE
  ///number of time steps between each attempt to remove/add a point
  ///also number of states over which curvature information is averaged
  static const int MTCUT_PERIOD = 7;

  ///error due to the cutting at different steps
  real           mtcuterror;
  
  ///index into the mtcuterror[]
  int            mtcuterrorindx;
#endif
  
  ///abscisse of the minus-end
  real           mtabminus;            
	
  ///rigidity scaling factor
  real           mtrigidkm; 
  
  ///array of rods, used in Attachments algorithm
  Array<PointMicrotub> mtrods;
  
#ifdef PROJECT_WITH_MATRIX
  //--------------------- members used for projecting with a matrix:
  // these members are used in file microtub_projmat.cc

  ///allocation size for projection Matrix;
  int            allocatedM;     
  
  ///projection matrix
  real    *      mtP;  
  
  ///differential of projection matrix
  real    *      mtDiffP;     
  
  ///intermediate of calculus
  real    *      mtJJtiJ;     
  
#else
  //--------------------- members used without explicit matrix calculation:
  // these members are used in file microtub_proj.cc
 
  ///allocation size for proj. without matrix;
  int            allocatedP;      
  
  ///precalculated differential of points of size DIM*nbSegments
  real    *      mtdpts;          
  
  ///J*J', a nbSegments * nbSegments matrix. We store the diagonal and one off-diagonal
  real    *      mtJJt;           
  
  ///vector for the projection correction of size nbSegments
  real    *      mtJJtiJforce;
  
  ///a flag to say that mtJJtiJforce is entirely null
  bool           mtJJtiJforce_null;
  
#endif
  
  ///vector of lagrange multipliers associated with the constraints
  real    *      lagrangeMult;
  
  ///true if the values in lagrangeMult are accurate
  bool           lagrangeValid;
  
public:
  //--------------------- methods:

  ///resets the pointers to arrays, should be called only by the constructors
  void        fiberConstructorBase();
  
  ///resets the values of member variables
  void        reset();

  ///get the type
  int         getType()                                const  { return mtType; }
  
  ///set the type
  void        setType( const int t )                          { mtType = t; }
    
  ///sets the number of points in the rods
  void        setNbPoints(int);
  
  ///set microtubules as a line, MINUS_END at where, direction dir (will be normalized), keeping original length
  void        setStraight(const Vecteur where, const Vecteur dir);

  ///set microtubules, 'which' end (+, -, or middle) at where, keeping the original length
  void        setStraight(const Vecteur where, const Vecteur dir, const MTEnd which);

  ///set microtubules as a line, MINUS_END at where, direction dir, of length newLength
  void        setStraight(const Vecteur where, const Vecteur dir, const real newLength);
  
  ///set microtubules, 'which' end (+, -, or middle) at where, length newLength
  void        setStraight(const Vecteur where, const Vecteur dir, const real newLength, const MTEnd which);


  //---------------------Constructor/Destructor
  
  ///Constructor - base
  Fiber()    { fiberConstructorBase(); }
    
  ///Constructor - creates a new microtubule with a center of gravity at position w
  Fiber(Vecteur w);
  
  ///Constructor -creates a microtubule with given specifications (position w, ...)
  Fiber(Vecteur w, Vecteur dir, real len, MTEnd which = MT_MINUS_END);

  ///Destructor
  virtual ~Fiber();

  ///allocate function
  int  allocate(const int size);
  
  ///Returns a Microtubule site at position p. 
  PointMicrotub*    getRod(const int p); 

  //---------------------abscisse
  //linear point-set support the concept of abscisse:
  ///abscisse = the distance from the origin (or from an end) taken along the tube
  //---------------------
  
  ///the distance from the origin to the minus-end of the fiber
  real        abscissaMinusEnd() const { return mtabminus; }
  
  ///returns the distance of the specified end from the origin.
  real        abscissa(const MTEnd end)                                                const;
  
  ///returns the distance of the specified end from the specified origin/end
  real        abscissa(const MTEnd end, const MTEnd from)                              const;

  ///returns the distance from the origin of a point specified by rod
  real        abscissa(const real rod)                                                 const;
  
  ///returns the distance of the rod, from the specified origin/end
  real        abscissa(const real rod, const MTEnd from)                               const;

  ///set an interpolated point to represent a given end/origin
  void        setInterpolation(PointInterpolated *, const MTEnd ) const;
  
  ///set an interpolated point to a site given by a distance from a given origin, return if the point in at the end
  MTEnd       setInterpolation(PointInterpolated *, const real ab, const MTEnd from=MT_ORIGIN ) const;
  
  //---------------------Distance from a point to rod
  
  ///calculate the distance from a given point to a given rod on the microtubule
  /** The PointInterpolated given should have two values set: mPS and mPoint1
  to define on which MT segment the distance will be calculated
  the value si->mCoef will be set here as the projection of w on the segment */
  real        distanceToRod(PointInterpolated *, const Vecteur &)                      const;  
  
  ///distanceToRodEnd also considers the end points, not only projected points
  real        distanceToRodWithEnds(PointInterpolated *, const Vecteur &)              const; 

  ///distanceToMT considers all the rods in the Microtub
  real        distanceToMT(PointInterpolated *, const Vecteur &)                       const; 

  ///optimized version of distance
  real        distanceToRodFast(PointInterpolated *, const Vecteur &)                  const; 
    
  //---------------------  
  ///returns the position of the interpolated point on the microtubule
  Vecteur     where(const PointInterpolated *)       const;
  
  ///returns the position of one of the microtubule's end / origin
  Vecteur     whereEnd(const MTEnd which)            const;
  
  ///returns the position of a point specified by abscissa on the microtubule
  Vecteur     where(const real ab, const MTEnd from) const;
  
  ///returns the deterministic forces (without Brownian) on the end of the microtubule
  real        forceOnEnd(const MTEnd which)          const;
    
  //---------------------
  
  ///the current sectionning length: length of the segments
  real        rodLength()                const { return mtcut; };
  
  ///the total length of the fiber
  real        length()                   const { return nbSegments() * mtcut; }
  
  //---------------------
  ///difference of two consecutive points, p and p+1
  Vecteur     dpts(const int p)                        const;
  
  ///tangent vector to the fiber, ie. the normalized difference, at an exact point
  Vecteur     dirP(const int p)                        const;

  ///tangent vector to the fiber, at an Interpolated point
  Vecteur     dirMT(const PointInterpolated *)         const;
  
  ///tangent vector to the fiber, at an given abscissa
  Vecteur     dirMT(const real ab, const MTEnd from)   const;
  
  ///tangent vector to the fiber, at an end
  Vecteur     dirEnd(MTEnd which)                      const;
  
  //--------------------- Growing/Shrinking
  ///the maximum sinus of the consecutive rod angles: indicate the errors due to curvature
  real        minCosRodAngles();
  
  ///Cuts the microtubule to imposet nbrod+1 points
  void        recut(int nbrod);
  
  ///set the number of points to minimize abs( mtcut - mtcutbest )
  void        optimalCut();
  
  /// increase length of microtubule by 'lon', at the minus end
  void        grow1(const real lon);
  
  ///increase length of microtubule by 'lon', at the plus end
  void        grow2(const real lon);
  
  /// increase length of microtubule at specified 'which' end 
  void        growAtEnd(const MTEnd which, const real lon);

  //---------------------
  ///move points to satisfy length constraints while conserving center of gravity of tubule
  void        reshape();
  
  ///set the mobility
  void        setMobility();
  //---------------------  
  
  void        constructorProjection();
  void        allocateProjection();
  void        deallocateProjection();
  void        computeProjectionMat();
  void        prepareProjectionDiff(const real *);
  void        addProjectionDiff(const real *, real * ) const;
  void        prepareProjection();
  void        setProjectedForces(const real *, real *) const;
  
  ///function to access the Lagrange multipliers, which show the tension in the Tubes
  real        getLagrange(const int ii) const;

  ///function to reset the lagrangeMult multipliers (they are set in the solving step)
  void        resetLagrange();
  
  //---------------------
  
  void        addRigidity(const real *, real * ) const;
  void        setRigidityUp(Matrix &, const int offset ) const;
    
  ///precondition()=1 tells solve() to use preconditionning on this block:
  bool        precondition() const { return true; }
   
  ///basic consistency check (debug)
  int         looksWrong() const;        
};


#endif
