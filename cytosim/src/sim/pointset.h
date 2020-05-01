//RCS: $Id: pointset.h,v 2.26 2005/04/22 12:13:17 nedelec Exp $
// -------------------------------pointset.h--------------------------------



#ifndef POINTSET_H
#define POINTSET_H

#include "object.h"
#include "space.h"
#include "transformation.h"
#include "assert_macro.h"

class Matrix;

///A set of points forming a geometrical object: Fiber, Solid, Nucleus...
class PointSet : public Object
{

  /// allocation of memory is done by chuncks, for efficiency
  static const int PSCHUNCK = 4;

 private:

  ///allocation size of array pspts[]
  int             psallocated;
  
  ///number of points in the set
  int             psmaxpts;
  
  ///Constructor base
  void            pointsetConstructorBase();

 protected:

  ///pspts[] is of size psallocated, and contain maxpts valid points of size DIM
  real         *  pspts;                 //usable range = [ 0 to DIM*maxpts [
  
  ///The center of gravity is calculed in calculatePosition() and stored here:
  Vecteur         pscenter;

  //--------------------------------------------------------------
public:

  ///Constructor
  PointSet();
  
  ///Copy constructor
  PointSet(const PointSet &);
  
  ///Copy assignement
  PointSet & operator =(const PointSet &);
  
  ///Destructor
  virtual ~PointSet()                          { deallocate(); }

  //--------------------------------------------------------------

  ///allocate() is called by addPoints for ex.,  returns size if memory was allocated
  virtual int   allocate(const int size);     
  
  ///deallocate() free all memory allocated by allocate()
  virtual void  deallocate();
  
  //--------------------------------------------------------------

  ///No more points
  void          clearPoints()                  { psmaxpts = 0; }  
  
  /// set all coordinates to zero (to get better debug/testing)
  void          resetPoints();          

  ///Adds a point, returning the point array index allocated
  int           addPoint(const Vecteur & w);

  ///sets the number of points in the list
  virtual void  setNbPoints(const int np)      { allocate(np); psmaxpts = np; }
  
  ///Position of an exact point
  Vecteur       getPoint(const int p)   const  { assert((p>=0) && (p<psmaxpts)); return Vecteur( pspts + DIM*p ); }

  ///Sets the position of the point p with vecteur w
  void          setPoint(const int p, const Vecteur & w);
  
  ///moves the point p by adding vecteur w
  void          movePoint(const int p, const Vecteur & w);
  
  ///add random noise unniformly to every coordinate (use for testing)
  void          addNoise(const real howmuch);

  //--------------------------------------------------------------

  /// number of points
  int           nbPoints()             const { return psmaxpts; }

  /// the index of the last point = nbPoints() - 1
  int           lastPoint()            const { return psmaxpts - 1; }
  
  /// the index of the last point = nbPoints() - 1
  int           nbSegments()           const { return psmaxpts - 1; }
  
  /// the index of the last segment = nbPoints - 2
  int           lastSegment()          const { return psmaxpts - 2; }

  ///Returns the modifiable address of the array of points
  real       *  pts()                        { return pspts; }

  ///Returns a modifiable address to the coordinate of one points
  real       &  pts(const int p)             { assert((p>=0) && (p<DIM*psmaxpts)); return pspts[p]; }

  ///Returns the address of the array of points
  real const *  getPts()               const { return pspts; }

  ///the address where point ii starts, as a constant pointer
  const real *  pointAddr(const int p) const { assert((p>=0) && (p<psmaxpts)); return pspts + DIM*p; }
  
  ///Position of an exact point
  Vecteur       whereP(const int p)    const { assert((p>=0) && (p<psmaxpts)); return Vecteur( pspts + DIM*p ); }
    
  ///Returns a modifiable coordinate from the array of points
  real          coord(const int p)     const { assert((p>=0) && (p<DIM*psmaxpts)); return pspts[p]; }
  
  //---------------------------------------------------------
  //           Position-related functions
  //---------------------------------------------------------
  
  ///Position of centre of gravity
  virtual Vecteur getPosition() const;
  
  ///calculate the center of gravity in member pscenter
  Vecteur       calculatePosition()          { pscenter = getPosition(); return pscenter; }
  
  ///translate object (moves all the points by w)
  virtual void  translatePosition(const Vecteur & T);
  
  ///rotate object
  void          transformPosition(const Transformation & T);
  
  ///space related function: 
  virtual int   insidePosition(const Space * s) const;
  
  ///this uses the first point, skipping center of gravity calculation
  virtual void  moduloPosition();
  
  ///call modulo around the center of gravity 
  virtual void  moduloPositionG();
  
  ///true if the corresponding point is in the Space
  virtual bool  insidePoint(const int point_index)               const = 0;
  
  ///project the corresponding point in the Space border 
  virtual void  projectPoint(const int point_index, real p[DIM]) const = 0;
  
  
  //--------------------------------------------------------------
  
  ///write to file
  void          write();
  
  ///reads file
  void          read();

};

#endif
