//RCS: $Id: solid.h,v 2.15 2005/04/25 08:27:00 nedelec Exp $
// ----------------------------------solid.h--------------------------------
// ------------------------UNDEFORMABLE SET OF POINTS-----------------------

#ifndef SOLID_H
#define SOLID_H

#include "array.h"
#include "node.h"
#include "space.h"
#include "simobject.h"
class Grafted;

///A Solid is an undeformable set of points build from PointSet
class Solid : public SimObject
{
 private:

  ///Type of the solid, at the moment this is unused
  int            soType;
  
  ///array to store the reference shape of the solid, as coordinates
  real     *     soShape;
  
  ///the number of points when fixShape() was last called, used for verifications.
  int            soShapeSize;
  
  ///this a temporary used for solids with 2 points only
  Vecteur        soAxis; 

  ///matrix containing the momentum of inertia
  real           soMom[ DIM*DIM ];
  
  ///second momentum of the reference shape
  real           soShapeMomentum;
    
  ///The list of attached motors (Grafted)
  Array<Grafted *> soGrafteds;

 public:

  ///Space where solids are happy to live
  static Space   *    space;

  ///Base constructor
  void        solidConstructorBase();
  
  ///allocate memory to hold 'size' points
  int         allocate(const int size);

  ///free all memory allocated by allocate()
  void        deallocate();
  
  ///default creator
  Solid()                                    { solidConstructorBase(); }
  
  ///creator according to parameter's specifications
  Solid(int);
  
  ///destructor
  virtual ~Solid();
  
  ///we override the next() derived from Node to fix the type
  Solid *         next()           const     { return static_cast<Solid*>( son ); }
  ///we override the prev() derived from Node to fix the type
  Solid *         prev()           const     { return static_cast<Solid*>( dad ); }
  
  ///sets the mobility
  void        setMobility();
  
  ///prepare for constrained projection
  void        prepareProjection();
  
  ///constrained projection 
  void        setProjectedForces(const real *, real *) const;

  ///set the Projection correction term, P', from the given forces
  void        prepareProjectionDiff(const real * forces);
  
  ///addProjectionDiff(real *X, real *Y) should perform Y <- Y + P' * X;
  void        addProjectionDiff(const real * X, real * Y) const;
  
  ///monte-carlo step
  void        step();
  
  ///
  void        diffuse();

  ///set the reference shape as a copy of the current one
  void        fixShape();
  
  ///rescale the current shape to the same size as the reference one
  void        rescale();
  
  ///restore the reference shape in the place and orientation of the current one
  void        reshape_really();
  
  ///reshape() alternatively calls rescale() or reshape_really()
  void        reshape();

  
  ///confinement of the space
  bool        isConfined()                const   { return space->isConfined(); }
  
  ///inside for a point in the space
  bool        insidePoint(const int pp)   const   { return space->isInside( whereP(pp) ); }

  ///project on the space's edge
  void        projectPoint(const int pp, real r[DIM] )
                                          const   { space->project( whereP(pp), r ); }
  
  ///a pointer to the space where the solid is living
  const Space* getSpace()                 const   { return space; }

  ///add a new point to the solid, and attach a grafted to it
  void        addPointWithGrafted(Vecteur w, int ty);
  
  ///link a new grafted into the solid list of grafted
  void        linkGrafted(Grafted * gh)       { soGrafteds.pushFront( gh ); }
  
  ///return the number of grafted attached to this solid
  int         nbGrafteds()              const { return soGrafteds.size(); }
  
  ///returns the ii-th grafted
  Grafted *   getGrafted(const int ii)  const { return soGrafteds[ ii ]; }

  ///write the solid to IO
  void        write();
  
  ///read the solid from IO
  void        read();
};

#endif
