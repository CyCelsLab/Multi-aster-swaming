//RCS: $Id: grafted.h,v 2.15 2005/04/22 18:54:43 nedelec Exp $
//==============================grafted.h=========================
#ifndef GRAFTED_H
#define GRAFTED_H

#include "object.h"
#include "hand.h"
#include "space.h"
#include "point_exact.h"

class PointSet;
class Nucleus;
class Solid;

//TODO: make the basetype of Grafted more automatic:
// the basetypes of the grafted:
enum GraftedType { 
  GH_BASE    ='X',   ///< attached to a fix position
  GH_SOLID   ='S',   ///< attached to a Solid
  GH_NUCLEUS ='N',   ///< attached to a Nucleus
  GH_CORTEX  ='C',   ///< attached to the cell boundaries
  GH_DIFFUSE ='D'    ///< a grafted which can freely diffuse
};


///a Grafted is one Hand and a Vecteur describing the attachment position
class Grafted : public Object
{
  //we are not associal: we have some friends, but not too many
  friend class GraftedList;

 public:

  //the space in which Grafteds live happily
  static Space  * space;
    
 private:
  
  ///the Motor in the grafted
  Hand            ghHand;
    
  ///the type of the object on which the Grafted is attached
  GraftedType     ghBaseType;
  
  ///the position of attachment for GH_BASE and GH_CORTEX
  Vecteur         ghPos;
  
  ///the old position (only used in play)
  Vecteur         ghPos_old; 

  ///the object on which attached
  ///\todo: use a PointInterpolated for ghBase
  PointExact      ghBase;

  ///resets all member's variables
  void            graftedConstructorBase(const int type, const Vecteur w = VZERO, const GraftedType kind = GH_BASE );

public:

  Grafted(const int type, const Vecteur w, const GraftedType kind=GH_BASE);

  //this creates a Grafted attached to a point of the given Solid:
  Grafted(const int type, const Solid * so, const int point_index );
  
  //this creates a Grafted attached to a point of the given Nucleus:
  Grafted(const int type, const Nucleus * nu, const int point_index );

  //TODO: allow interpolated positions between points on solids

  ///default constructor
  Grafted()                                         { graftedConstructorBase(0); };

  ///constructor according to MP.
  Grafted(int type);
  
  ///destructor
  virtual ~Grafted();

  ///we override the next() derived from Node to fix the type
  Grafted *      next()                      const  { return static_cast<Grafted*>( son ); }
  ///we override the prev() derived from Node to fix the type
  Grafted *      prev()                      const  { return static_cast<Grafted*>( dad ); }
  
  ///copy constructor
  Grafted&   operator = (const Grafted &);

  
  ///return the position in space of the object
  Vecteur getPosition()                      const  { return ghPos; }
  ///set the object position to the given vecteur
  void setPosition(const Vecteur & w)               { ghPos = w; };
  ///translate object's position by the given vecteur
  void translatePosition(const Vecteur & T)         { ghPos += T; }
  ///apply given transformation in space to the object
  void transformPosition(const Transformation & T)  { T.rightMultOverride( ghPos ); }    
  ///test that grafted is inside the space
  int insidePosition(const Space * s)        const  { return s->isInside( ghPos.getXYZ() ); }
  /// modulo the position of the grafted
  void moduloPosition()                             { Space::modulo( ghPos ); }
  ///a pointer to the space where the grafted is living
  const Space* getSpace()                    const  { return space; }

  /// the position when setPositionOld() was called
  const Vecteur& getPositionOld()            const  { return ghPos_old;}  
  ///copy the current position in and 'old' register
  void           setPositionOld()                   { ghPos_old = ghPos; }
  
  
  ///true if the Grafted is in fact linking two PointSet
  bool           isLink()                    const  { return ghBase.getPS() != 0; }
  ///true if the Grafted simulates a cortical motor
  int            isCortical()                const  { return ghBaseType == GH_CORTEX; }
  ///true if the Grafted simulates a diffusible one-motor
  int            isDiffusing()               const  { return ghBaseType == GH_DIFFUSE; }
  
  ///returns the state of the Hand
  int            getState()                  const  { return ghHand.getState(); }
  ///update the State of the Hand
  void           setState()                         { ghHand.setState(); }
  ///update the State of the Hand, with a specified reason
  void           updateState( const int why );
  ///true if the state was updated at this time step
  bool           isStateNew() const                 { return ghHand.isStateNew(); }

  ///attach the Hand to the given site
  int     attach(PointMicrotub & site)              { return ghHand.attach(site); }
  ///detach the Hand for the specified reason
  void    detach(const ReasonToUpdate why)          { ghHand.detach(why); }
  ///attach the Hand to a give Microtub
  void    attachTo(const Microtub * mt, const real ab, const MTEnd from) { ghHand.attachTo(mt, ab, from); }
  ///move the Hand to a specified position
  void    moveTo(const real ab, const MTEnd from)   { ghHand.moveTo(ab, from); }
  
  
  ///position of the Hand
  Vecteur  whereHand()                       const  { return ghHand.whereHand(); }
  ///the position of the Graft
  Vecteur  whereGrafted();
  ///the distance (Hand - Graft)
  Vecteur  getStretch();

  ///return the base of the Grafted, ie. where it is grafted
  const PointExact &   getBase()             const  { return ghBase; }
  ///the Hand
  const Hand &   getHand()                   const  { return ghHand; }
  
  ///the type of the Hand
  int     getType()                          const  { return ghHand.getType(); }  
  ///true if Hand attached
  bool    isAttached()                       const  { return ghHand.isAttached(); }
  ///true if Hand free
  bool    isFree()                           const  { return ghHand.isFree(); }
  ///the haEnd of the Hand
  MTEnd   getEnd()                           const  { return ghHand.getEnd(); }
  ///the link stiffness used in sim_solve.cc
  real    getStiffness()                     const  { if (ghHand.hasStaticType()) return MP.askmratio; else return 1.0; }
  
  ///Monte-Carlo step for a free Grafted
  void    stepFree();
  ///Monte-Carlo step for a bound Grafted
  void    stepBound();
  
  ///read from IO
  void    read();
  ///write to IO
  void    write();
};


#endif
