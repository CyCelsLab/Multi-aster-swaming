//RCS: $Id: hand.h,v 2.14 2005/04/22 18:55:03 nedelec Exp $
//==============================hand.h=================================

///partial declaration of class Microtub
#ifndef HAND_H
#define HAND_H

class Microtub;
#include "sim_param.h"
#include "object.h"
#include "point_microtub.h"

///Variable of reasons of microtubules to detach.
enum  ReasonToUpdate { 
  ATTACH            = 0,
  DETACH_SPONTAN    = 2,      ///< spontaneous / pend or poff
  DETACH_MAXFORCE   = 3,      ///< maximum Force
  DETACH_INTRALINK  = 4,      ///< intra mt-bridg e, i.e. a link on adjacent rods
  DETACH_MTSHRINKS  = 5,      ///< induced by filament shortening
  DETACH_HAMAXEND   = 6,      ///< limitation induced by hamaxperend
  DETACH_LIMAXDUP   = 7,      ///< limitation induced by comaxperpair
  DETACH_DESTRUCT   = 9,      ///< destruct or
  DETACH_FILEREAD   = 100     ///< a new state has been read, updating the current one
};

///type of the object to which a Hand can belong
enum  HandBaseType {
  HA_BASE_NONE     = 0,       ///< Initial value, 
  HA_BASE_CX_SIDE1 = 1,       ///< the Hand is the first side of a Complex
  HA_BASE_CX_SIDE2 = 2,       ///< the Hand is the second side of a Complex
  HA_BASE_GRAFTED  = 3        ///< the Hand belongs to a Grafted
};

enum  HandAction {
  PLUSSTEP = 1,               ///< step forward
  MINUSSTEP = 2,              ///< step backward
  STAY = 3,                   ///< do nothing
  DETACH = 4                  ///< unbind from microtubule
};


//----------------------------------------------------------------------------
/**the Node is used here to link into the microtubule's attached hands lists
not into the objects list in sim, as usual.
*/

/** \TODO: we should maybe have derived class for each type of hand ?
*/

///Hand can bind Microtub and move on them: they represent molecular motors such as Kinesin
class Hand : public Node, public PointMicrotub
{
  friend class Grafted;
  friend class Complex;

 private:

    ///type of Hand
    int                haType; 
  
    ///used for reporting transitions
    int                haState;            
  
    ///complex; box; particle
    HandBaseType       haBaseType;  
	 
    ///the object to which I belong
    Object         *   haBase;             
    
    ///time of binding/unbinding
    Step               haStateUpdateTime;           
    
    ///the key can be set to prevent certain kinds of attachments
    BindingKey         haBindingKey;

    // time since last step was taken - used in stepUnderForce
    real               hadtLastStep;
    
    // time for one step - calculated from hand speed and mt dimer distance at binding
    real               haStepdt;

public:

    ///we override the next() derived from Node to fix the type
    Hand *         next()        { return static_cast<Hand*>( son ); }
    ///we override the prev() derived from Node to fix the type
    Hand *         prev()        { return static_cast<Hand*>( dad ); }
        
    ///Resets parameters of Hand
    void  handConstructorBase(const HandBaseType, Object *, const int );
    
    ///Constructor
    Hand( HandBaseType , Object *, int type = 0 );
    
    ///Constructor - without object or type
    Hand()                       { handConstructorBase(HA_BASE_NONE, 0, 0); }
    
    ///Destructor
    virtual ~Hand();
    
    ///Returns type of Hand
    int            getType()       const { return haType; }
  
    
    ///true if the binding of this hand to the given rod is not allowed
    bool           bindingIsForbidden( const PointMicrotub & site ) const;
    
    ///attach the hand at the position described in site
    bool           attach( PointMicrotub & site );
    
    ///attach the hand to a abscissa on a Microtub
    bool           attachTo( const Microtub *, const real abscis, const MTEnd from=MT_ORIGIN );
    
    ///detach the hand from its current Microtub
    void           detach( const ReasonToUpdate );
    
    
    ///Hands of type MAGIC_TYPE are special: they never detach, and bind immediately
    bool           hasStaticType() const { return ( haType == MAGIC_TYPE ); }
            
    ///returns 2 if at end of microubule, 1 if not at end but attached, and 0 if not attached. 
    int            getState()      const { return isAttached()?(getEnd()?2:1):0; }
    
    ///sets the binding/unbinding state of the hand
    void           setState();
    
    ///Changes the bound/unbound state of each hand
    void           updateState(const ReasonToUpdate);
    
    ///true if the state was updated at this time step
    bool           isStateNew() const;
    
    
    
    ///Calls PointMicrotub::moveBy and updates getMT() -> lattice
    void           moveByCheckLattice( real deltaAbs );

    ///location of self
    Vecteur        whereHand() const { return PointInterpolated::where(); }   
    
    //place where the other side is attached
    Vecteur        whereOtherSide() const;
    
    ///Shortest vector between the grafted object and itself (the 'stretch' of the bond) in the periodic space
    Vecteur        getStretch() const;
        
    ///Returns the other Hand if this is part of a complexes, zero otherwise
    const Hand *   otherHand() const;
    
    ///Returns the other Microtubule, when this is part of a complex, and other hand is attached
    const Microtub * otherMicrotub() const;
    
    ///Returns 1 if hand1 or hand2 are bound, or if grafted, 0 otherwise
    int            otherBound() const;
    
    ///the binding key, used to prevent illegal attachments
    BindingKey     bindingKey() const { return haBindingKey; }
    
    ///the Monte-Carlo step for a hand which is not under any force
    void           step();
    
    ///the Monte-Carlo step for a hand which is pulled by an force
    void           stepUnderForce( const Vecteur & );
    
    ///read from IO
    void           read();
    
    ///checkBoundaries, and write to IO
    void           write();
    
    ///checking the consistency of the data
    int            looksWrong() const;
    //Force dependent stepping behaviour
    enum HandAction getAction(Vecteur force);     			//plus ended motor
	enum HandAction getActionMinus(Vecteur force);			//minus ended motor

	real getActionMinusNK(Vecteur force);			        // Neha Khetan, //minus ended motor for Gliding Assay discrete model of dynein
};

#endif

