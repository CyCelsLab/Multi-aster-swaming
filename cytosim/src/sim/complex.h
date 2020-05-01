//RCS: $Id: complex.h,v 2.16 2005/04/22 18:54:08 nedelec Exp $
//---------------------------------------------------------------------
//===============================complex.h=============================
//---------------------------------------------------------------------

#ifndef COMPLEX_H
#define COMPLEX_H

#include "object.h"
#include "hand.h"
#include "space.h"

///A Complex is a set of two Hands linked by an elastic element
class Complex : public Object
{
    friend class Hand;
    friend class ComplexList;
    
public:
        
    ///the space where the complexes live
    static Space   * space;
    
private:
        
    ///type of complex, if = -1, then the complex is static
    int              cxType;

    ///position of complex
    Vecteur          cxPos;
    
    ///position of complex in previous frame (used in play)
    Vecteur          cxPos_old;
    
    ///first Hand
    Hand             cxHand1;
    
    ///second Hand
    Hand             cxHand2;
    
    ///state of the complex: free / bound / bridge
    int              cxState;
    
    ///step at which last change of state occured
    Step             cxStateUpdateTime;

    ///contructorBase
    void             complexConstructorBase(const int type, const int typeha1=0, const int typeha2=0, const Vecteur w = VZERO);

  public:

    ///we override the next() derived from Node to fix the type
    Complex *         next()    const   { return static_cast<Complex*>(son); }
    ///we override the prev() derived from Node to fix the type
    Complex *         prev()    const   { return static_cast<Complex*>(dad); }
    
    ///default constructor
    Complex()                           { complexConstructorBase(0); }
    ///constructor following the parameters instructions in MP.
    Complex(int type);
    ///constructor for a Complex of type type at position w
    Complex(int type, Vecteur w);
    ///constructor for a Complex with specified hand-type and position
    Complex(int type, int typeha1, int typeha2, Vecteur w=VZERO);
    
    ///destructor
    virtual ~Complex();

    ///copy operator
    Complex &      operator =(const Complex &);
    ///swap the two hands of the complex, keeping them attached where they are
    void       swapHands();
    ///swap the binding sites of the two hands of the complex,
    void       swapHandTypes();
        
    ///return the type of the complex
    int          getType()      const { assert( cxType >= 0 ); return cxType; }
    ///the link stiffness used in sim_solve.cc
    real    getStiffness()      const  { if (cxType != MAGIC_TYPE) return 1; else return MP.askmratio; }
    ///true if the complex is not binding anything
    bool          isFree()      const { return ( cxHand1.isFree()  && cxHand2.isFree()  ); }
    ///true if the complex is bound to one and only one microtub
    bool      isAttached()      const { return ( cxHand1.isAttached() || cxHand2.isAttached() ); }
    ///true if the complex is bound to two microtub
    bool        isBridge()      const { return ( cxHand1.isAttached() && cxHand2.isAttached() ); }
    ///number of bound hands
    int     nbBoundHands()      const { return ( cxHand1.isAttached() + cxHand2.isAttached() ); }
    ///number of free hands
    int      nbFreeHands()      const { return ( cxHand1.isFree()  + cxHand2.isFree()  ); }
    ///compute the distance between the two hands: cxHand2.whereHand() - cxHand1.whereHand()
    Vecteur   getStretch()      const;
    ///cosinus of the angle between the two microtubules attached by the hands
    real cosAngle()             const { return cxHand1.dirMT() * cxHand2.dirMT() ; }

    ///get the state, a combination of the states of the hands
    int         getState()      const { return cxHand1.getState() + 3 * cxHand2.getState(); }
    ///update cxState and cxStateUpdateTime
    void       setState(); 
    ///update cxState and cxStateUpdateTime, with a reason for this update
    void     updateState( int why );
    ///true if the state whas just updated
    bool      isStateNew() const;
        
    ///returns the current position of the complex
    Vecteur calculatePosition();
    ///returns the vecteur where position is stored
    Vecteur getPosition()                             const { return cxPos; }
    ///set the current position to a given vecteur
    void    setPosition(const Vecteur & T)                  { cxPos = T; }
    ///translate object's position by the given vecteur
    void    translatePosition(const Vecteur & T)            { cxPos += T; }
    ///apply given transformation in space to the object
    void    transformPosition(const Transformation & T)     { T.rightMultOverride(cxPos); }
    ///true if the complex is inside the space
    int     insidePosition(const Space * s = space)   const { return s->isInside(cxPos.getXYZ()); }
    ///modulo the current position vecteur in the space
    void    moduloPosition()                                { Space::modulo( cxPos ); }
    ///project the current position vecteur on the Edge
    void    projectPosition(real p[DIM])                    { space->project( cxPos, p ); }
    ///copies the current position to cxPos_old
    void    setPositionOld()                                { cxPos_old = cxPos; }
    ///returns the position of the complex when setPositionOld was called
    const Vecteur& getPositionOld()                   const { return cxPos_old; }
    ///a pointer to the space where the complex is living
    const Space* getSpace()                           const { return space; }
    
    ///perform one step for a free complex: diffusion
    void   stepFree();
    ///perform one step for a bound complex
    void   stepBound();
    ///perform one step for a Bridge complex
    void   stepBridge();    
    
    ///try to attach the Hands for a free complex
    void   tryToAttachFree();
    ///try to attach the Hand for a bound Complex
    void   tryToAttachBound();
            
   
    //propagate some functions for the cxHand1:
    const Hand &       getHand1()      const { return cxHand1; }
    Vecteur              where1()      const { return cxHand1.where(); }
    bool            isAttached1()      const { return cxHand1.isAttached(); }
    const Microtub *     getMT1()      const { return cxHand1.getMT(); }
    MTEnd               getEnd1()      const { return cxHand1.getEnd(); }
    int                getType1()      const { return cxHand1.getType(); }
    int               getState1()      const { return cxHand1.getState(); }

    //propagate some functions for the cxHand2:
    const Hand &       getHand2()      const { return cxHand2; }
    Vecteur              where2()      const { return cxHand2.where(); }
    bool            isAttached2()      const { return cxHand2.isAttached(); }
    const Microtub *     getMT2()      const { return cxHand2.getMT(); }
    MTEnd               getEnd2()      const { return cxHand2.getEnd(); }
    int                getType2()      const { return cxHand2.getType(); }
    int               getState2()      const { return cxHand2.getState(); }

    ///attach Hand1 at the given position
    void    attachTo1(const Microtub * mt, const real ab, const MTEnd from)  { cxHand1.attachTo(mt, ab, from); }
    ///attach Hand2 at the given position
    void    attachTo2(const Microtub * mt, const real ab, const MTEnd from)  { cxHand2.attachTo(mt, ab, from); }
    ///attach Hand1 at the given PointMicrotub
    int     attach1(PointMicrotub & site)                  { return cxHand1.attach(site); }
    ///attach Hand2 at the given PointMicrotub
    int     attach2(PointMicrotub & site)                  { return cxHand2.attach(site); }

        
    ///write to IO
    void    write();
    ///read from IO
    void    read();
};


#endif

