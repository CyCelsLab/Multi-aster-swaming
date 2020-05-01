//RCS: $Id: aster.h,v 2.13 2005/04/22 12:12:23 nedelec Exp $

#ifndef ASTER_H
#define ASTER_H

#include "microtub.h"
#include "microtub_organizer.h"
class PointExact;
class PointInterpolated;
class Solid;

//-------------------------------------------------------------------
//---------------------------- ASTER --------------------------------
//-------------------------------------------------------------------



///Aster uses a Solid and Clamp to join a bunch of Microtub in a radial configuration
class Aster : public Object, public MicrotubOrganizer
{
  friend class AsterList;
  
private:
  
  ///accessory class to store information useful to make the links Solid - Microtub
  class         Clamp;

  ///number of Microtub in the aster
  int           asMTmax;
  
  ///this store the information needed to make the links between the solid and the Microtub
  Array<Clamp>  asClamp; ///< the list of Clamps to link the Solid and Microtub
  
  ///which kind of Microtub end is in the center: MT_PLUS_END or MT_MINUS_END
  MTEnd         asFocus;

  ///as solid is used to mechanically hold the Microtub at the center
  Solid       * asSolid;
  
  ///radius of the solid sustainning the central links
  real          asRadius;
    
  //-------------------------------------------------------------------
  
  ///base constructor
  void       asterConstructorBase(const MTEnd focus=MT_MINUS_END);
  
  ///set the points on asSolid, DIM specific
  void       setSolidAndClamps(int nbmts);
  
  ///creates and set the attached microtubs
  void       setMicrotubs(int nbMTs);
  
  ///add some grafted to asSolid
  void       addGrafted();

  ///we immediately re-create a new Microtub where one has disapeared
  int        deregisterMicrotub(Microtub * mt);
  
public:
  
  ///default constructor
  Aster()                             { asterConstructorBase(); }
  
  ///constructor for nb_mts Microtubs, at position where, polarity focus
  Aster( int nb_mts, Vecteur where, MTEnd focus = MT_MINUS_END );
  
  ///constructor according to MP, for inverted or normal asters
   /// Aster(MTEnd focus );  /// ca
  
  Aster( MTEnd focus ) ;  
  
  ///destructor  
  virtual ~Aster();
  
  
  ///we override the next() derived from Node to fix the type
  Aster *         next()        const  { return static_cast<Aster*>( son ); }
  ///we override the prev() derived from Node to fix the type
  Aster *         prev()        const  { return static_cast<Aster*>( dad ); }
  
  
  ///the polarity of the aster: which end are in the center
  MTEnd      getFocus()         const  { return asFocus;  }
  
  ///return the scaffolding solid
  Solid    * getSolid()         const  { return asSolid; }
  
  
  
  ///return the center of gravity from all MT central ends
  Vecteur    getPosition()      const  { return MicrotubOrganizer::getPosition( asFocus ); }
  
	int getAsmtmax(){ return asMTmax;}
  
  real getAnisotropy(); // this function calculates the anisotropy sum(len1)-sum(len2)/sum(len1 +len2)
  real getAsymmetry(); // this function sum(len1)/sum(len2)
	
  ///move the aster with its associated solid and microtubs
  virtual void       translatePosition( const Vecteur & T );
  
  ///transform the aster with its associated solid and microtubs
  virtual void       transformPosition( const Transformation & T );
  
  ///inside function for space: returns how many MT and Solid points are inside
  virtual int        insidePosition( const Space * s = Microtub::space ) const;
  
  ///modulo function for periodic space
  void       moduloPosition();
  
  ///set the Points for the first clamp of Microtub at indx, return false if clamp should not be done
  bool       setClamp1(PointExact *, PointExact *, const int indx) const;
  
  ///set the Points for the second clamp of Microtub at indx, return false if clamp should not be done
  bool       setClamp2(PointExact *, PointInterpolated *, const int indx) const;

  ///perform on Monte-carlo step: this could be used to implement nucleation
  void       step() {}


  ///read from IO
  void       read();
  
  ///write to IO
  void       write();
  
  ///debug function
  int        looksWrong()  const;
};



//-------------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------



///Aster::Clamp defines an attachment between the Solid and a Microtub
class Aster::Clamp {
public:
  int  clamp1;  ///< index of the Solid-point corresponding to MT end
  int  clamp2;  ///< index of the Solid-point corresponding to MT intermediate point
  real clampA;  ///< abscissa of the MT intermediate point
  
  ///clear member values
  void clear()
  {
    clamp1 = 0;
    clamp2 = 0;
    clampA = 0;
  }
  
  ///creator calls clear(), just for beauty
  Clamp() 
  {
    clear();
  }
  
  ///set member values
  void set( int c1, int c2, real ca )
  {
    clamp1 = c1;
    clamp2 = c2;
    clampA = ca;
  }
  
};

#endif //_ASTER_H
