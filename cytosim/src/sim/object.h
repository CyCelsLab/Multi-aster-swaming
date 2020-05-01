//RCS: $Id: object.h,v 1.6 2005/04/10 21:47:30 nedelec Exp $

#ifndef OBJECT_H
#define OBJECT_H


#include "node.h"
#include "vecteur.h"
#include "transformation.h"
#include "space.h"

///Common (and virtual) base class for objects in the simulation
class Object : public Node
{

public:

  ///Creator
  Object() {}

  ///destructor
  virtual ~Object() {}

  //---------------------------------------------------------
  //           Position-related functions
  //---------------------------------------------------------

  ///return the position in space of the object
  virtual Vecteur getPosition() const = 0;

  ///translate object's position by the given vecteur
  virtual void translatePosition(const Vecteur &) = 0;

  ///apply given transformation in space to the object
  virtual void transformPosition(const Transformation &) = 0;


  ///set the object position to the given vecteur
  virtual void setPosition(const Vecteur & w)
  {
    translatePosition( w - getPosition());
  }

  ///apply given transformation in space to the object
  virtual void transformPositionG(const Transformation & T)
  {
    Vecteur G = getPosition();
    translatePosition( -G );
    transformPosition(  T );
    translatePosition(  G );
  }

  ///number of points inside Space s
  virtual int insidePosition(const Space * s) const = 0;

  ///perform modulo for all points
  virtual void moduloPosition() = 0;

  ///initial configurations within given space
  static Vecteur    initPosition(const Space * space, const int mode);

  ///CHAITANYA: initial configurations within given space, position can be explicitly set
  static Vecteur    initPosition(const Space * space, const int mode, const Vecteur initPos );

  ///initial configurations within given space, type can be set as a function of position
  static Vecteur    initPosition(const Space * space, const int mode, int & type);

  ///pointer to the space where the object is living
  virtual const Space* getSpace() const { return 0; }

  ///true if the Space is confined
  virtual bool  isConfined() const { return false; }

  //---------------------------------------------------------
  //                  I/O functions
  //---------------------------------------------------------

  ///read object from the file
  virtual void read() = 0;

  ///write object to file
  virtual void write() = 0;

};


#endif

