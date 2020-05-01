//RCS: $Id: scalar_field.h,v 2.3 2005/05/04 15:34:34 clausen Exp $
// -------------------------------lattice.h--------------------------------

#ifndef SCALAR_FIELD_H
#define SCALAR_FIELD_H

#include "vecteur.h"
#include "assert_macro.h"
#include "map.h"


class ScalarField
{

 private:
  int sizeI;//x-direction size of scalarField
  int sizeJ;//y-direction size of scalarField
  int sizeK;//z-direction size of scalarField
//the i, j, k address
  int cellI; //location scalar field value to be accessed
  int cellJ; //location of scalar field
  int cellK; //location of scalar field
//CHAITANYA: missing- the appropriate x, y, z address

  real *field;


 protected:

 public:




   void setFieldValues(int initCode);

   //Alternatively set field around some "SOURCE" object/objects
   void setFieldValues(int initCode, Vecteur sourcePos );


  //--------------------------------------------------------------

  ///Constructor
  ScalarField();

  ///Destructor
  virtual ~ScalarField()                 { ; }

  real value(Vecteur v);

  void readFromFile( char *fileName );

  void writeToFile( char *fileName );

  void test(void) { printf("Scalarfield Test header worked\n"); }

  void stepField();

};

#endif
