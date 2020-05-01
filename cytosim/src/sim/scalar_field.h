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

  int cellX; //location scalar field value to be accessed
  int cellY; //location of scalar field
  int cellZ; //location of scalar field

  real *field;

 protected:

 public:


  //--------------------------------------------------------------

  ///Constructor
  ScalarField();

  ///Destructor
  virtual ~ScalarField()                 { ; }

  int sizeX;//x-direction size of scalarField
  int sizeY;//y-direction size of scalarField
  int sizeZ;//z-direction size of scalarField

  void setFieldValues(int initCode = 1);//initialize the field
  //void setFmatrix(int initCode = 1); //initialize the field matrix

  real value(Vecteur v);

  void readFromFile(char *fileName);

  void writeToFile(char *fileName);

  void test(void) { printf("Scalarfield Test header worked\n"); }

  void stepField();

  real fmatrix[0][0][0];//field matrix


};

#endif
