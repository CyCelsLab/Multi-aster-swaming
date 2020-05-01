//
// C++ Interface: poly_roots
//
// Description: 
//
//
// Author: Dietrich Foethke <foethke@embl.de>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
// CVS: $Id: poly_roots.h,v 1.3 2005/04/08 10:42:02 foethke Exp $


#ifndef POLY_ROOTS_H
#define POLY_ROOTS_H


#include "types.h"


/// PolyRoots contains functions for finding the roots of cubic and quartic equations
namespace PolyRoots {


  /// complex data type
  typedef struct {
    double dat[2];
  } complex;


  /// multiplication of two complex numbers
  inline complex operator*(complex a, complex b) {
    complex c;
    c.dat[0] = a.dat[0]*b.dat[0]-a.dat[1]*b.dat[1];
    c.dat[1] = a.dat[0]*b.dat[1]+b.dat[0]*a.dat[1];
    return c;
  }

  /// division of a real number by a complex
  inline complex operator/(double a, complex b) {
    complex c;
    double denom = b.dat[0]*b.dat[0]+b.dat[1]*b.dat[1];
    c.dat[0] =  a*b.dat[0]/denom;
    c.dat[1] = -a*b.dat[1]/denom;
    return c;
  }

  /// find the real roots of a 3rd degree polynomial (cubic) equation
  int solveCubic(double a, double b, double c, double d, double sol[3]);


  /// find the complex roots of a 3rd degree polynomial (cubic) equation
  int zSolveCubic(double a, double b, double c, double d, complex sol[3]);


  /// find the real roots of a 4th degree polynomial (quartic) equation
  int solveQuartic(double a, double b, double c, double d, double e, double sol[4]);


  /// find the complex roots of a 4th degree polynomial (quartic) equation
  int zSolveQuartic(double a, double b, double c, double d, double e, complex sol[4]);
  
}

#endif // POLY_ROOTS_H
