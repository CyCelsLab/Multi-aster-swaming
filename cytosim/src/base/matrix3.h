//RCS: $Id: matrix3.h,v 2.0 2004/08/16 16:36:14 nedelec Exp $
//----------------------------------matrix3.h----------------------------

#ifndef MATRIX3_H
#define MATRIX3_H

#define MATRIX Matrix33
#define SZ 3
#include "matrixbase.h"
#undef MATRIX
#undef SZ

//-------------------------------------------------------    
#include "vecteur3.h"

inline Vecteur3 operator * (const Matrix33 & a, const real v[3]) { 
  Vecteur3 r;
  r.XX= a[0]*v[0] + a[3]*v[1] + a[6]*v[2]; 
  r.YY= a[1]*v[0] + a[4]*v[1] + a[7]*v[2]; 
  r.ZZ= a[2]*v[0] + a[5]*v[1] + a[8]*v[2]; 
  return r;
}

inline Vecteur3 operator * (const Matrix33 & a, const Vecteur3 & v) { 
  Vecteur3 r;
  r.XX= a[0]*v[0] + a[3]*v[1] + a[6]*v[2]; 
  r.YY= a[1]*v[0] + a[4]*v[1] + a[7]*v[2]; 
  r.ZZ= a[2]*v[0] + a[5]*v[1] + a[8]*v[2]; 
  return r;
}

#endif  // MATRIX3_H
