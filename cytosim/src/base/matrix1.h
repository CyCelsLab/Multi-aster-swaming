//RCS: $Id: matrix1.h,v 2.1 2004/09/25 20:45:17 nedelec Exp $
//--------------------------------matrix1.h---------------------------------
#ifndef MATRIX1_H
#define MATRIX1_H

#define SZ 1
#define MATRIX Matrix11
#include "matrixbase.h"
#undef SZ
#undef MATRIX

//-------------------------------------------------------    
#include "vecteur1.h"

inline Vecteur1 operator * (const Matrix11 & a, const real v[1]) { 
  return Vecteur1( a[0] * v[0], 0, 0);
}

inline Vecteur1 operator * (const Matrix11 & a, const Vecteur1 & v) { 
  return Vecteur1( a[0] * v[0], 0, 0);
}

#endif //_MATRIX_1
