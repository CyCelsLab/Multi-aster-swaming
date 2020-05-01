//RCS: $Id: matrix2.h,v 2.1 2004/09/25 20:35:41 nedelec Exp $
//----------------------------------Matrix2.h----------------------------
#ifndef MATRIX2_H
#define MATRIX2_H

#define SZ 2
#define MATRIX Matrix22
#include "matrixbase.h"
#undef SZ
#undef MATRIX

//-------------------------------------------------------    
#include "vecteur2.h"

inline Vecteur2 operator * (const Matrix22 & a, const real v[2]) { 
  return Vecteur2( a[0] * v[0] + a[2] * v[1],
                   a[1] * v[0] + a[3] * v[1], 0);
}

inline Vecteur2 operator * (const Matrix22 & a, const Vecteur2 & v) { 
  return Vecteur2( a[0] * v[0] + a[2] * v[1],
                   a[1] * v[0] + a[3] * v[1], 0);
}


#endif // _MATRIX_2
