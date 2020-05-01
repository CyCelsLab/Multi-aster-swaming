//RCS: $Id: vecteur.h,v 2.0 2004/08/16 16:39:25 nedelec Exp $
//===========================vecteur.h=============================
//Here we define the class Vecteur to be the one of Vecteur1,2 or 3
//with the dimension DIM. Vecteur1, 2 or 3 are still available


#ifndef VECTEUR_H
#define VECTEUR_H

#include "main.h"

#ifndef DIM
     #error DIM is not defined
#endif

#include "vecteur1.h"
#include "vecteur2.h"
#include "vecteur3.h"

#if (DIM==1)
   typedef Vecteur1 Vecteur;
   typedef real CrossProduct;
   #define MZERO  0
#endif


#if (DIM==2)
   typedef Vecteur2 Vecteur;
   typedef real CrossProduct;
   #define MZERO  0
#endif


#if (DIM==3)
   typedef Vecteur3 Vecteur;
   typedef Vecteur CrossProduct;
   #define MZERO  VZERO
#endif


/// constant vectors: zero, and the three axial unit vectors
const Vecteur VZERO(0,0,0);
const Vecteur VX(1,0,0);
const Vecteur VY(0,1,0);
const Vecteur VZ(0,0,1);


#endif // VECTEUR_H
