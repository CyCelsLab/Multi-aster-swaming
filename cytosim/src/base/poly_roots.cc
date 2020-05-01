//
// C++ Implementation: poly_roots
//
// Description: 
//
//
// Author: Dietrich Foethke <foethke@embl.de>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
// CVS: $Id: poly_roots.cc,v 1.4 2005/04/11 08:42:32 foethke Exp $


#include "poly_roots.h"
#include "smath.h"


#define SWAP(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)
#define REAL(z)  ((z).dat[0])
#define IMAG(z)  ((z).dat[1])


// the code for finding the real roots of a cubic function was copied and
// modified from the sources of the GNU Scientific Library (GSL)
// http://www.gnu.org/software/gsl/
int PolyRoots::solveCubic(double a, double b, double c, double d, double sol[3])
{
  if( a != 0 ) {
    b /= a;
    c /= a;
    d /= a;
  } else return -1;
  
  double q = (b * b - 3 * c);
  double r = (2 * b * b * b - 9 * b * c + 27 * d);

  double Q = q / 9;
  double R = r / 54;

  double Q3 = Q * Q * Q;
  double R2 = R * R;

  double CR2 = 729 * r * r;
  double CQ3 = 2916 * q * q * q;

  if (R == 0 && Q == 0)
    {
      sol[0] = - b / 3 ;
      sol[1] = - b / 3 ;
      sol[2] = - b / 3 ;
      return 3 ;
    }
  else if (CR2 == CQ3) 
    {
      /* this test is actually R2 == Q3, written in a form suitable
         for exact computation with integers */

      /* Due to finite precision some double roots may be missed, and
         considered to be a pair of complex roots z = x +/- epsilon i
         close to the real axis. */

      double sqrtQ = sqrt (Q);

      if (R > 0)
        {
          sol[0] = -2 * sqrtQ  - b / 3;
          sol[1] = sqrtQ - b / 3;
          sol[2] = sqrtQ - b / 3;
        }
      else
        {
          sol[0] = - sqrtQ  - b / 3;
          sol[1] = - sqrtQ - b / 3;
          sol[2] = 2 * sqrtQ - b / 3;
        }
      return 3 ;
    }
  else if (CR2 < CQ3) /* equivalent to R2 < Q3 */
    {
      double sqrtQ = sqrt (Q);
      double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
      double theta = acos (R / sqrtQ3);
      double norm = -2 * sqrtQ;
      sol[0] = norm * cos (theta / 3) - b / 3;
      sol[1] = norm * cos ((theta + 2.0 * PI) / 3) - b / 3;
      sol[2] = norm * cos ((theta - 2.0 * PI) / 3) - b / 3;
      
      /* Sort sol[0], sol[1], sol[2] into increasing order */

      if (sol[0] > sol[1])
        SWAP(sol[0], sol[1]) ;
      
      if (sol[1] > sol[2])
        {
          SWAP(sol[1], sol[2]) ;
          
          if (sol[0] > sol[1])
            SWAP(sol[0], sol[1]) ;
        }
      
      return 3;
    }
  else
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0/3.0);
      double B = Q / A ;
      sol[0] = A + B - b / 3;
      return 1;
    }
}


// the code for finding the real roots of a cubic function was copied and
// modified from the sources of the GNU Scientific Library (GSL)
// http://www.gnu.org/software/gsl/
int PolyRoots::zSolveCubic( double a, double b, double c, double d, complex sol[3])
{
  if( a != 0 ) {
    b /= a;
    c /= a;
    d /= a;
  } else return -1;

  double q = (b * b - 3 * c);
  double r = (2 * b * b * b - 9 * b * c + 27 * d);

  double Q = q / 9;
  double R = r / 54;

  double Q3 = Q * Q * Q;
  double R2 = R * R;

  double CR2 = 729 * r * r;
  double CQ3 = 2916 * q * q * q;

  if (R == 0 && Q == 0)
    {
      REAL (sol[0]) = -b / 3;
      IMAG (sol[0]) = 0;
      REAL (sol[1]) = -b / 3;
      IMAG (sol[1]) = 0;
      REAL (sol[2]) = -b / 3;
      IMAG (sol[2]) = 0;
      return 3;
    }
  else if (CR2 == CQ3) 
    {
      /* this test is actually R2 == Q3, written in a form suitable
         for exact computation with integers */

      /* Due to finite precision some double roots may be missed, and
         will be considered to be a pair of complex roots z = x +/-
         epsilon i close to the real axis. */

      double sqrtQ = sqrt (Q);

      if (R > 0)
        {
          REAL (sol[0]) = -2 * sqrtQ - b / 3;
          IMAG (sol[0]) = 0;
          REAL (sol[1]) = sqrtQ - b / 3;
          IMAG (sol[1]) = 0;
          REAL (sol[2]) = sqrtQ - b / 3;
          IMAG (sol[2]) = 0;
        }
      else
        {
          REAL (sol[0]) = -sqrtQ - b / 3;
          IMAG (sol[0]) = 0;
          REAL (sol[1]) = -sqrtQ - b / 3;
          IMAG (sol[1]) = 0;
          REAL (sol[2]) = 2 * sqrtQ - b / 3;
          IMAG (sol[2]) = 0;
        }
      return 3;
    }
  else if (CR2 < CQ3)  /* equivalent to R2 < Q3 */
    {
      double sqrtQ = sqrt (Q);
      double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
      double theta = acos (R / sqrtQ3);
      double norm = -2 * sqrtQ;
      double r0 = norm * cos (theta / 3) - b / 3;
      double r1 = norm * cos ((theta + 2.0 * M_PI) / 3) - b / 3;
      double r2 = norm * cos ((theta - 2.0 * M_PI) / 3) - b / 3;

      /* Sort r0, r1, r2 into increasing order */

      if (r0 > r1)
        SWAP (r0, r1);

      if (r1 > r2)
        {
          SWAP (r1, r2);

          if (r0 > r1)
            SWAP (r0, r1);
        }

      REAL (sol[0]) = r0;
      IMAG (sol[0]) = 0;

      REAL (sol[1]) = r1;
      IMAG (sol[1]) = 0;

      REAL (sol[2]) = r2;
      IMAG (sol[2]) = 0;

      return 3;
    }
  else
    {
      double sgnR = (R >= 0 ? 1 : -1);
      double A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), 1.0 / 3.0);
      double B = Q / A;

      if (A + B < 0)
        {
          REAL (sol[0]) = A + B - b / 3;
          IMAG (sol[0]) = 0;

          REAL (sol[1]) = -0.5 * (A + B) - b / 3;
          IMAG (sol[1]) = -(sqrt (3.0) / 2.0) * fabs(A - B);

          REAL (sol[2]) = -0.5 * (A + B) - b / 3;
          IMAG (sol[2]) = (sqrt (3.0) / 2.0) * fabs(A - B);
        }
      else
        {
          REAL (sol[0]) = -0.5 * (A + B) - b / 3;
          IMAG (sol[0]) = -(sqrt (3.0) / 2.0) * fabs(A - B);

          REAL (sol[1]) = -0.5 * (A + B) - b / 3;
          IMAG (sol[1]) = (sqrt (3.0) / 2.0) * fabs(A - B);

          REAL (sol[2]) = A + B - b / 3;
          IMAG (sol[2]) = 0;
        }

      return 3;
    }
}


int PolyRoots::solveQuartic(double a, double b, double c, double d, double e, double sol[4])
{
  complex zSol[4];
  int numReal = 0;
  int numImag;
    
  numImag = zSolveQuartic( a, b, c, d, e, zSol );
  
  if( numImag == 4 ) {
    // count the number of real solutions
    for( int ii=0; ii<4; ii++ ) {
      if( IMAG(zSol[ii]) == 0 ) {
        sol[numReal] = REAL(zSol[ii]);
       numReal++;
     }
   }

   return numReal;
 }
 
 return -1;
}


int PolyRoots::zSolveQuartic(double a, double b, double c, double d, double e, complex sol[4])
{  
  double  f, g, h, r, s;
  complex zSolCubic[3];
  int     numSol;
  int     numReal = 0;
  int     zSol1, zSol2;

  // divide all terms by a which yields an equation of the form X^4 + bX^3 + cX^2 + dX + e = 0
  if( a != 0 ) {
    b = b/a;
    c /= a;
    d /= a;
    e /= a;
  } else return -1;
  
  // define new coefficients that lead to the equation Y^4 + fY^2 + gY + h = 0,
  // where Y = X - b/4;
  f = c - (b*b*0.375);
  g = d + (b*b*b*0.125) - (b*c*0.5);
  h = e - (b*b*b*b*3./256.) + (b*b*c*0.0625) - (b*d*0.25);
//  printf("f: %f, g: %f, h: %f\n", f, g, h);
  
  // if h=0, x=0 is a trivial solultion of the quartic and we only need to solve
  // the cubic equation Y^3 + fY + g = 0
  //
  // still has to be written...
  
  // if g=0, the quartic is actually a quadratic in Z=Y^2:  Z^2 + fZ + h = 0
  if( g == 0 ) {
    // the solutions of the quadratic equation are: Z = +-sqrt(f^2/4 - h) -f/2
    s = (f*f/4. - h);
    
    if( s > 0 ) {
      // the two solutions of the quadratic are real
      REAL(zSolCubic[0]) =  sqrt(s) - f/2.;
      IMAG(zSolCubic[0]) =  0.;
      REAL(zSolCubic[1]) = -sqrt(s) - f/2.;
      IMAG(zSolCubic[1]) =  0.;
    } else if( s == 0 ) {
      // the two solutions of the quadratic are real and the same
      REAL(zSolCubic[0]) = -f/2.;
      IMAG(zSolCubic[0]) =  0.;
      REAL(zSolCubic[1]) = -f/2.;
      IMAG(zSolCubic[1]) =  0.;      
    } else {
      // one real (x=0), and two complex solutions
      REAL(zSolCubic[0]) = -f/2.;
      IMAG(zSolCubic[0]) =  sqrt(-s);
      REAL(zSolCubic[1]) = -f/2.;
      IMAG(zSolCubic[1]) = -sqrt(-s);
    }
    
/*    printf( "Z1: %e+%ei, Z2: %e+%ei\n", \
            zSolCubic[0].dat[0], zSolCubic[0].dat[1], \
            zSolCubic[1].dat[0], zSolCubic[1].dat[1] );*/
    
    // take the square root of the two solutions for Z to get four solutions for Y
    if( IMAG(zSolCubic[0]) != 0 ) {
      real phi, norm;
      // square root of complex numbers
      phi  = atan2( IMAG(zSolCubic[0]), REAL(zSolCubic[0]) );
      norm = sqrt(REAL(zSolCubic[0])*REAL(zSolCubic[0]) + IMAG(zSolCubic[0])*IMAG(zSolCubic[0]));
      REAL(sol[0]) =  sqrt(norm) * cos(phi/2.);
      IMAG(sol[0]) =  sqrt(norm) * sin(phi/2.);
    } else {
      if( REAL(zSolCubic[0]) > 0 ) {
        REAL(sol[0]) = sqrt(REAL(zSolCubic[0]));
        IMAG(sol[0]) = 0;
      } else if( REAL(zSolCubic[0]) == 0 ) {
        REAL(sol[0]) = 0;
        IMAG(sol[0]) = 0;
      } else {
        REAL(sol[0]) = 0;
        IMAG(sol[0]) = sqrt(-REAL(zSolCubic[0]));
      }
    }
    REAL(sol[1])  = -REAL(sol[0]) - b/4.;
    IMAG(sol[1])  = -IMAG(sol[0]);
    REAL(sol[0]) -=  b/4.;
    
    if( IMAG(zSolCubic[1]) != 0 ) {
      real phi, norm;
      // square root of complex numbers
      phi  = atan2( IMAG(zSolCubic[1]), REAL(zSolCubic[1]) );
      norm = sqrt(REAL(zSolCubic[1])*REAL(zSolCubic[1]) + IMAG(zSolCubic[1])*IMAG(zSolCubic[1]));
      REAL(sol[2]) =  sqrt(norm) * cos(phi/2.);
      IMAG(sol[2]) =  sqrt(norm) * sin(phi/2.);      
    } else {
      if( REAL(zSolCubic[1]) > 0 ) {
        REAL(sol[2]) = sqrt(REAL(zSolCubic[1]));
        IMAG(sol[2]) = 0;
      } else if( REAL(zSolCubic[1]) == 0 ) {
        REAL(sol[2]) = 0;
        IMAG(sol[2]) = 0;
      } else {
        REAL(sol[2]) = 0;
        IMAG(sol[2]) = sqrt(-REAL(zSolCubic[1]));
      }
    }
    REAL(sol[3])  = -REAL(sol[2]) - b/4;
    IMAG(sol[3])  = -IMAG(sol[2]);
    REAL(sol[2]) -= b/4;
  
    return 4;
  
  } else {
    // solve the associated cubic: X^3 + (f/2)X^2 + ((f^2 - 4*h)/16)X - g*g/64 = 0
    // all solutions are non zero, since g !=0
    
    numSol = zSolveCubic( 1., (f/2.), ((f*f-4.*h)/16.), (-g*g/64.), zSolCubic );

/*    printf( "zSolCubic1: %e+%ei, zSolCubic2: %e+%ei, zSolCubic3: %e+%ei\n", \
            zSolCubic[0].dat[0], zSolCubic[0].dat[1], \
            zSolCubic[1].dat[0], zSolCubic[1].dat[1], \
            zSolCubic[2].dat[0], zSolCubic[2].dat[1] );*/
  
    // we always find three solutions, but maybe only one of them is real
    for( int ii=0; ii<3; ii++ ) {
      if( IMAG(zSolCubic[ii]) == 0 ) numReal++;
    }
//    printf("numReal: %d\n", numReal);
    
    if( numReal == 1 ) {
      // one real, two complex, we use the two complex to go on
      double phi, norm;
    
      // the two complex solutions must be complex conjugates
      // therefore it's enough to find one of them
      if( IMAG(zSolCubic[0]) == 0 ) {
        zSol1 = 1;
        // there must be two complex solutions, so this must be one of them
        if( IMAG(zSolCubic[zSol1]) == 0 ) return -1;
      } else {
        zSol1 = 0;
      }
    
      // square root of complex numbers
      phi  = atan2( IMAG(zSolCubic[zSol1]), REAL(zSolCubic[zSol1]) );
      norm = sqrt(REAL(zSolCubic[zSol1])*REAL(zSolCubic[zSol1]) + IMAG(zSolCubic[zSol1])*IMAG(zSolCubic[zSol1]));
      REAL(zSolCubic[zSol1]) = sqrt(norm) * cos(phi/2.);
      IMAG(zSolCubic[zSol1]) = sqrt(norm) * sin(phi/2.);
       
      r = -g / (8.*norm);
      s =  b / 4.;
//      printf("r: %f, s: %f\n", r, s);
    
      REAL(sol[0]) =  2*REAL(zSolCubic[zSol1]) + r - s;
      REAL(sol[1]) =                               - r - s;
      REAL(sol[2]) =                               - r - s;
      REAL(sol[3]) = -2*REAL(zSolCubic[zSol1]) + r - s;
      IMAG(sol[0]) =  0;
      IMAG(sol[1]) = -2*IMAG(zSolCubic[zSol1]);
      IMAG(sol[2]) =  2*IMAG(zSolCubic[zSol1]);
      IMAG(sol[3]) =  0;
    
      return 4;
      
    } else {
      // all real and all non zero, since g != 0
      
      if( REAL(zSolCubic[2]) != 0 ) zSol1 = 2;
      else if ( REAL(zSolCubic[1]) != 0 ) zSol1 = 1;
      else {
        printf("Not enough nonzero solutions to solve quatric equation!\n");
        exit(0);
      }
      if( zSol1 == 2 ) {
        if( REAL(zSolCubic[1]) != 0 ) zSol2 = 1;
        else if( REAL(zSolCubic[0]) != 0 ) zSol2 = 0;
        else {
          printf("Not enough nonzero solutions to solve quatric equation!\n");
          exit(0);
        }
      } else {
        if( REAL(zSolCubic[0]) != 0 ) zSol2 = 0;
        else {
          printf("Not enough nonzero solutions to solve quatric equation!\n");
          exit(0);
        }
      }

      // calculate square roots
      if( REAL(zSolCubic[zSol1]) < 0 ) {
        IMAG(zSolCubic[zSol1]) = sqrt( -REAL(zSolCubic[zSol1]) );
        REAL(zSolCubic[zSol1]) = 0.;
      } else
        REAL(zSolCubic[zSol1]) = sqrt(  REAL(zSolCubic[zSol1]) );
      if( REAL(zSolCubic[zSol2]) < 0 ) {
        IMAG(zSolCubic[zSol2]) = sqrt( -REAL(zSolCubic[zSol2]) );
        REAL(zSolCubic[zSol2]) = 0.;
      } else
        REAL(zSolCubic[zSol2]) = sqrt(  REAL(zSolCubic[zSol2]) );
      
/*      printf( "p: %e+%ei, q: %e+%ei\n", \
            REAL(zSolCubic[zSol1]), IMAG(zSolCubic[zSol1]), \
            REAL(zSolCubic[zSol2]), IMAG(zSolCubic[zSol2]) );*/
      
      complex pq = zSolCubic[zSol1]*zSolCubic[zSol2];
      complex r  = ( -g/8.) / pq;
      s = b/4.;
//      printf("r: %f+%fi, s: %f\n", REAL(r), IMAG(r), s);
      
      REAL(sol[0]) =  REAL(zSolCubic[zSol1]) + REAL(zSolCubic[zSol2]) + REAL(r) - s;
      REAL(sol[1]) =  REAL(zSolCubic[zSol1]) - REAL(zSolCubic[zSol2]) - REAL(r) - s;
      REAL(sol[2]) = -REAL(zSolCubic[zSol1]) + REAL(zSolCubic[zSol2]) - REAL(r) - s;
      REAL(sol[3]) = -REAL(zSolCubic[zSol1]) - REAL(zSolCubic[zSol2]) + REAL(r) - s;
      IMAG(sol[0]) =  IMAG(zSolCubic[zSol1]) + IMAG(zSolCubic[zSol2]) + IMAG(r);
      IMAG(sol[1]) =  IMAG(zSolCubic[zSol1]) - IMAG(zSolCubic[zSol2]) - IMAG(r);
      IMAG(sol[2]) = -IMAG(zSolCubic[zSol1]) + IMAG(zSolCubic[zSol2]) - IMAG(r);
      IMAG(sol[3]) = -IMAG(zSolCubic[zSol1]) - IMAG(zSolCubic[zSol2]) + IMAG(r);
      
      return 4;
    }
  }
}
