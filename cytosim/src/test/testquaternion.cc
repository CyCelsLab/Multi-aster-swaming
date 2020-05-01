//RCS: $Id: testquaternion.cc,v 2.0 2004/08/16 16:44:16 nedelec Exp $

#include <cstdio>
#include "quaternion.h"
#include "matrix3.h"

int main()
{
  Matrix33 m, n;
  Quaternion q, p;
 
  real angle = PI/6.0, b;
  real v[]={0,0,0};
  real w[]={0,0,0};
  real error, error_max;

  printf("------------------- rotations of PI/6 -----------------\n");

  for( int ii = 0; ii<3; ++ii )
    {
      v[0]=0; v[1]=0; v[2]=0;
      v[ii] = 1.;

      q.setFromAxisAngle( v, angle );
      q.computeAxisAngle( w, &b );

      q.print();    
      printf(" norm = %.2f\n", q.norm() );

      printf(" %.2f %.2f %.2f  angle %.2f", v[0],v[1],v[2], angle);
      printf(" ?=? %.2f %.2f %.2f  angle %.2f\n", w[0],w[1],w[2], b);

      q.setThisMatrix33( m );
      m.print();
      printf("rotation error = %f\n", m.maxDeviationFromRotation());

      n.setRotationAroundAxis( v, angle );
      n -= m;
      printf("difference with mat = %f\n", n.maxNorm());
    }

  printf("------------------- identity ---------------------------\n");

  m.setToOne();
  m.print();
  q.setFromMatrix33( m );
  q.println(stdout);

  printf("-------------- quat-quat multiplication ----------------\n");

  for( int ii = 0; ii < 4; ++ii )
    {
      for( int jj = 0; jj < 4; ++jj )
	{
	  p = Quaternion(0); 
	  p[ii] = 1;
	  q = Quaternion(0); 
	  q[jj] = 1;
	  printf(" p = "); p.print();
	  printf(" q = "); q.print();
	  printf(" p * q = "); (p*q).println();
	}
      printf("\n");
    }

  printf("------------- many conversion quat-mat-quat -------------\n");

  for( int ii = -1; ii < 2; ++ii )
  for( int jj = -1; jj < 2; ++jj )
  for( int kk = -1; kk < 2; ++kk )
    {
      error_max=0;
      v[0]=ii; v[1]=jj; v[2]=kk;
      for( angle = 0; angle < PI; angle += 1 )
        {
        p.setFromAxisAngle( v, angle );
        p.setThisMatrix33( m );
        q.setFromMatrix33( m );
        //p.print(); q.print();
        if ( q[0] * p[0] < 0 ) q = -q;
        q -= p;
        error = q.norm();
        if ( error > error_max ) error_max = error;
        }
    }
    printf("  max error = %f\n", error_max);
  
  printf("------------rotation mult. is not commutative-----------\n");
  
  for( int ii = 0; ii<3; ++ii )
  for( int jj = 0; jj<3; ++jj )
    {
      q.setFromMainAxisAngle( ii, angle );
      p.setFromMainAxisAngle( jj, angle );
      
      (q*p).print(stdout);
      (p*q).println(stdout);
    }

  printf("------------ rotation around main axis------------------\n");
  angle = PI/6;
  for( int ii = 0; ii<3; ++ii )
    {
      q.setFromMainAxisAngle( ii, angle );
      q.setThisMatrix33( n );
      n.print(stdout);
    }
      

  double openGLmatrix[16];
  q.setThisMatrix16(openGLmatrix);
  
  return 0;
}
