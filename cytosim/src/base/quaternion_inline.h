//RCS: $Id: quaternion_inline.h,v 2.4 2005/01/14 17:43:56 nedelec Exp $
// Class Quaternion:
// Copyright F. Nedelec, EMBL, Oct 2002   Email: nedelec@embl.de
//

#ifndef QUATERNION_H
#define QUATERNION_H

#include <cstdio>
#include "types.h"
#include "smath.h"

/** Quaternions are extensions of complex number in dimension 4.
a quaternion is =  q[0] + i * q[1] + j * q[2] + k * q[3],
where q[?] are four real numbers. q[0] is the real part, and
the other are imaginary.
other implementations may have different conventions!

i, j and k are the three unit quaternions extending the real axis.
They multiply as: i*i = -1, i*j = k and permutations of these rules
ATTENTION: quaternions multiplication is NOT conmutative : j*i = -k

unit quaternions are mostly useful to represent rotations in space:
The group of rotation in 3D space has three degrees of freedom.
Symmetric matrices can represent them using 6 scalar numbers, while
unit quaternion represent them using only 4 scalar numbers.
besides the economy, using quaternion reduces numerical drifts due
to round-off errors.
It is also easier to normalise a quaternion than a matrix, 
to restore a perfect rotation, or to interpolate two rotations.

The class Quaternion below implements the standard mathematical
operations, plus conversions to and from 3x3 real matrices, 
and to 4x4 transformation matrices used in OpenGL.
*/

/// Quaternion Class
class Quaternion 
{
private:
    
    //quaternions have four dimensions:
    real q[4];
    
    //q[] represents the quaternion  =  q[0] + i * q[1] + j * q[2] + k * q[3]
    //the real part is q[0], check that before adding anything to this code,
    //as public implementations on the web have various conventions!
    
public:
    
    Quaternion() {}
    virtual ~Quaternion() {}
    
    
    Quaternion(const real a, const real b=0, const real c=0, const real d=0) {
        q[0] = a; 
        q[1] = b; 
        q[2] = c; 
        q[3] = d; 
    }

    void set(const real a, const real b = 0, const real c = 0, const real d = 0)
    {
        q[0] = a; 
        q[1] = b;
        q[2] = c;
        q[3] = d; 
    }
 
    void reset()
    {
	 set( 0, 0, 0, 0);
    }

    operator real * () 
    {
        return q; 
    }
    
    real  & operator [] ( const int n )    
    {
        return q[ n ]; 
    }
    
    real    operator [] ( const int n ) const
    { 
        return q[n]; 
    }
    
    
    real normSquare() const
    {
        return q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    }
    
    real norm() const
    {
        return sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
    }
    
    void normalizeit()
    {
        real s = norm();
        q[ 0 ] /= s;
        q[ 1 ] /= s;
        q[ 2 ] /= s;
        q[ 3 ] /= s;
    }
    
    Quaternion normalized()
    {
        real s = norm();
        return Quaternion( q[ 0 ]/s, q[ 1 ]/s, q[ 2 ]/s, q[ 3 ]/s );
    }
    
    void conjugateit()
    {
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
    }    
    
    Quaternion conjugated() const
    { 
        return Quaternion( q[0], -q[1], -q[2], -q[3] );
    }
    
    void inverseit()
    {
        real x = -normSquare();
        q[0] /= -x;
        q[1] /=  x;
        q[2] /=  x;
        q[3] /=  x;
    }
    
    Quaternion inverted() const
    { 
        real x = -normSquare();
        return Quaternion( -q[0]/x, q[1]/x, q[2]/x, q[3]/x );
    }
    
    void opposeit()
    {
        q[0] = -q[0];
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
    }
    
    Quaternion  operator - () const
    { 
        return Quaternion(-q[0], -q[1], -q[2], -q[3]); 
    }
    
    
    void squareit()
    {
        real a = q[0], b = q[1], c = q[2], d = q[3];
        q[0] = a*a - b*b - c*c - d*d;
        a += a;
        q[1] = a*b;
        q[2] = a*c;
        q[3] = a*d;
    }    
    
    Quaternion squared() const
    {
        return Quaternion(q[0]*q[0] - q[1]*q[1] - q[2]*q[2] - q[3]*q[3],
                          2*q[0]*q[1], 
                          2*q[0]*q[2], 
                          2*q[0]*q[3]);
    }
    
    Quaternion operator * ( const real f ) const 
    {
        return Quaternion(q[0]*f, q[1]*f, q[2]*f, q[3]*f); 
    }
    
    Quaternion operator / ( const real f ) const 
    { 
        return Quaternion(q[0]/f, q[1]/f, q[2]/f, q[3]/f);
    }
    
    void operator += ( const real f )
    { 
        q[0] += f; 
    }
    
    void operator -= ( const real f )
    { 
        q[0] -= f; 
    }
    
    void operator *= ( const real f )
    { 
        q[0] *= f; 
        q[1] *= f; 
        q[2] *= f;
        q[3] *= f; 
    }
    
    void operator /= ( const real f )
    { 
        q[0] /= f; 
        q[1] /= f; 
        q[2] /= f;
        q[3] /= f; 
    } 
    
    Quaternion operator + ( const Quaternion &a ) const
    { 
        return Quaternion(q[0]+a.q[0], q[1]+a.q[1], q[2]+a.q[2], q[3]+a.q[3]); 
    }
    
    Quaternion operator - ( const Quaternion &a ) const 
    {
        return Quaternion(q[0]-a.q[0], q[1]-a.q[1], q[2]-a.q[2], q[3]-a.q[3]);
    }
    
    
    void operator += ( const Quaternion & a ) 
    {
        q[0] += a.q[0];
        q[1] += a.q[1];
        q[2] += a.q[2];
        q[3] += a.q[3];
    }
    
    void operator -= ( const Quaternion & a ) 
    {
        q[0] -= a.q[0];
        q[1] -= a.q[1];
        q[2] -= a.q[2];
        q[3] -= a.q[3];
    }
    
    
    //right multiplication: this <- this * a
    void rightMult( const Quaternion & a )
    {
        real q0 = q[0],   q1 = q[1],   q2 = q[2],   q3 = q[3];
        
        q[0] = q0 * a.q[0] - q1 * a.q[1] - q2 * a.q[2] - q3 * a.q[3];
        q[1] = q0 * a.q[1] + q1 * a.q[0] + q2 * a.q[3] - q3 * a.q[2];
        q[2] = q0 * a.q[2] - q1 * a.q[3] + q2 * a.q[0] + q3 * a.q[1];
        q[3] = q0 * a.q[3] + q1 * a.q[2] - q2 * a.q[1] + q3 * a.q[0];
    }
    
    //left multiplications: this <- a * this;
    void leftMult( const Quaternion & a )
    {
        real q0 = q[0],   q1 = q[1],   q2 = q[2],   q3 = q[3];
        
        q[0] = q0 * a.q[0] - q1 * a.q[1] - q2 * a.q[2] - q3 * a.q[3];
        q[1] = q0 * a.q[1] + q1 * a.q[0] - q2 * a.q[3] + q3 * a.q[2];
        q[2] = q0 * a.q[2] + q1 * a.q[3] + q2 * a.q[0] - q3 * a.q[1];
        q[3] = q0 * a.q[3] - q1 * a.q[2] + q2 * a.q[1] + q3 * a.q[0];
    }
    
    
    // so-called 'fast' implementation of the normal multiplication:
    // a few multiplications supressed, but a lot of additions added
    // since additions and multiplications cost about the same,
    // the gain is not clear to me. To be tested!
    void rightMult_fast( const Quaternion & a )
     {
        real E = (q[3] + q[1]) * (a.q[1] + a.q[2]);
        real F = (q[3] - q[1]) * (a.q[1] - a.q[2]);
        real G = (q[0] + q[2]) * (a.q[0] - a.q[3]);
        real H = (q[0] - q[2]) * (a.q[0] + a.q[3]);
        real A = F - E;
        real B = F + E;
        real C = (q[0] + q[1]) * (a.q[0] + a.q[1]);
        real D = (q[0] - q[1]) * (a.q[2] + a.q[3]);
        E = (q[3] + q[2]) * (a.q[0] - a.q[1]);
        F = (q[3] - q[2]) * (a.q[2] - a.q[3]);
        q[0] = F + (A + G + H) * 0.5;
        q[1] = C + (A - G - H) * 0.5;
        q[2] = D + (B + G - H) * 0.5;
        q[3] = E + (B - G + H) * 0.5;
     }
     
    void leftMult_fast( const Quaternion & a )
    {
        real E = (a.q[3] + a.q[1])*(q[1] + q[2]);
        real F = (a.q[3] - a.q[1])*(q[1] - q[2]);
        real G = (a.q[0] + a.q[2])*(q[0] - q[3]);
        real H = (a.q[0] - a.q[2])*(q[0] + q[3]);
        real A = F - E;
        real B = F + E;
        real C = (a.q[0] + a.q[1])*(q[0] + q[1]);
        real D = (a.q[0] - a.q[1])*(q[2] + q[3]);
        E = (a.q[3] + a.q[2])*(q[0] - q[1]);
        F = (a.q[3] - a.q[2])*(q[2] - q[3]);
        q[0] = F + (A + G + H) * 0.5;
        q[1] = C + (A - G - H) * 0.5;
        q[2] = D + (B + G - H) * 0.5;
        q[3] = E + (B - G + H) * 0.5;
    }
    
    void operator *= ( const Quaternion & a )
    {
        rightMult(a);
    }
    
    void operator /= ( const Quaternion & a )
    {
        rightMult( a.inverted() );
    }
    
    Quaternion operator * ( const Quaternion & a ) const
    {
        Quaternion result = *this;
        result.rightMult( a );
        return result;
    }
    
    Quaternion operator / ( const Quaternion & a ) const
    {
        Quaternion result = *this;
        result.rightMult( a.inverted() );
        return result;
    }
    
    void setThisMatrix33( real m[3*3] ) const
    //generate the rotation matrix from the quaternion
    //assuming the quaternion is of norm = 1
    {
        real rx, ry, rz, xx, yy, yz, xy, xz, zz, x2, y2, z2; 
        
        x2 = q[1] + q[1];
        y2 = q[2] + q[2];
        z2 = q[3] + q[3];
        
        rx = q[0] * x2; ry = q[0] * y2; rz = q[0] * z2;
        xx = q[1] * x2; xy = q[1] * y2; xz = q[1] * z2;
        yy = q[2] * y2; yz = q[2] * z2; zz = q[3] * z2;
        
        m[0 + 3*0] = 1.0 - (yy + zz);
        m[1 + 3*0] = xy + rz;
        m[2 + 3*0] = xz - ry;
        
        m[0 + 3*1] = xy - rz;
        m[1 + 3*1] = 1.0 - (xx + zz);
        m[2 + 3*1] = yz + rx;
        
        m[0 + 3*2] = xz + ry;
        m[1 + 3*2] = yz - rx;
        m[2 + 3*2] = 1.0 - (xx + yy);
    }
    
    void setFromMatrix33( const real m[3*3] )
    {
        real  s;
        real trace = m[0 + 3*0] + m[1 + 3*1] + m[2 + 3*2];
        
        // check the diagonal
        if (trace > 0) {
            s = sqrt( trace + 1.0 );
            q[0] = s * 0.5;
            s = 0.5 / s;
            q[1] = (m[2 + 3*1] - m[1 + 3*2]) * s;
            q[2] = (m[0 + 3*2] - m[2 + 3*0]) * s;
            q[3] = (m[1 + 3*0] - m[0 + 3*1]) * s;
        } else {                
            // trace is negative
            // find biggest coefficient on diagonal:
            int i = 0;
            if (m[1 + 3*1] > m[0 + 3*0]) i = 1;
            if (m[2 + 3*2] > m[i + 3*i]) i = 2;
            
            s = sqrt( 1.0 + 2*m[i + 3*i] - trace );
            q[i+1] = s * 0.5;
            if (s != 0) s = 0.5 / s;
            int j = (i+1) % 3;
            int k = (j+1) % 3;
            q[j+1] = s * ( m[j + 3*i] + m[i + 3*j] );
            q[k+1] = s * ( m[i + 3*k] + m[k + 3*i] );
            q[0]   = s * ( m[k + 3*j] - m[j + 3*k] );
        }
    }
    
    void setThisMatrix16( double m[16], const real trans[3]=0 ) const
    //produces a OpenGL transformation matrix of translation trans
    //followed by the rotation represented by the quaternion *this
    //assuming the quaternion is of norm = 1,
    {
        real rx, ry, rz, xx, yy, yz, xy, xz, zz, x2, y2, z2; 
        
        x2 = q[1] + q[1];
        y2 = q[2] + q[2];
        z2 = q[3] + q[3];
        
        rx = q[0] * x2; ry = q[0] * y2; rz = q[0] * z2;
        xx = q[1] * x2; xy = q[1] * y2; xz = q[1] * z2;
        yy = q[2] * y2; yz = q[2] * z2; zz = q[3] * z2;
        
        m[0 + 4*0] = 1.0 - (yy + zz);
        m[1 + 4*0] = xy + rz;
        m[2 + 4*0] = xz - ry;
        
        m[0 + 4*1] = xy - rz;
        m[1 + 4*1] = 1.0 - (xx + zz);
        m[2 + 4*1] = yz + rx;
        
        m[0 + 4*2] = xz + ry;
        m[1 + 4*2] = yz - rx;
        m[2 + 4*2] = 1.0 - (xx + yy);
        
        if ( trans ) {
            m[0 + 4*3] = trans[0];
            m[1 + 4*3] = trans[1];
            m[2 + 4*3] = trans[2];
        } else {
            m[0 + 4*3] = 0;
            m[1 + 4*3] = 0;
            m[2 + 4*3] = 0;		
        }
        
        m[3 + 4*0] = 0;
        m[3 + 4*1] = 0;
        m[3 + 4*2] = 0;
        m[3 + 4*3] = 1.0;
    }
    
    
    void setFromPolar( const real r, 
                         const real phi, const real theta, const real psi )
    {
        q[0] = r * cos(phi);                   //r*cos(phi)
        q[1] = (q[3]=r*sin(phi)) * cos(theta); //r*sin(phi)*cos(theta)
        q[2] = (q[3]*=sin(theta)) * cos(psi);  //r*sin(phi)*sin(theta)*cos(psi)
        q[3] *= sin(psi);                      //r*sin(phi)*sin(theta)*sin(psi)
    }
    
    void setThisPolar( real * r, real * phi, real * theta, real * psi )
    {
        *r     = norm();
        *phi   = acos( q[0] / *r );
        *theta = acos( q[1] / (*r * sin(*phi)) );
        *psi   = atan2( q[3], q[2] );
    }
    
    static Quaternion newFromPolar( const real r, 
                                    const real phi, const real theta, const real psi )
    {
        Quaternion result;
        result.setFromPolar(r, phi, theta, psi);
        return result;
    }
    
    void setFromAxisAngle( const real v[3], const real angle )
    {
        //we normalize v for more security
        real n = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
        real d = angle * 0.5;
        real sd = sin( d ) / n;
        q[0] = cos( d );
        q[1] = v[0] * sd;
        q[2] = v[1] * sd;
        q[3] = v[2] * sd;
    }
    
    void setFromMainAxisAngle( const int axis, const real angle )
    //along the unit axis of dimension 'axis' ( 0:x, 1:y, 2:z )
    {
        real d = angle * 0.5;
        q[0] = cos( d );
        q[1] = 0;
        q[2] = 0;
        q[3] = 0;
        q[axis+1] = sin( d );
    }
    
    void setFromAxisAngle( const real v[3] )
    //the angle is given by the norm of the vector v
    //safely returns 1 if v = { 0, 0, 0 }
    {
        real n = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
        real d = n * 0.5;
        real sd = sin( d );
        if ( n ) sd /= n;
        q[0] = cos( d );
        q[1] = v[0] * sd;
        q[2] = v[1] * sd;
        q[3] = v[2] * sd;
    }
    
    void computeAxisAngle( real v[3], real * angle ) const 
    {
        real sd = sqrt( q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
        * angle = 2 * atan2( sd, q[0] );
        if ( sd != 0 ) sd = 1.0 / sd;
        v[0]   = q[1] * sd;
        v[1]   = q[2] * sd;
        v[2]   = q[3] * sd;
    }
    
    real getAngle() const 
    {
        real sd = sqrt( q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
        return  2 * atan2( sd, q[0] );
    }
    
    
    Quaternion slerp( const Quaternion &a, const Quaternion &b, const real t )
    //interpolation between two rotations 'a' and 'b'.
    {
        // Calculate the cosine of the angle between the two
        real scale0;
        real scale1;
        real co = a.q[0]*b.q[0] + a.q[1]*b.q[1] + a.q[2]*b.q[2] + a.q[3]*b.q[3];
        
        // If the angle is significant, use the spherical interpolation
        if( (1 - fabs(co)) > 1e-6 ){
            real temp  = acos( fabs(co) );
            real sinus = sin ( temp );
            scale0 = sin( ( - t) * temp) / sinus;
            scale1 = sin(     t  * temp) / sinus;
        } else {  //  Else use the cheaper linear interpolation
            scale0 = 1 - t;
            scale1 = t;
        }
        
        if( co < 0 ) scale1 = -scale1;
        
        // Return the interpolated result
        return (a * scale0) + (b * scale1);
    }
    
    
    
    void print( FILE * out = stdout )
    {
        fprintf( out, "( %7.3f %7.3f %7.3f %7.3f )", q[0], q[1], q[2], q[3]);
    }  
    
    void println( FILE * out = stdout )
    {
        fprintf( out, "( %7.3f %7.3f %7.3f %7.3f )\n", q[0], q[1], q[2], q[3]);
    }
};

#endif
