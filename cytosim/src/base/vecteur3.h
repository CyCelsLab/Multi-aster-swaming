//RCS: $Id: vecteur3.h,v 2.5 2005/02/14 21:52:04 nedelec Exp $
//----------------------------------Vecteur3.h------------------------
#ifndef VECTEUR3_H
#define VECTEUR3_H

#define VECTEUR Vecteur3
#include "vecteurbase.h"
#undef VECTEUR


Vecteur3::Vecteur3(const real & xx, const real & yy, const real & zz) { 
  XX = xx; 
  YY = yy; 
  ZZ = zz; 
}

Vecteur3::Vecteur3(const real b[]) { 
  XX = b[0]; 
  YY = b[1]; 
  ZZ = b[2]; 
}

void Vecteur3::set(const real & xx, const real & yy, const real & zz) {
  XX = xx; 
  YY = yy; 
  ZZ = zz; 
}

void Vecteur3::operator =(const real b[]) {
  XX = b[0]; 
  YY = b[1];
  ZZ = b[2];   
}

void Vecteur3::set(const real b[]) {
  XX = b[0]; 
  YY = b[1];
  ZZ = b[2]; 
}


void Vecteur3::clear() {
  XX=0;
  YY=0;
  ZZ=0;
}


void Vecteur3::setRandom() {
  XX = RNG.sreal();
  YY = RNG.sreal();
  ZZ = RNG.sreal();
}

void Vecteur3::setRandom(const real & a) {
  XX = a*RNG.sreal();
  YY = a*RNG.sreal();
  ZZ = a*RNG.sreal();
}

void Vecteur3::addRandom(const real & a) {
  XX += a*RNG.sreal();
  YY += a*RNG.sreal();
  ZZ += a*RNG.sreal();
}

Vecteur3 Vecteur3::random() {
  return Vecteur3( RNG.sreal(), RNG.sreal(), RNG.sreal());
}

Vecteur3 Vecteur3::random(const real & a) {
  return Vecteur3( a*RNG.sreal(), a*RNG.sreal(), a*RNG.sreal());
}

real Vecteur3::normSquare()  const { 
  return     (XX*XX+YY*YY+ZZ*ZZ); 
}

real Vecteur3::norm()    const {
  return sqrt(XX*XX+YY*YY+ZZ*ZZ); 
}

void Vecteur3::normalize() {
  real s = 1.0 / norm(); 
  XX *= s; 
  YY *= s; 
  ZZ *= s; 
}

void Vecteur3::normalizeAndScale(const real & s) {
  real ss = s / norm(); 
  XX *= ss; 
  YY *= ss;
  ZZ *= ss; 
}

void Vecteur3::normalizeAndDivide(const real & s) {
  real ss = 1.0 / ( s * norm() ); 
  XX *= ss; 
  YY *= ss;
  ZZ *= ss; 
}

void Vecteur3::normalizeAndVmult(const real & s, const real a[]) {
  real ss = s / norm(); 
  XX *= ss * a[0]; 
  YY *= ss * a[1]; 
  ZZ *= ss * a[2]; 
}

Vecteur3 Vecteur3::normalized() const {
  real s=norm();
  return Vecteur3(XX/s, YY/s, ZZ/s);
}

void Vecteur3::vmult(const real b[]) {
  XX *= b[0];
  YY *= b[1];
  ZZ *= b[2];
}

void Vecteur3::vdiv(const real b[]) { 
  XX /= b[0];
  YY /= b[1]; 
  ZZ /= b[2]; 
}

void Vecteur3::oppose() { 
  XX=-XX;
  YY=-YY;
  ZZ=-ZZ;
}
  

//------------------------------------------------------------------

Vecteur3 operator +(const Vecteur3 & a, const Vecteur3 & b) { 
  return Vecteur3(a.XX+b.XX, a.YY+b.YY, a.ZZ+b.ZZ);
}

Vecteur3 operator -(const Vecteur3 & a, const Vecteur3 & b) { 
  return Vecteur3(a.XX-b.XX,a.YY-b.YY,a.ZZ-b.ZZ);
}

Vecteur3 operator -(const Vecteur3 & b) {
  return Vecteur3(-b.XX,-b.YY,-b.ZZ);
}


Vecteur3 operator ^(const Vecteur3 & a, const Vecteur3 & b) { 
  return Vecteur3(a.YY * b.ZZ - a.ZZ * b.YY,
                  a.ZZ * b.XX - a.XX * b.ZZ,
                  a.XX * b.YY - a.YY * b.XX);
}

Vecteur3 operator ^(const Vecteur3 & a, const real * b) { 
  return Vecteur3(a.YY * b[2] - a.ZZ * b[1],
                  a.ZZ * b[0] - a.XX * b[2],
                  a.XX * b[1] - a.YY * b[0]);
}

real operator *(const Vecteur3 & a, const Vecteur3 & b) {  
  return (a.XX*b.XX+ a.YY*b.YY+ a.ZZ*b.ZZ);
}
 
Vecteur3 operator *(const Vecteur3 & a, const real & s) { 
  return Vecteur3(s*a.XX, s*a.YY, s*a.ZZ); 
}

Vecteur3 operator *(const real & s, const Vecteur3 & a) { 
  return Vecteur3(s*a.XX, s*a.YY, s*a.ZZ);
}

Vecteur3 operator /(const Vecteur3 & a, const real & s) { 
  return Vecteur3(a.XX/s, a.YY/s, a.ZZ/s);
}

Vecteur3 interpolate(const Vecteur3 & a, const real & x, const Vecteur3 & b) {
  return Vecteur3(a.XX+x*b.XX, a.YY+x*b.YY, a.ZZ+x*b.ZZ); 
}
  
real distanceSquare(const Vecteur3 & a, const Vecteur3 & b) {
  return (a.XX-b.XX)*(a.XX-b.XX) + (a.YY-b.YY)*(a.YY-b.YY) + (a.ZZ-b.ZZ)*(a.ZZ-b.ZZ);
}

real distance(const Vecteur3 & a, const Vecteur3 & b) {
  return sqrt((a.XX-b.XX)*(a.XX-b.XX) + (a.YY-b.YY)*(a.YY-b.YY) + (a.ZZ-b.ZZ)*(a.ZZ-b.ZZ));
}


void Vecteur3::operator +=(const Vecteur3 & b) {
  XX+=b.XX; 
  YY+=b.YY;
  ZZ+=b.ZZ;
}

void Vecteur3::operator +=(const real * b) {
  XX+=*b; 
  YY+=*(b+1);
  ZZ+=*(b+2);
}

void Vecteur3::operator -=(const Vecteur3 & b) {
  XX-=b.XX; 
  YY-=b.YY;
  ZZ-=b.ZZ;
}

void Vecteur3::operator *=(const real & b) {
  XX*=b;
  YY*=b;
  ZZ*=b;
}

void Vecteur3::operator /=(const real & b) {
  XX/=b;
  YY/=b;
  ZZ/=b;
}


int operator ==(const Vecteur3 & a, const Vecteur3 & b) {
  return ((a.XX==b.XX) && (a.YY==b.YY) && (a.ZZ==b.ZZ));
}

int operator !=(const Vecteur3 & a, const Vecteur3 & b) {
  return ((a.XX!=b.XX) || (a.YY!=b.YY) || (a.ZZ!=b.ZZ));
}

void Vecteur3::print( FILE * out) {
  fprintf( out, "( %7.3f %7.3f %7.3f )", XX, YY, ZZ); 
}

void Vecteur3::println( FILE * out) {
  fprintf( out, "( %7.3f %7.3f %7.3f )\n", XX, YY, ZZ); 
}

#endif // VECTEUR3_H
