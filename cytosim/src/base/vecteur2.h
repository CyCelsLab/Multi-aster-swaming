//RCS: $Id: vecteur2.h,v 2.4 2005/02/14 21:51:42 nedelec Exp $
//----------------------------------Vecteur2.h------------------------
#ifndef VECTEUR2_H
#define VECTEUR2_H

#define CROSS_PRODUCT_IS_REAL
#define VECTEUR Vecteur2
#include "vecteurbase.h"
#undef VECTEUR
#undef CROSS_PRODUCT_IS_REAL

Vecteur2::Vecteur2(const real & xx, const real & yy, const real &) { 
  XX = xx; 
  YY = yy; 
}

Vecteur2::Vecteur2(const real b[]) { 
  XX = b[0]; 
  YY = b[1]; 
}

void Vecteur2::set(const real & xx, const real & yy, const real &) {
  XX = xx; 
  YY = yy; 
}

void Vecteur2::operator =(const real b[]) {
  XX = b[0]; 
  YY = b[1]; 
}

void Vecteur2::set(const real b[]) {
  XX = b[0]; 
  YY = b[1];
}


void Vecteur2::clear() {
  XX=0;
  YY=0;
}

void Vecteur2::setRandom() {
  XX = RNG.sreal();
  YY = RNG.sreal();
}

void Vecteur2::setRandom(const real & a) {
  XX = a*RNG.sreal();
  YY = a*RNG.sreal();
}

void Vecteur2::addRandom(const real & a) {
  XX += a*RNG.sreal();
  YY += a*RNG.sreal();
}

Vecteur2 Vecteur2::random() {
  return Vecteur2( RNG.sreal(), RNG.sreal(), 0);
}

Vecteur2 Vecteur2::random(const real & a) {
  return Vecteur2( a*RNG.sreal(), a*RNG.sreal(), 0);
}

real Vecteur2::normSquare()  const { 
  return     (XX*XX+YY*YY); 
}

real Vecteur2::norm()    const {
  return sqrt(XX*XX+YY*YY); 
}

void Vecteur2::normalize() {
  real s = 1.0 / norm(); 
  XX *= s; 
  YY *= s; 
}

void Vecteur2::normalizeAndScale(const real & s) {
  real ss = s / norm(); 
  XX *= ss; 
  YY *= ss;
}

void Vecteur2::normalizeAndDivide(const real & s) {
  real ss = 1.0 / ( s * norm()); 
  XX *= ss; 
  YY *= ss;
}

void Vecteur2::normalizeAndVmult(const real & s, const real a[]) {
  real ss = s / norm(); 
  XX *= ss * a[0]; 
  YY *= ss * a[1]; 
}

Vecteur2 Vecteur2::normalized() const {
  real s=norm();
  return Vecteur2(XX/s, YY/s, 0);
}

void Vecteur2::vmult(const real b[]) {
  XX *= b[0];
  YY *= b[1];
}

void Vecteur2::vdiv(const real b[]) { 
  XX /= b[0];
  YY /= b[1]; 
}

void Vecteur2::oppose() { 
  XX=-XX;
  YY=-YY;
}
  

//------------------------------------------------------------------

Vecteur2 operator +(const Vecteur2 & a, const Vecteur2 & b) { 
  return Vecteur2(a.XX+b.XX, a.YY+b.YY, 0);
}

Vecteur2 operator -(const Vecteur2 & a, const Vecteur2 & b) { 
  return Vecteur2(a.XX-b.XX, a.YY-b.YY, 0);
}

Vecteur2 operator -(const Vecteur2 & b) {
  return Vecteur2(-b.XX, -b.YY, 0);
}

real operator ^(const Vecteur2 & a, const Vecteur2 & b) { 
  return a.XX * b.YY - a.YY * b.XX;
}

Vecteur2 operator ^(const real & a, const Vecteur2 & b) {
  return Vecteur2( -a*b.YY, a*b.XX, 0); 
}

Vecteur2 operator ^(const Vecteur2 & a, const real & b) {
  return Vecteur2( a.YY*b, -a.XX*b, 0); 
}


real operator *(const Vecteur2 & a, const Vecteur2 & b) {  
  return (a.XX*b.XX + a.YY*b.YY);
}
 
Vecteur2 operator *(const Vecteur2 & a, const real & s) { 
  return Vecteur2(s*a.XX, s*a.YY, 0); 
}

Vecteur2 operator *(const real & s, const Vecteur2 & a) { 
  return Vecteur2(s*a.XX, s*a.YY, 0);
}

Vecteur2 operator /(const Vecteur2 & a, const real & s) { 
  return Vecteur2(a.XX/s, a.YY/s, 0);
}


Vecteur2 interpolate(const Vecteur2 & a, const real & x, const Vecteur2 & b) {
  return Vecteur2(a.XX+x*b.XX, a.YY+x*b.YY, 0); 
}
  
real distanceSquare(const Vecteur2 & a, const Vecteur2 & b) {
  return (a.XX-b.XX)*(a.XX-b.XX) + (a.YY-b.YY)*(a.YY-b.YY);
}

real distance(const Vecteur2 & a, const Vecteur2 & b) {
  return sqrt( (a.XX-b.XX)*(a.XX-b.XX) + (a.YY-b.YY)*(a.YY-b.YY));
}


void Vecteur2::operator +=(const Vecteur2 & b) {
  XX+=b.XX;
  YY+=b.YY;
}

void Vecteur2::operator +=(const real * b) {
  XX+=*b;
  YY+=*(b+1);
}

void Vecteur2::operator -=(const Vecteur2 & b) {
  XX-=b.XX; 
  YY-=b.YY;
}

void Vecteur2::operator *=(const real & b) {
  XX*=b;
  YY*=b;
}

void Vecteur2::operator /=(const real & b) {
  XX/=b;
  YY/=b;
}


int operator ==(const Vecteur2 & a, const Vecteur2 & b) {
  return ((a.XX==b.XX) && (a.YY==b.YY));
}

int operator !=(const Vecteur2 & a, const Vecteur2 & b) {
  return ((a.XX!=b.XX) || (a.YY!=b.YY) );
}


void Vecteur2::print( FILE * out) {
  fprintf( out, "( %7.3f %7.3f 0 )", XX, YY); 
}

void Vecteur2::println( FILE * out) {
  fprintf( out, "( %7.3f %7.3f 0 )\n", XX, YY); 
}

#endif // VECTEUR3_H
