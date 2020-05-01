//RCS: $Id: vecteur1.h,v 2.5 2005/02/14 21:51:14 nedelec Exp $
//----------------------------------Vecteur1.h------------------------
#ifndef VECTEUR1_H
#define VECTEUR1_H

#define CROSS_PRODUCT_IS_REAL
#define VECTEUR Vecteur1
#include "vecteurbase.h"
#undef VECTEUR
#undef CROSS_PRODUCT_IS_REAL

Vecteur1::Vecteur1(const real & xx, const real &, const real &) { 
  XX = xx; 
}

Vecteur1::Vecteur1(const real b[]) { 
  XX = b[0]; 
}

void Vecteur1::set(const real & xx, const real &, const real &) {
  XX = xx; 
}

void Vecteur1::operator =(const real b[]) {
  XX = b[0]; 
}

void Vecteur1::set(const real b[]) {
  XX = b[0]; 
}


void Vecteur1::clear() {
  XX=0.;
}

void Vecteur1::setRandom() {
  XX = RNG.sreal();
}

void Vecteur1::setRandom(const real & a) {
  XX = a*RNG.sreal();
}

void Vecteur1::addRandom(const real & a) {
  XX += a*RNG.sreal();
}

Vecteur1 Vecteur1::random() {
  return Vecteur1( RNG.sreal(), 0, 0);
}

Vecteur1 Vecteur1::random(const real & a) {
  return Vecteur1( a*RNG.sreal(), 0, 0);
}

real Vecteur1::normSquare()  const { 
  return     (XX*XX); 
}

real Vecteur1::norm()    const {
  return sqrt(XX*XX); 
}

void Vecteur1::normalize() {
  real s = norm(); 
  XX /= s; 
}

void Vecteur1::normalizeAndScale(const real & s) {
  real ss = s / norm(); 
  XX *= ss; 
}

void Vecteur1::normalizeAndDivide(const real & s) {
  real ss = s * norm(); 
  XX /= ss; 
}

void Vecteur1::normalizeAndVmult(const real & s, const real a[]) {
  real ss = s / norm(); 
  XX *= ss * a[0]; 
}

Vecteur1 Vecteur1::normalized() const {
  real s=norm();
  return Vecteur1(XX/s, 0, 0);
}

void Vecteur1::vmult(const real b[]) {
  XX *= b[0];
}

void Vecteur1::vdiv(const real b[]) { 
  XX /= b[0];
}

void Vecteur1::oppose() { 
  XX=-XX;
}
  

//------------------------------------------------------------------

Vecteur1 operator +(const Vecteur1 & a, const Vecteur1 & b) { 
  return Vecteur1(a.XX+b.XX, 0, 0);
}

Vecteur1 operator -(const Vecteur1 & a, const Vecteur1 & b) { 
  return Vecteur1(a.XX-b.XX, 0, 0);
}

Vecteur1 operator -(const Vecteur1 & b) {
  return Vecteur1(-b.XX, 0, 0);
}


real operator ^(const Vecteur1 & a, const Vecteur1 & b) { 
  return 0.;
}

Vecteur1 operator ^(const real & a, const Vecteur1 & b) {
  return Vecteur1( 0, 0, 0); 
}

Vecteur1 operator ^(const Vecteur1 & a, const real & b) {
  return Vecteur1( 0, 0, 0); 
}


real operator *(const Vecteur1 & a, const Vecteur1 & b) {  
  return (a.XX*b.XX);
}
 
Vecteur1 operator *(const Vecteur1 & a, const real & s) { 
  return Vecteur1(s*a.XX, 0, 0); 
}

Vecteur1 operator *(const real & s, const Vecteur1 & a) { 
  return Vecteur1(s*a.XX, 0, 0);
}

Vecteur1 operator /(const Vecteur1 & a, const real & s) { 
  return Vecteur1(a.XX/s, 0, 0);
}


Vecteur1 interpolate(const Vecteur1 & a, const real & x, const Vecteur1 & b) {
  return Vecteur1(a.XX+x*b.XX, 0, 0); 
}
  
real distanceSquare(const Vecteur1 & a, const Vecteur1 & b) {
  return (a.XX-b.XX)*(a.XX-b.XX);
}

real distance(const Vecteur1 & a, const Vecteur1 & b) {
  return fabs(a.XX-b.XX);
}


void Vecteur1::operator +=(const Vecteur1 & b) {
  XX+=b.XX; 
}

void Vecteur1::operator +=(const real * b) {
  XX+=*b; 
}

void Vecteur1::operator -=(const Vecteur1 & b) {
  XX-=b.XX; 
}

void Vecteur1::operator *=(const real & b) {
  XX*=b;
}

void Vecteur1::operator /=(const real & b) {
  XX/=b;
}


int operator ==(const Vecteur1 & a, const Vecteur1 & b) {
  return (a.XX==b.XX);
}

int operator !=(const Vecteur1 & a, const Vecteur1 & b) {
  return (a.XX!=b.XX);
}


void Vecteur1::print( FILE * out) {
  fprintf( out, "( %7.3f 0 0 )", XX); 
}

void Vecteur1::println( FILE * out) {
  fprintf( out, "( %7.3f 0 0 )\n", XX); 
}

#endif // VECTEUR1_H
