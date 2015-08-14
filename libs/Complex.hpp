/* 
   Copyright (C) 1988 Free Software Foundation
   written by Doug Lea (dl@rocky.oswego.edu)

   This file is part of the GNU C++ Library.  This library is free
   software; you can redistribute it and/or modify it under the terms of
   the GNU Library General Public License as published by the Free
   Software Foundation; either version 2 of the License, or (at your
   option) any later version.	This library is distributed in the hope
   that it will be useful, but WITHOUT ANY WARRANTY; without even the
   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   PURPOSE.  See the GNU Library General Public License for more details.
   You should have received a copy of the GNU Library General Public
   License along with this library; if not, write to the Free Software
   Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
*/

#ifndef __Complex_h__
#define __Complex_h__ 1

#define __ATT_complex__

#include <iostream>
#include <cmath>

using std::istream;
using std::ostream;
using std::ws;

class Complex
{
#ifdef __ATT_complex__
 public:
#else
 protected:
#endif

  float_tt re;
  float_tt im;

 public:

  Complex() {}
  Complex(float_tt r, float_tt i=0) : re(r), im(i) {}
  Complex(const Complex& y) : re(y.re), im(y.im) {}
	
  ~Complex() {}

  float_tt real() const {return re;}
  float_tt imag() const {return im;}

  const Complex& operator = (const Complex& y);
	
  const Complex& operator += (const Complex& y);
  const Complex& operator += (float_tt y);
  const Complex& operator -= (const Complex& y);
  const Complex& operator -= (float_tt y);
  const Complex& operator *= (const Complex& y);
  const Complex& operator *= (float_tt y);
  const Complex& operator /= (const Complex& y); 
  const Complex& operator /= (float_tt y);
	
  void error(char* msg) const;
};

// inline members

inline const Complex& Complex::operator = (const Complex& y) 

{ 
  re = y.re; im = y.im; return *this; 
} 

inline const Complex& Complex::operator += (const Complex& y)
{ 
  re += y.re;  im += y.im; return *this; 
}

inline const Complex& Complex::operator += (float_tt y)
{ 
  re += y; return *this; 
}

inline const Complex& Complex::operator -= (const Complex& y)
{ 
  re -= y.re;  im -= y.im; return *this; 
}

inline const Complex& Complex::operator -= (float_tt y)
{ 
  re -= y; return *this; 
}

inline const Complex& Complex::operator *= (const Complex& y)
{  
  float_tt r = re * y.re - im * y.im;
  im = re * y.im + im * y.re; 
  re = r; 
  return *this; 
}

inline const Complex& Complex::operator *= (float_tt y)
{  
  re *= y; im *= y; return *this; 
}

inline const Complex& Complex::operator /= (const Complex& y)
{
  register float_tt t1,t2,t3;
  t2=1.0/(y.re*y.re+y.im*y.im);
  t1=t2*y.re; t2 *= y.im; t3=re;
  re *= t1; re += im*t2;
  im *= t1; im -= t3*t2;
  return *this;
}

inline const Complex& Complex::operator /= (float_tt y)
{
  re /= y;
  im /= y;
  return *this;
}

//	functions

inline int operator == (const Complex& x, const Complex& y)
{
  return x.re == y.re && x.im == y.im;
}

inline int operator == (const Complex& x, float_tt y)
{
  return x.im == 0.0 && x.re == y;
}

inline int operator != (const Complex& x, const Complex& y)
{
  return x.re != y.re || x.im != y.im;
}

inline int operator != (const Complex& x, float_tt y)
{
  return x.im != 0.0 || x.re != y;
}

inline Complex operator - (const Complex& x)
{
  return Complex(-x.re, -x.im);
}

inline Complex conj(const Complex& x)
{
  return Complex(x.re, -x.im);
}

inline Complex operator + (const Complex& x, const Complex& y)
{
  return Complex(x.re+y.re, x.im+y.im);
}

inline Complex operator + (const Complex& x, float_tt y)
{
  return Complex(x.re+y, x.im);
}

inline Complex operator + (float_tt x, const Complex& y)
{
  return Complex(x+y.re, y.im);
}

inline Complex operator - (const Complex& x, const Complex& y)
{
  return Complex(x.re-y.re, x.im-y.im);
}

inline Complex operator - (const Complex& x, float_tt y)
{
  return Complex(x.re-y, x.im);
}

inline Complex operator - (float_tt x, const Complex& y)
{
  return Complex(x-y.re, -y.im);
}

inline Complex operator * (const Complex& x, const Complex& y)
{
  return Complex(x.re*y.re-x.im*y.im, x.re*y.im+x.im*y.re);
}

inline Complex multconj(const Complex& x, const Complex& y)
{
  return Complex(x.re*y.re+x.im*y.im,x.im*y.re-x.re*y.im);
}

inline Complex operator * (const Complex& x, float_tt y)
{
  return Complex(x.re*y, x.im*y);
}

inline Complex operator * (float_tt x, const Complex& y)
{
  return Complex(x*y.re, x*y.im);
}

inline Complex operator / (const Complex& x, const Complex& y)
{
  register float_tt t1,t2;
  t2=1.0/(y.re*y.re+y.im*y.im);
  t1=t2*y.re; t2 *= y.im;
  return Complex(x.im*t2+x.re*t1, x.im*t1-x.re*t2);
}

inline Complex operator / (const Complex& x, float_tt y)
{
  return Complex(x.re/y,x.im/y);
}

inline Complex operator / (float_tt x, const Complex& y)
{
  register float_tt factor;
  factor=1.0/(y.re*y.re+y.im*y.im);
  return Complex(x*y.re*factor,-x*y.im*factor);
}

inline float_tt real(const Complex& x)
{
  return x.re;
}

inline float_tt imag(const Complex& x)
{
  return x.im;
}

inline float_tt abs2(const Complex& x)
{
  return x.re*x.re+x.im*x.im;
}

inline float_tt abs(const Complex& x)
{
  return sqrt(abs2(x));
}

inline float_tt arg(const Complex& x)
{
  return x.im != 0.0 ? atan2(x.im, x.re) : 0.0;
}

// Return the principal branch of the square root (non-negative float_tt part).
inline Complex sqrt(const Complex& x)
{
  float_tt mag=abs(x);
  if(mag == 0.0) return Complex(0.0,0.0);
  else if(x.re > 0) {
    float_tt re=sqrt(0.5*(mag+x.re));
    return Complex(re,0.5*x.im/re);
  } else {
    float_tt im=sqrt(0.5*(mag-x.re));
    if(x.im < 0) im=-im;
    return Complex(0.5*x.im/im,im);
  }
}

inline Complex polar(float_tt r, float_tt t)
{
  return Complex(r*cos(t), r*sin(t));
}

// Complex exponentiation
inline Complex pow(const Complex& z, const Complex& w)
{
  float_tt u=w.re;
  float_tt v=w.im;
  if(z == 0.0) return w == 0.0 ? 1.0 : 0.0;
  float_tt logr=0.5*log(abs2(z));
  float_tt th=arg(z);
  float_tt phi=logr*v+th*u;
  return exp(logr*u-th*v)*Complex(cos(phi),sin(phi));
}

inline Complex pow(const Complex& z, float_tt u)
{
  if(z == 0.0) return u == 0.0 ? 1.0 : 0.0;
  float_tt logr=0.5*log(abs2(z));
  float_tt theta=u*arg(z);
  return exp(logr*u)*Complex(cos(theta),sin(theta));
}

inline istream& operator >> (istream& s, Complex& y)
{
  char c;
  s >> ws >> c;
  if(c == '(') {
    s >> y.re >> c;
    if(c == ',') s >> y.im >> c;
    else y.im=0.0;
  } else {
    s.putback(c);
    s >> y.re; y.im=0.0;
  }
  return s;
}

inline ostream& operator << (ostream& s, const Complex& y)
{
  s << "(" << y.re << "," << y.im << ")";
  return s;
}

//inline bool isfinite(const Complex& z)
//{
//#ifdef _WIN32
//  return _finite(z.re) && _finite(z.im);
//#else
//  return !(isinf((double)z.re) || isnan((double)z.re) || isinf((double)z.im) || isnan((double)z.im));
//#endif
//}

#endif
