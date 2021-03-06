/*
QSTEM - image simulation for TEM/STEM/CBED
    Copyright (C) 2000-2010  Christoph Koch
	Copyright (C) 2010-2013  Christoph Koch, Michael Sarahan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MATRIXLIB_H
#define MATRIXLIB_H

#include "stemtypes_fftw3.hpp"

#include <boost/format.hpp>
#include <boost/log/trivial.hpp>
using boost::format;

#define PI 3.1415926535898
#define PI180 1.7453292519943e-2

namespace slicepp
{

void ludcmp(float_tt **a, int n, int *indx, float_tt *d);
void lubksb(float_tt **a, int n, int *indx, float_tt b[]);
float_tt det_3x3 ( const FloatArray2D& a);
void inverse_3x3 (FloatArray2D& res, const FloatArray2D& a);
void trans_3x3 (float_tt *Mt, const float_tt *Ms);

float_tt fsin(float_tt x);
float_tt fcos(float_tt x);

// svdcmp1 uses the NR unit-offset vectors :-(
void svdcmp1(float_tt **a, int m, int n, float_tt w[], float_tt **v);
float_tt pythag(float_tt a, float_tt b);

/* vector functions:
 */
void crossProduct(const float_tt *a, const float_tt *b, float_tt *c);
float_tt dotProduct(const float_tt *a, const float_tt *b);
float_tt findLambda(plane *p,const float_tt *point, int revFlag);
void showMatrix(FloatArray2D M,int Nx, int Ny,const char *name);
void vectDiff_f(float_tt *a, float_tt *b, float_tt *c,int revFlag);
float_tt vectLength(float_tt *vect);
void makeCellVect(grainBox& grain, std::vector<float_tt>& vax, std::vector<float_tt>& vby, std::vector<float_tt>& vcz);
void rotateVect(float_tt *vectIn,float_tt *vectOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);
void rotateMatrix(const FloatArray2D& matrixIn, FloatArray2D& matrixOut, float_tt phi_x, float_tt phi_y, float_tt phi_z);

/* |vect| */
float_tt vectLength(float_tt *vect);

/* c = a*b */
void matrixProduct(const FloatArray2D& a,int Nxa, int Nya,const FloatArray2D& b,int Nxb, int Nyb, FloatArray2D& c);
void matrixProductInt(float_tt **a,int Nxa, int Nya, int **b,int Nxb, int Nyb, float_tt **c);

} // end namespace slicepp

#endif
