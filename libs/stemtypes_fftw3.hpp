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

#ifndef STEMTYPES_H
#define STEMTYPES_H

#define FLOAT_PRECISION 1
#define DBL_STD_MAX 1.0e+308

#if(FLOAT_PRECISION == 1)

typedef float float_tt;
#define armarowvec arma::frowvec
#define armavec arma::fvec
#define armamat arma::fmat
#define afcfloat af::cfloat
#define STD_MAX 1.0e+38
#define REAL_MIN FLT_MIN
#define REAL_MAX FLT_MAX
#define REAL_EPSILON FLT_EPSILON
#define REAL_DIG FLT_DIG

#else  // FLOAT_PRECISION

typedef double float_tt;
#define armarowvec arma::rowvec
#define armavec arma::vec
#define armamat arma::mat
#define afcfloat af::cdouble
#define STD_MAX DBL_STD_MAX
#define REAL_MIN DBL_MIN
#define REAL_MAX DBL_MAX
#define REAL_EPSILON DBL_EPSILON
#define REAL_DIG DBL_DIG

#endif  // FLOAT_PRECISION


#include "Complex.hpp"
#include "fftw3.h"
#include <boost/shared_ptr.hpp>
#include "boost/multi_array.hpp"
#include "elTable.hpp"
#include "Complex.hpp"
#include <complex>
#include <cmath>
#include <iomanip>
#include <map>
#include <list>
#include <armadillo>
#include <boost/rational.hpp>
using namespace std;


namespace QSTEM
{
typedef boost::rational<unsigned> ratio;
typedef Complex complex_tt;
typedef std::vector<float_tt> RealVector;
typedef std::vector<complex_tt> ComplexVector;

typedef boost::multi_array<complex_tt, 3> ComplexArray3D;
typedef boost::multi_array<complex_tt, 2> ComplexArray2D;
typedef boost::multi_array<float_tt, 3> FloatArray3D;
typedef boost::multi_array<float_tt, 2> FloatArray2D;
typedef ComplexArray3D::index ComplexArray3DIndex;
typedef boost::multi_array_types::index_range range;
typedef ComplexArray3D::array_view<2>::type ComplexArray2DView;
typedef ComplexArray3D::array_view<3>::type ComplexArray3DView;
typedef FloatArray2D::array_view<2>::type FloatArray2DView;
typedef FloatArray3D::array_view<3>::type FloatArray3DView;
typedef boost::multi_array_ref<complex_tt,3> ComplexArray3DPtr;
typedef boost::multi_array_ref<complex_tt,2> ComplexArray2DPtr;

#define SHOW_SINGLE_POTENTIAL 0

#define BW (2.0F/3.0F)	/* bandwidth limit */

// Helpful mathematical constants
#define RAD2DEG 57.2958
#define SQRT_2 1.4142135
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

////////////////////////////////////////////////////////////////////////
// Define physical constants
////////////////////////////////////////////////////////////////////////
#define ELECTRON_CHARGE (1.6021773e-19)
#define PICO_AMPERE (1e-12/ELECTRON_CHARGE)
#define MILLISEC_PICOAMP (1e-3*PICO_AMPERE)
#define EM_SIGMA 7288400

////////////////////////////////////////////////////////////////

const float_tt PI = 2*acos((float_tt)0.0);

// FFTW constants
const int k_fftMeasureFlag = FFTW_ESTIMATE;

struct atom {
public:
	armavec r = armavec(3);
	// float dx,dy,dz;  // thermal displacements
	float_tt dw;      // Debye-Waller factor
	float_tt occ;     // occupancy
	float_tt q;       // charge
	unsigned Znum;
	float_tt mass;
	atom();
	atom(std::vector<atom>::iterator a){
		dw=a->dw;
		occ = a->occ;
		q =a->q;
		Znum = a->Znum;
		mass =a->mass;
		r = armavec(a->r);
	}
	atom(const atom& a){
		dw=a.dw;
		occ = a.occ;
		q =a.q;
		Znum = a.Znum;
		mass =a.mass;
		r = armavec(a.r);
	}
	// constructor (used for testing)
	atom(float_tt _mass, std::string _symbol,
			float_tt _x, float_tt _y, float_tt _z,
			float_tt _dw, float_tt _occ, float_tt _charge);

};

inline atom::atom()
: dw(0)
, occ(1)
, q(0)
, Znum(0)
, mass(0)
, r(armavec(3))
{
}

inline atom::atom(float_tt _mass, std::string _symbol,
		float_tt _x, float_tt _y, float_tt _z,
		float_tt _dw, float_tt _occ, float_tt _charge)
: dw(_dw)
, occ(_occ)
, q(_charge)
, mass(_mass)
{
	Znum = getZNumber(_symbol.c_str());
}

/* Planes will be defined by the standard equation for a plane, i.e.
 * a point (point) and 2 vectors (vect1, vect2)
 */
struct plane {
	armavec v1 = armavec(3);
	armavec v2 = armavec(3);
	armavec n = armavec(3);
	armavec p = armavec(3);
};

struct grainBox {
	int amorphFlag;
	float_tt density;
	float_tt rmin;
	float_tt rFactor;
	/* density, atomic distance, reduced atomic distance
	 * for amorphous material.  Red. r is for making a hex.
	 * closed packed structure, which will fill all space,
	 * but will only be sparsely filled, and later relaxed.
	 */
	string name;
	std::vector<atom> unitCell = std::vector<atom>();
	float_tt ax,by,cz;
	float_tt alpha=0, beta=0, gamma=0;
	float_tt tiltx=0,tilty=0,tiltz=0;
	float_tt shiftx=0,shifty=0,shiftz=0;
	std::vector<plane> planes;   /* pointer to array of bounding planes */
	float_tt sphereRadius, sphereX,sphereY,sphereZ; /* defines a sphere instead of a grain with straight edges */
	armamat M;
} ;

struct superCellBox {
	float_tt cmx,cmy,cmz;  /* fractional center of mass coordinates */
	float_tt ax,by,cz;
	std::vector<atom> atoms; /* contains all the atoms within the super cell */
	std::vector<int> uniqueatoms;
	std::vector<int> znums;
	std::vector<float_tt> xyzPos;
	std::vector<float_tt> occupancy;
	std::map<int,std::vector<int>> zNums;
};



struct atomBox {
	bool used;   /* indicate here whether this atom is used in the
		 particular problem */
	unsigned nx,ny,nz;
	float_tt dx,dy,dz;
	float_tt B;
	FloatArray3D rpotential;
	ComplexArray3D potential;
} ;

typedef boost::shared_ptr<plane> planePtr;
typedef boost::shared_ptr<grainBox> grainBoxPtr;
typedef boost::shared_ptr<superCellBox> superCellBoxPtr;
typedef boost::shared_ptr<atomBox> atomBoxPtr;
typedef boost::shared_ptr<atom> atomPtr;


static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 80)
{
	//	if ( (x != n) && (x % (n/100+1) != 0) ) return;

	float ratio  =  x/(float)n;
	int   c      =  ratio * w;

	std::cout << std::setw(3) << (int)(ratio*100) << "% [";
	for (int x=0; x<c; x++) std::cout << "=";
	for (int x=c; x<w; x++) std::cout << " ";
	//    std::cout << "]\r" << std::flush;
	std::cout << "]" << std::endl << std::flush;
}

// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
#define DLL_IMPORT __declspec(dllimport)
#define DLL_EXPORT __declspec(dllexport)
#define DLL_LOCAL
#else
#if __GNUC__ >= 4
#define DLL_IMPORT __attribute__ ((visibility ("default")))
#define DLL_EXPORT __attribute__ ((visibility ("default")))
#define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#else
#define DLL_IMPORT
#define DLL_EXPORT
#define DLL_LOCAL
#endif
#endif

}

#endif // STEMTYPES_H
