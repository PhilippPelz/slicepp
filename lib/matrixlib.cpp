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

#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include "matrixlib.hpp"

#define _CRTDBG_MAP_ALLOC
#include <stdio.h>	/* ANSI C libraries */
#include <stdlib.h>
#ifdef WIN32
#if _DEBUG
#include <crtdbg.h>
#endif
#endif

#define TINY 0//1.0e-20; //A small number.

static float sqrarg = 0.0f; 
#define IMIN(a,b) ((a) < (b) ? (a) : (b))    // minimum of 2 int values
#define FMAX(a,b) ((a) > (b) ? (a) : (b))    // maximum of a and b
#define SIGN(a,b) ((b) >= 0 ? (a) : -(a))    // Magnitude of a times sign of b.
#define SQR(a) (((sqrarg=(a)) == 0.0) ? 0.0 : sqrarg*sqrarg)

namespace slicepp
{

void ludcmp(float_tt **a, int n, int *indx, float_tt *d) {
	/* Given a matrix a[1..n][1..n], this routine replaces it by the
     LU decomposition of a rowwise permutation of itself. 
     a and n are input. a is output, arranged as in equation (2.3.14) above; 
     indx[1..n] is an output vector that records the row permutation e ected
     by the partial pivoting; 
     d is output as   1 depending on whether the number of row interchanges 
     was even or odd, respectively. This routine is used in combination with 
     lubksb to solve linear equations or invert a matrix. 
	 */
	int i,imax=0,j,k;
	float_tt big,dum,sum,temp;
	/*
    vv stores the implicit scaling of each row.
	 */
	std::vector<float_tt> vv(n);

	*d=1.0; //No row interchanges yet.
	for (i=0;i<n;i++) {
		// Loop over rows to get the implicit scaling information.
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big)
				big=temp;
		if (big == 0.0) {
			printf("Singular matrix in routine ludcmp\n");
			exit(0);
		}
		//No nonzero largest element.
		vv[i]=1.0/big; // Save the scaling.
	}
	for (j=0;j<n;j++) {
		// This is the loop over columns of Crout's method.
		for (i=0;i<j;i++) {  //This is equation (2.3.12) except for i = j.
			sum=a[i][j];
			for (k=0;k<i;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		//Initialize for the search for largest pivot element.
		for (i=j;i<n;i++) {
			//This is i = j of equation (2.3.12) and i = j +1: ::N of equation (2.3.13).
			sum=a[i][j];
			for (k=0;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				//Is the  gure of merit for the pivot better than the best so far?
				big=dum; imax=i;
			}
		}
		if (j != imax) {
			//Do we need to interchange rows?
			for (k=0;k<n;k++) {
				//Yes, do so...
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d); //...and change the parity of d.
			vv[imax]=vv[j];
			// Also interchange the scale factor.
		}
		indx[j]=imax;
		if (a[j][j] == 0.0)
			a[j][j]=TINY;
		/*If the pivot element is zero the matrix is singular
      (at least to the precision of the algorithm). 
      For some applications on singular matrices, it is 
      desirable to substitute TINY for zero. 
		 */
		if (j != n-1) { // Now,  nally, divide by the pivot element.
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
	// Go back for the next column in the reduction.
	//free(vv);
	// free_vector(vv,1,n);
}

void lubksb(float_tt **a, int n, int *indx, float_tt b[]) { 
	/*
    Solves the set of n linear equations A   X = B. 
    Herea[0..n-1][0..n-1] is input, not as the matrix A but rather 
    as its LU decomposition, determined by the routine ludcmp. 
    indx[0..n-1] is input as the permutation vector returned by ludcmp.
    b[0..n-1] is input as the right-hand side vector B, and returns 
    with the solution vector X. a, n, andindx are not modi ed by 
    this routine and can be left in place for successive calls with 
    di erent right-hand sides b. This routine takes into account 
    the possibility that b will begin with many zero elements, 
    so it is e cient for use in matrix inversion. 
	 */
	int i,ii=-1,ip,j;
	float_tt sum;
	for (i=0;i<n;i++) {
		/* When ii is set to a positive value, it will become the index
       of the  rst nonvanishing element of b. Wenow do the forward 
       substitution, equation (2.3.6). The only new wrinkle is to 
       unscramble the permutation as we go. 
		 */
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii+1)
			for (j=ii;j<=i-1;j++)
				sum -= a[i][j]*b[j];
		else
			if (sum) ii=i;
		/*
      A nonzero element was encountered, so from now on we will have 
      to do the sums in the loop above. 
		 */
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) { //Now we do the backsubstitution, equation (2.3.7).
		sum=b[i];
		for (j=i+1;j<n;j++)
			sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i]; // Store a component of the solution vector X.
	}
	//All done!
}

/******************************************************************************
 Routine:   det_3x3
 Input:     m - matrix (3x3) address
 Output:    returns the determinant of 'm'
 ******************************************************************************/
float_tt det_3x3 ( const FloatArray2D& mat)
{
	float_tt det;
	const float_tt* mat1 = mat.data();
	det = mat1[0] * (mat1[4] * mat1[8] - mat1[7] * mat1[5])
				- mat1[1] * (mat1[3] * mat1[8] - mat1[6] * mat1[5])
				+ mat1[2] * (mat1[3] * mat1[7] - mat1[6] * mat1[4]);

	return det;
}

/******************************************************************************
 Routine:   trans_3x3
 Input:     Mstarget,Msource - matrix (3x3)
 ******************************************************************************/
void trans_3x3 (float_tt *Mt, const float_tt *Ms)
{
	int i,j;
	for (i=0;i<2;i++) for (j=0;j<3;j++) Mt[i*3+j] = Ms[j*3+i];
	return;
}


/******************************************************************************
 Routine:   inverse_3x3
 Input:     m - matrix (3x3) address
 Output:    returns the inverse matrix of 'm'
 ******************************************************************************/
void inverse_3x3(FloatArray2D& res, const FloatArray2D& a){

	float_tt det = det_3x3 (a);

	// printf("det: %g\n",det);

	if (fabs (det) < 0.0005f)
	{
		float_tt* res1 = res.data();
		res1[0] = 1.0f;
		res1[1] = 0.0f;
		res1[2] = 0.0f;
		res1[3] = 0.0f;
		res1[4] = 1.0f;
		res1[5] = 0.0f;
		res1[6] = 0.0f;
		res1[7] = 0.0f;
		res1[8] = 1.0f;
		return;
	}

	det = 1.0f / det;
	float_tt* res1 = res.data();
	const float_tt* a1 = a.data();
	res1[0] =  (a1[4]*a1[8] - a1[5]*a1[7]) * det;
	res1[1] = -(a1[1]*a1[8] - a1[7]*a1[2]) * det;
	res1[2] =  (a1[1]*a1[5] - a1[4]*a1[2]) * det;

	res1[3] = -(a1[3]*a1[8] - a1[5]*a1[6]) * det;
	res1[4] =  (a1[0]*a1[8] - a1[6]*a1[2]) * det;
	res1[5] = -(a1[0]*a1[5] - a1[3]*a1[2]) * det;

	res1[6] =  (a1[3]*a1[7] - a1[6]*a1[4]) * det;
	res1[7] = -(a1[0]*a1[7] - a1[6]*a1[1]) * det;
	res1[8] =  (a1[0]*a1[4] - a1[1]*a1[3]) * det;
}

/* This routine calculates inv(M)*M for any square N x M matrix (M >= N).
 * This function is useful for determining, whether a vector is represented by M or not,
 * even if M is not orthogonal. For any rowvector vec
 * vec == invMM*vec, if vec is represented by M, or vec!=invMM*vec if not.
 * the matrix M will be preserved
 * The invMM matrix will be of size M X M
 */
//float_tt** invMM(float_tt **Mmatrix, int N, int M) {
//	static int Mold = 0; // Nold = 0,
//	static float_tt **invMMmatrix = NULL;
//	// int i;
//
//	if (N > M) return NULL;
//
//	if (Mold < M) {
//		Mold = M;
//		if (invMMmatrix != NULL) {
//			free(invMMmatrix[0]); free(invMMmatrix);
//		}
//		invMMmatrix = float2D(M,M,"invMMmatrix");
//	}
//
//	return invMMmatrix;
//}


void svdcmp1(float_tt **a, int m, int n, float_tt w[], float_tt **v) {
	/* Given a matrix a[1..m][1..n], this routine computes its singular value decomposition,
     A = U   W   V T . 
     The matrix U replaces a on output. The diagonal matrix of singular values 
     W is out- put as a vector w[1..n]. The matrix V (not the transpose V T ) is output
     as v[1..n][1..n].
	 */
	// uses:  float_tt pythag(float_tt a, float_tt b);
	int flag,i,its,j,jj,k,l=0,nm=0;
	float_tt anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	// rv1=vector(1,n); we will just alloacte memory for n+1 elements
	rv1 = (float_tt *)malloc((n+1)*sizeof(float_tt));
	g=scale=anorm=0.0; // Householder reduction to bidiagonal form.
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++)
				scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if(i<= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {

				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++)
						s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++)
						a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	} 
	for (i=n;i>=1;i--) { // Accumulation of right-hand transformations.
		if(i < n) {
			if (g) {
				for (j=l;j<=n;j++) // Float_Tt division to avoid possible under ow.
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	} 
	for (i=IMIN(m,n);i>=1;i--) { // Accumulation of left-hand transformations.
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		}
		else for (j=i;j<=m;j++) a[j][i]=0.0; ++a[i][i];
	} 
	for (k=n;k>=1;k--) {
		// Diagonalization of the bidiagonal form: Loop over singular values,
		// and over allowed iterations.
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) { // Test for splitting.
				nm=l-1; // Note that rv1[1] is always zero.
				if ((fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0; // Cancellation of rv1[l], if l>1.
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((fabs(f)+anorm) == anorm)
						break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) { //Convergence.
				if (z < 0.0) { // Singular value is made nonnegative.
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) {
				printf("no convergence in 30 svdcmp iterations");
				exit(0);
			}
			x=w[l]; // Shiftfrom bottom 2-by-2minor.
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			// Next QR transformation:
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				// Rotation can be arbitrary if z = 0.
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	} 
	//  free_vector(rv1,1,n);
	free(rv1);
}


float_tt pythag(float_tt a, float_tt b) {
	/* Computes (a 2 + b 2 ) 1=2 without destructive under ow or over ow.
	 */
	float_tt absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb)
		return absa*sqrt(1.0+SQR(absb/absa));
	else
		return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


/* This function will calculate the 3-dimensional vector cross product 
 * c = [a x b]
 */
void crossProduct(const float_tt *a, const float_tt *b, float_tt *c) {
	c[0] = a[1]*b[2]-a[2]*b[1];
	c[1] = a[2]*b[0]-a[0]*b[2];
	c[2] = a[0]*b[1]-a[1]*b[0];
	return;
}

float_tt dotProduct(const float_tt *a, const float_tt *b) {
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

/* vector difference c = a - b 
 * If revFlag < 0 we will treat the first parameter's (a) coordinates
 * in reversed order, i.e. a[0]=z, a[1]=y, a[2]=x.
 */

void vectDiff_f(const float_tt *a,const float_tt *b, float_tt *c,int revFlag) {
	if (revFlag > 0) {
		c[0] = a[0]-b[0];
		c[1] = a[1]-b[1];
		c[2] = a[2]-b[2];
	}
	else {
		c[0] = a[2]-b[0];
		c[1] = a[1]-b[1];
		c[2] = a[0]-b[2];
	}
}

#define HALF_PI    (0.5f * PI)
#define TWO_PI     (2.0f * PI)
#define TWO_PI_INV (1.0f / TWO_PI)

inline float_tt Hill(float_tt x)
{
	const float_tt a0 = 1.0f;
	const float_tt a2 = 2.0f / PI - 12.0f / (PI * PI);
	const float_tt a3 = 16.0f / (PI * PI * PI) - 4.0f / (PI * PI);
	const float_tt xx = x * x;
	const float_tt xxx = xx * x;

	return a0 + a2 * xx + a3 * xxx;
}

float_tt fsin(float_tt x)
{
	// wrap x within [0, TWO_PI)
	const float_tt a = x * TWO_PI_INV;
	x -= static_cast<int>(a) * TWO_PI;
	if (x < 0.0f)
		x += TWO_PI;

	// 4 pieces of hills
	if (x < HALF_PI)
		return Hill(HALF_PI - x);
	else if (x < PI)
		return Hill(x - HALF_PI);
	else if (x < 3.0f * HALF_PI)
		return -Hill(3.0f * HALF_PI - x);
	else
		return -Hill(x - 3.0f * HALF_PI);
}
float_tt fcos(float_tt x)
{
  return fsin(x + HALF_PI);
}



void showMatrix(FloatArray2D M,int Nx, int Ny,const char *name) {
	int i,j;

	printf("%s:\n",name);
	for (i=0;i<Nx;i++) {
		for (j=0;j<Ny;j++) printf("%10f ",M[i][j]);
		printf("\n");
	}
}

/* This function will determine whether the point (point) lies
 * above or below the plane (p)
 * The plane defining vectors vect1, vect2, norm, and point will be used
 * revFlag indicates whether the parameter point is given in 
 * forward(1) (x,y,z), or reversed(-1) (z,y,x) order.
 * This is important for using the reversed order in the atom struct.
 */
float_tt findLambda(plane *p,const float_tt *point, int revFlag) {
	FloatArray2D M(boost::extents[3][3]),
			Minv(boost::extents[3][3]);
	//static float_tt *diff=NULL;
	std::vector<float_tt> diff(3);
	float_tt lambda; /* dummy variable */

	/*
  printf("hello x=(%g %g %g)\n",p->pointX,p->pointY,p->pointZ);
  printf("hello p->norm=(%g %g %g), (%d %d %d)\n",p->normX,p->normY,p->normZ,
	 (int)M[0],(int)M[1],(int)M[2]);
	 */

	//	M[0][0] = -((*p).normX);
	//	M[0][3] = -((*p).normY);
	//	M[0][6] = -((*p).normZ);
	//	M[0][1] = p->vect1X;
	//	M[1][1] = p->vect1Y;
	//	M[2][1] = p->vect1Z;
	//	M[0][2] = p->vect2X;
	//	M[1][2] = p->vect2Y;
	//	M[2][2] = p->vect2Z;
	M[0][0] = -((*p).n[0]);
	M[1][0] = -((*p).n[1]);
	M[2][0] = -((*p).n[2]);
	M[0][1] = p->v1[0];
	M[1][1] = p->v1[1];
	M[2][1] = p->v1[2];
	M[0][2] = p->v2[0];
	M[1][2] = p->v2[1];
	M[2][2] = p->v2[2];

	//	BOOST_LOG_TRIVIAL(info)<< format("(%3.3f,%3.3f,%3.3f) - (%3.3f,%3.3f,%3.3f)") % point[0]% point[1]% point[2]% p->p[0]% p->p[1]% p->p[2];

	vectDiff_f(point,p->p.mem,&diff[0],revFlag);

	//	BOOST_LOG_TRIVIAL(info)<< format("diff=(%3.3f,%3.3f,%3.3f)") % diff[0]% diff[1]% diff[2];

	inverse_3x3 (Minv,M);
	//	showMatrix(M,3,3,"M");
	//	showMatrix(Minv,3,3,"Minv");
	//	showMatrix(&diff[0],1,3,"diff");
	// ludcmp(M,3,index,&d);
	// lubksb(M,3,index,point);

	lambda = dotProduct( Minv.data(),&diff[0]);
	// printf("lambda: %g\n",lambda);
	return lambda;
}




/* This function will perform a 3D rotation about the angles
 * phi_x,phi_y,phi_z (in that sequence and about the respective orthogonal axes)
 * of the vector vectIn, and store the result in vectOut.
 * The angles will be written in a matrix, which will be saved for the next rotation,
 * in case the angles don't change.
 * The same vector can be specified as input, as well as output vector.
 */
void rotateVect(float_tt *vectIn,float_tt *vectOut, float_tt phi_x, float_tt phi_y, float_tt phi_z) {
	static FloatArray2D Mrot(boost::extents[3][3]);
	std::vector<float_tt>vectOutTemp(3);
	static float_tt sphi_x=0, sphi_y=0, sphi_z=0;
	if ((phi_x!=sphi_x) || (phi_y!=sphi_y) || (phi_z!=sphi_z)) {
		Mrot[0][0] = cos(phi_z)*cos(phi_y);
		Mrot[0][1] = cos(phi_z)*sin(phi_y)*sin(phi_x)-sin(phi_z)*cos(phi_x);
		Mrot[0][2] = cos(phi_z)*sin(phi_y)*cos(phi_x)+sin(phi_z)*sin(phi_x);

		Mrot[1][0] = sin(phi_z)*cos(phi_y);
		Mrot[1][1] = sin(phi_z)*sin(phi_y)*sin(phi_x)+cos(phi_z)*cos(phi_x);
		Mrot[1][2] = sin(phi_z)*sin(phi_y)*cos(phi_x)-cos(phi_z)*sin(phi_x);

		Mrot[2][0] = -sin(phi_y);
		Mrot[2][1] = cos(phi_y)*sin(phi_x);
		Mrot[2][2] = cos(phi_y)*cos(phi_x);

		sphi_x = phi_x;
		sphi_y = phi_y;
		sphi_z = phi_z;
	}
	vectOutTemp[0] = Mrot[0][0]*vectIn[0]+Mrot[0][1]*vectIn[1]+Mrot[0][2]*vectIn[2];
	vectOutTemp[1] = Mrot[1][0]*vectIn[0]+Mrot[1][1]*vectIn[1]+Mrot[1][2]*vectIn[2];
	vectOutTemp[2] = Mrot[2][0]*vectIn[0]+Mrot[2][1]*vectIn[1]+Mrot[2][2]*vectIn[2];
	memcpy(vectOut,&vectOutTemp[0],3*sizeof(float_tt));

	return;
}
static float_tt sphi_x=0, sphi_y=0, sphi_z=0;
void rotateMatrix(const FloatArray2D& matrixIn, FloatArray2D& matrixOut, float_tt phi_x, float_tt phi_y, float_tt phi_z) {
	int i,j,k;
	static FloatArray2D Mrot(boost::extents[3][3]);
	std::fill(Mrot.data(), Mrot.data() + Mrot.size(), 0);
	Mrot[0][0] = Mrot[1][1] = Mrot[2][2] = 1;
	std::vector<float_tt>matrixOutTemp(9);

	if ((phi_x!=sphi_x) || (phi_y!=sphi_y) || (phi_z!=sphi_z)) {
		Mrot[0][0] = cos(phi_z)*cos(phi_y);
		Mrot[0][1] = cos(phi_z)*sin(phi_y)*sin(phi_x)-sin(phi_z)*cos(phi_x);
		Mrot[0][2] = cos(phi_z)*sin(phi_y)*cos(phi_x)+sin(phi_z)*sin(phi_x);

		Mrot[1][0] = sin(phi_z)*cos(phi_y);
		Mrot[1][1] = sin(phi_z)*sin(phi_y)*sin(phi_x)+cos(phi_z)*cos(phi_x);
		Mrot[1][2] = sin(phi_z)*sin(phi_y)*cos(phi_x)-cos(phi_z)*sin(phi_x);

		Mrot[2][0] = -sin(phi_y);
		Mrot[2][1] = cos(phi_y)*sin(phi_x);
		Mrot[2][2] = cos(phi_y)*cos(phi_x);

		sphi_x = phi_x;
		sphi_y = phi_y;
		sphi_z = phi_z;
	}
	memset(&matrixOutTemp[0],0,9*sizeof(float_tt));
	for (i=0;i<3;i++) for (j=0;j<3;j++) for (k=0;k<3;k++) {
		matrixOutTemp[i*3+j] += Mrot[i][k]*matrixIn[k][j];
	}
	memcpy(matrixOut.data(),&matrixOutTemp[0],9*sizeof(float_tt));

	return;
}

/* This function will create a basis vector set for the axes of the unit cell
 * which are defined in lengths and angles.
 * The rotationall degree of freedom will be fixed by the fact that 
 * vax = (ax,0,0)
 * vby = (by*cos(gamma),by*sin(gamma),0)
 * vcz = (c*cos(beta),(c*cos(alpha)-cos(beta)*cos(gamma))/sin(gamma),
 * c*sqrt(1-cos(alpha)^2-cos(beta)^2-cos(gamma)^2+2*cos(alpha)*cos(beta)*cos(gamma))/sin(gamma))
 * I assume that 
 * alpha = angle between by and cz, 
 * beta = angle between ax and cz, and
 * gamma = angle between by and ax.
 */
// #define SQR(x) ((x)*(x))
void makeCellVect(grainBox& grain, std::vector<float_tt>& vax, std::vector<float_tt>& vby, std::vector<float_tt>& vcz) {

	if ((grain.alpha == 90) && (grain.beta == 90) && (grain.gamma == 90)) {
		printf("Orthogonal unit cell\n");
		vax[0] = grain.ax; vax[1] = 0; vax[2] = 0;
		vby[0] = 0; vby[1] = grain.by;  vby[2] = 0;
		vcz[0] = 0; vcz[1] = 0; vcz[2] = grain.cz;
	}

	else {

		vax[0] = grain.ax; vax[1] = 0; vax[2] = 0;

		vby[0] = grain.by*cos(grain.gamma*PI180);
		vby[1] = grain.by*sin(grain.gamma*PI180);
		vby[2] = 0;

		vcz[0] = grain.cz*cos(grain.beta*PI180);
		vcz[1] = grain.cz*(cos(grain.alpha*PI180)-cos(grain.beta*PI180)*
				cos(grain.gamma*PI180))/sin(grain.gamma*PI180);
		if ((fabs(grain.alpha-90.0) < 1e-4) && (fabs(grain.beta-90.0) < 1e-4))
			vcz[2] = grain.cz;
		else  {// the general function still seems to have a bug in it.
			// printf("alpha: %g, beta: %g\n",grain.alpha-90.0,grain.beta-90.0);
			vcz[2] = grain.cz*(sqrt(1-SQR(cos(grain.alpha*PI180))-SQR(cos(grain.beta*PI180))
					+2*cos(grain.alpha*PI180)*cos(grain.beta*PI180)*
					cos(grain.gamma*PI180))/sin(grain.gamma*PI180));
		}
	}
}

float_tt vectLength(float_tt *vect) {
	return sqrt(vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]);
}

/* c = a*b */
void matrixProduct(const FloatArray2D& a,int Nxa, int Nya,const FloatArray2D& b,int Nxb, int Nyb, FloatArray2D& c) {
	int i,j,k;

	if (Nya != Nxb) {
		printf("multiplyMatrix: Inner Matrix dimensions do not agree!\n");
		return;
	}
	for (i=0;i<Nxa;i++)
		for (j=0;j<Nyb;j++) {
			c[i][j] = 0.0;
			for (k=0;k<Nya;k++)
				c[i][j] += a[i][k]*b[k][j];
		}
}


/* c = a*b */
void matrixProductInt(float_tt **a,int Nxa, int Nya, int **b,int Nxb, int Nyb, float_tt **c) {
	int i,j,k;

	if (Nya != Nxb) {
		printf("multiplyMatrix: Inner Matrix dimensions do not agree!\n");
		return;
	}

	for (i=0;i<Nxa;i++) for (j=0;j<Nyb;j++) {
		c[0][i*Nyb+j] = 0.0;
		for (k=0;k<Nya;k++) c[0][i*Nyb+j] += a[0][i*Nya+k] * (float_tt)(b[0][k*Nyb+j]);
	}
}

} // end namespace slicepp
