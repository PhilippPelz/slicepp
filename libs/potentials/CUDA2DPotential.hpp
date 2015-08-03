/*
 * CUDA2DPotential.hpp
 *
 *  Created on: Jul 29, 2015
 *      Author: wenxuan
 */
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <float.h>
#include "float.h"
#include <math.h>
#include <string.h>
#include <cufft.h>
#include <curand.h>
#include <cublas_v2.h>
#include "curand_kernel.h"
#include "projectedPotential.hpp"
#include "pot_base.hpp"

#ifndef CUDA2DPOTENTIAL_HPP_
#define CUDA2DPOTENTIAL_HPP_
#define iwCoord(s,i,m) if ((i) > (m)/2) {(s) = (i) - (m);} else {(s) = (i);}

#define owCoord(t,i,m) ((t) = ((i) - ((m)/2))) //m >> 1;

#define iwCoordIp(i,m) if ((i) > (m)/2) {(i) -= (m);}

#define owCoordIp(i,m) ((i) -= ((m)/2))

#define dbCoord(i1, i2, j, m1) ((i1) = ((j) % (m1))); ((i2) = ((j) - (i1)) /(m1))

#define trCoord( i1, i2, i3, j, m1, m2 ) ( (i1) = ( ( (j) % ( (m1) * (m2) ) ) % (m1) ) ); ( (i2) = ( ( ( (j) % ( (m1) * (m2) ) ) - (i1) ) / (m1) ) ); ( (i3) = ( ( (j) - (i1) - ( (m1) * (i2) ) ) / ( (m1) * (m2) ) ) )

#define sgCoord(j, i1, i2, m1) ((j) = (((i2) * (m1)) + (i1)))

#define sgCoord3D(j, i1, i2, i3, m1, m2) ( (j) = ( ( (m1)*(m2)*(i3) ) + ( (m1)*(i2) ) + (i1) ) )


namespace QSTEM {

class CUDA2DPotential: public CPotential {

public:
	CUDA2DPotential(const ConfigPtr& c,const PersistenceManagerPtr& persist);
	virtual ~CUDA2DPotential();
//	virtual void CenterAtomZ(atom& atom, float_tt &z);
//	virtual void DisplayParams();
	virtual void MakeSlices(superCellBoxPtr info);
//	virtual void Refresh();
//	virtual void ReadPotential(std::string &fileName, unsigned subSlabIdx);
//	virtual void SetStructure(StructurePtr structure);

//	inline ComplexArray2DView GetSlice(unsigned idx){return _t[boost::indices[idx][range(0,_c->Model.nx)][range(0,_c->Model.ny)]];}
//	virtual af::array GetSlice(unsigned idx);
//	virtual af::array GetSubPotential(int startx, int starty, int nx, int ny);
	void phaseGrating(struct cublasContext *cublasHandle, cufftComplex* V_d, int nAt, int nZ, std::vector<float_tt> xyzCoordFP_d, float_tt imPot, std::vector<int> Z_d, std::vector<int> Zlist, std::vector<float_tt> occ_d, cufftHandle cufftPlanBatch, int s );
protected:
//	virtual void SliceSetup();
//	virtual void AddAtomToSlices(atom& atom, float_tt atomX, float_tt atomY, float_tt atomZ)=0;
	virtual void ComputeAtomPotential(int znum){};
	virtual void SaveAtomicPotential(int znum)=0;
	void myGBSize(int *gbs, int size);
	int myGSize(int size);
	int myBSize(int size);
};
__global__ void squareAtoms_d (  cufftComplex* V, int nAt, int *Z, int Z0, float_tt *xyz, float_tt imPot, float_tt *occ, int s, int nx, int ny, int nslice, float_tt dx, float_tt dy, float_tt dz );
__global__ void divideBySinc ( cufftComplex* V, int nx, int ny, float_tt PI);
__global__ void multiplyWithProjectedPotential_d ( cufftComplex* V1, cufftComplex* V2, int nx, int ny);
__global__ void initialValues ( cuComplex* V, int size, float_tt initRe, float_tt initIm );
__device__ int mySignum_d( float x );

} /* namespace QSTEM */

#endif /* CUDA2DPOTENTIAL_HPP_ */
