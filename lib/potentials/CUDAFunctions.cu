/*
 * CUDAFunctions.cpp
 *
 *  Created on: Aug 3, 2015
 *      Author: wenxuan
 */

#include <stdlib.h>
#include <iostream>
#include "CUDAFunctions.hpp"
#include "projectedPotential.hpp"

namespace QSTEM {

CUDAFunctions::CUDAFunctions(superCellBoxPtr info, cModelConfPtr mc) {
	_info = info;
	_mc = mc;

	af::timer time = af::timer::start();

	// create the needed arrays on the gpu, with a few workarounds to get arrayfire working with cuda kernels
	xyzPos = af::array(info->xyzPos.size(), (float_tt *) info->xyzPos.data(), afHost);
	occupancy = af::array(info->occupancy.size(), (float_tt *) info->occupancy.data(), afHost);
	znums = af::array(info->znums.size(), (int *) info->znums.data(), afHost);

	xyzPos *= 1;
	occupancy *= 1;
	znums *= 1;

	xyzPos_d = xyzPos.device<float_tt>();
	occupancy_d = occupancy.device<float_tt>();
	znums_d = znums.device<int>();

	af::sync();

	slicePixels = _mc->nx * _mc->ny;
	gS = myGSize(slicePixels);
	bS = myBSize(slicePixels * _mc->nSlices);
	gS2D = myGSize(slicePixels);

    af_id = af::getDevice();
    afcu_get_stream(&af_stream,af_id);
//     = afcu::getStream(af_id);

//	BOOST_LOG_TRIVIAL(info)<< boost::format( "%g msec used copying data to gpu") % (af::timer::stop(time)*1000);
}

void CUDAFunctions::GetPhaseGrating(cufftComplex* V_slice, int slice, std::map<int, af::array> & atomPot) {
	int nAtom = _info->atoms.size();

	SetComplex( V_slice, 0.f, 0.f);
	SetComplex( _V_accum_ptr, 0.f, 0.f);
	for (int& Z : _info->uniqueZ) {
		SetComplex( _V_elem_ptr, 0.f, 0.f);
		SetComplex( _V_atom_ptr, 0.f, 0.f);

		GetAtomDeltaFunctions(_V_elem_ptr, Z, slice);

		cuda_assert(cudaDeviceSynchronize());

		_V_elem.unlock();
		_v_accum.unlock();

		_V_elem = af::fft(_V_elem);
		_V_elem *= atomPot[Z];
		_V_elem = af::ifft(_V_elem);

		_V_elem *= slicePixels;
		_v_accum = _V_elem + _v_accum;
		af::sync();

		PotentialToTransmission(V_slice, _V_accum_ptr);
	}
}
void CUDAFunctions::PotentialToTransmission(cufftComplex* pot, cufftComplex* trans){
	potential2Transmission<<< gS, bS , 0, af_stream>>> (pot, trans, slicePixels);
}
void CUDAFunctions::SetComplex(cufftComplex* a, float real, float imag){
	int slicePixels = _mc->nx * _mc->ny;
	const int gS = myGSize(slicePixels);
	const int bS = myBSize(slicePixels * _mc->nSlices);
	initialValues<<< gS * 2, bS , 0, af_stream>>> ( a, slicePixels, 0.f, 0.f);
}
void CUDAFunctions::GetAtomicPotential(cufftComplex* V, int Z) {
	int slicePixels = _mc->nx * _mc->ny;
	const int bS = myBSize(slicePixels * _mc->nSlices);
	const int gS2D = myGSize(slicePixels);
	createAtomicPotential<<< gS2D, bS , 0, af_stream>>> ( V, Z, _mc->nx, _mc->ny, _mc->dx, _mc->dy);
}
void CUDAFunctions::GetSincAtomicPotential(cufftComplex* V, int Z) {
	int slicePixels = _mc->nx * _mc->ny;
	const int bS = myBSize(slicePixels * _mc->nSlices);
	const int gS2D = myGSize(slicePixels);
	createAtomicPotential<<< gS2D, bS , 0, af_stream>>> ( V, Z, _mc->nx, _mc->ny, _mc->dx, _mc->dy);
	divideBySinc<<< gS2D, bS , 0, af_stream>>> ( V, _mc->nx, _mc->ny, PI);
}
void CUDAFunctions::GetAtomDeltaFunctions(cufftComplex* V, int Z, int slice) {
	int nAtom = _info->atoms.size();
	putAtomDeltas<<< myGSize( nAtom ), myBSize( nAtom ) , 0, af_stream>>> ( V, nAtom, znums_d, Z, xyzPos_d, _mc->imPot, occupancy_d, slice, _mc->nx, _mc->ny, _mc->nSlices, _mc->dx, _mc->dy, _mc->dz);
}
void CUDAFunctions::printPotArray(cufftComplex* V_d, int nx, int ny) {
	cufftComplex *V_host;
	V_host = (cufftComplex *) malloc(nx * ny * sizeof(cufftComplex));
	cuda_assert(cudaMemcpy(V_host, V_d, nx * ny * sizeof(cufftComplex), cudaMemcpyDeviceToHost));
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			float x = V_host[i * nx + j].x;
			float y = V_host[i * nx + j].y;
			if (x > 1e-10 || y > 1e-10)
				cout << "(" << x << ", " << y << ")" << endl;
		}
	}
	free(V_host);
}

void CUDAFunctions::printFloatArray(float_tt* f, int nx, int ny, int offset) {
	char *f_host;
	f_host = (char *) malloc(nx * ny * sizeof(float_tt));
	cuda_assert(cudaMemcpy((void * )f_host, (void * )f, nx * ny * sizeof(float_tt), cudaMemcpyDeviceToHost));
	f_host += offset;
	float_tt *ff = (float_tt *) f_host;
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			cout << ff[i * ny + j] << " ";
		}
		cout << endl;
	}
	free(f_host);
}

void CUDAFunctions::printIntArray(int* p, int size) {
	int *f_host;
	f_host = (int *) malloc(size * sizeof(int));
	cuda_assert(cudaMemcpy(f_host, p, size * sizeof(int), cudaMemcpyDeviceToHost));
	for (int i = 0; i < size; i++) {
		cout << f_host[i] << endl;
	}
	free(f_host);
}
void CUDAFunctions::initPotArrays(int slicePixels) {
	_V_elem = af::array(slicePixels, c32);
	_V_atom = af::array(slicePixels, c32);
	_v_accum = af::array(slicePixels, c32);
	af::sync();
	_V_elem_ptr = (cufftComplex *) _V_elem.device<af::af_cfloat>();
	_V_atom_ptr = (cufftComplex *) _V_atom.device<af::af_cfloat>();
	_V_accum_ptr = (cufftComplex *) _v_accum.device<af::af_cfloat>();
}
void CUDAFunctions::unlockArrays() {
	_V_elem.unlock();
	_V_atom.unlock();
	_v_accum.unlock();
}
__global__ void initialValues(cuComplex* V, int size, float_tt initRe, float_tt initIm) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < 2 * size) {
		if ((i % 2) == 0) {
			V[i / 2].x = initRe;
		}

		else {
			V[i / 2].y = initIm;
		}
	}
}

__global__ void putAtomDeltas(cufftComplex* V, int nAt, int *Z, int Z0, float_tt *xyz, float_tt imPot, float_tt *occ, int s, int nx, int ny,
		int nSlices, float_tt dx, float_tt dy, float_tt dz) //double linear interpolation
		{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nAt) {
		if (Z[i] == Z0) {
			const int m1 = nx;
			const int m2 = ny;
			float_tt x1, x2;
			x1 = xyz[i * 3 + 0] / dx - 0.5f;
			x2 = xyz[i * 3 + 1] / dy - 0.5f;
			int i3 = (int) (roundf(xyz[i * 3 + 2] / dz - 0.5f));

			if (((x1 > 1.f) && (x1 < ((float_tt) (m1 - 2)))) && ((x2 > 1.f) && (x2 < ((float_tt) (m2 - 2)))) && ((i3 > s - 1) && (i3 <= s))) {
				int i1 = (int) roundf(x1);
				int i2 = (int) roundf(x2);
				int j;
				float_tt r1 = x1 - ((float_tt) i1);
				float_tt r2 = x2 - ((float_tt) i2);
				float_tt temp;

				sgCoord(j, i1, i2, m1);
				temp = (1 - fabsf(r1)) * (1 - fabsf(r2)) * occ[i];
				//V[j].x += temp;
				//V[j].y += temp * imPot;
				atomicAdd(&(V[j].x), temp);
				atomicAdd(&(V[j].y), temp * imPot);

				i2 += sign(r2);
				sgCoord(j, i1, i2, m1);
				temp = (1 - fabsf(r1)) * fabsf(r2) * occ[i];
				//V[j].x += temp;
				//V[j].y += temp * imPot;
				atomicAdd(&(V[j].x), temp);
				atomicAdd(&(V[j].y), temp * imPot);

				i1 += sign(r1);
				sgCoord(j, i1, i2, m1);
				temp = fabsf(r1) * fabsf(r2) * occ[i];
				//V[j].x += temp;
				//V[j].y += temp * imPot;
				atomicAdd(&(V[j].x), temp);
				atomicAdd(&(V[j].y), temp * imPot);

				i2 -= sign(r2);
				sgCoord(j, i1, i2, m1);
				temp = fabsf(r1) * (1 - fabsf(r2)) * occ[i];
				//V[j].x += temp;
				//V[j].y += temp * imPot;
				atomicAdd(&(V[j].x), temp);
				atomicAdd(&(V[j].y), temp * imPot);

			}
		}
	}
}

__global__ void divideBySinc(cufftComplex* V, int nx, int ny, float_tt PI) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	const int m1 = nx;
	const int m2 = ny;

	if (i < m1 * m2) {
		int i1, i2;
		dbCoord(i1, i2, i, m1);
		iwCoordIp(i1, m1);
		iwCoordIp(i2, m2);

		float_tt y = PI;
		float_tt x = ((float_tt) i1) / ((float_tt) m1) * y;
		x = (x + FLT_EPSILON) / (sinf(x) + FLT_EPSILON);
		y *= ((float_tt) i2) / ((float_tt) m2);
		x *= (y + FLT_EPSILON) / (sinf(y) + FLT_EPSILON);

		V[i].x *= x;
		V[i].y *= x;
	}
}
__global__ void multiplyWithProjectedPotential_d(cufftComplex* V1, cufftComplex* V2, int nx, int ny) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	const int m1 = nx;
	const int m2 = ny;

	if (i < m1 * m2) {
		float_tt V2x = V2[i].x;
		V1[i].x *= V2x;
		V1[i].y *= V2x;
	}
}

__global__ void potential2Transmission(cufftComplex* t, cufftComplex* V, int size) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < size) {
		float_tt Vx = V[i].x;
		float_tt Vy = V[i].y;
		t[i].x = expf(-Vy) * cosf(Vx);
		t[i].y = expf(-Vy) * sinf(Vx);
	}
}

__device__ int sign(float x) {
	int i;
	if (x < 0.f) {
		i = -1;
	} else {
		i = 1;
	}

	return (i);
}

int CUDAFunctions::myGSize(int size) {
	int dev = 0;
	cudaSetDevice(dev);
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, dev);

	const int maxGS = deviceProp.maxGridSize[0] / 2; // HALF of max gridsize allowed by device, it is taken double elsewhere
	const int maxBS = deviceProp.maxThreadsDim[0]; // Maximum blocksize allowed by device.

	int bS = maxBS;
	int gS = size / bS + 1;

	if (gS > maxGS) {
		gS = maxGS;
	}

	if (bS > maxBS) {
		bS = maxBS;
	}

	if ((bS * gS) < size) {
		fprintf( stderr, "    WARNING: Dimensions of the object too large for the GPU.");
	}

	return gS;
}

int CUDAFunctions::myBSize(int size) {
	int dev = 0;
	cudaSetDevice(dev);
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, dev);

	const int maxGS = deviceProp.maxGridSize[0] / 2; // HALF of max gridsize allowed by device, it is taken double elsewhere
	const int maxBS = deviceProp.maxThreadsDim[0]; // Maximum blocksize allowed by device.

	int bS = maxBS;
	int gS = size / bS + 1;

	if (gS > maxGS) {
		gS = maxGS;
	}

	if (bS > maxBS) {
		bS = maxBS;
	}

	if ((bS * gS) < size) {
		fprintf( stderr, "    WARNING: Dimensions of the object too large for the GPU.");
	}

	return bS;
}

}