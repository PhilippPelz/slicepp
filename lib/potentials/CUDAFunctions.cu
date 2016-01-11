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
#include <boost/format.hpp>

#include <thrust/complex.h>
using boost::format;
namespace slicepp {
CUDAFunctions::CUDAFunctions(){}
CUDAFunctions::CUDAFunctions(superCellBoxPtr info, cModelConfPtr mc) {
	_info = info;
	_mc = mc;

	slicePixels = _mc->n[0] * _mc->n[1];
	_stream = afcu::getStream(af::getDevice());
	_gS = myGSize(slicePixels);
	_gS3D = myGSize(slicePixels* _mc->n[2]);
	_bS = myBSize(slicePixels);
	_bS3D = myBSize(slicePixels* _mc->n[2]);
}

void CUDAFunctions::PotentialToTransmission(cufftComplex* pot, cufftComplex* trans) {
	potential2Transmission<<< _gS, _bS, 0, _stream >>> (pot, trans, slicePixels);
}
void CUDAFunctions::cmul(cufftComplex* a1, cufftComplex* a2) {
	multiplyWithProjectedPotential_d<<< _gS, _bS, 0, _stream >>> (a1,a2,_mc->n[0],_mc->n[1]);
}
void CUDAFunctions::SetComplex2D(cufftComplex* a, float real, float imag) {
	initialValues<<< _gS, _bS , 0, _stream>>> ( a, slicePixels, 0.f, 0.f);
}
void CUDAFunctions::SetComplex3D(cufftComplex* a, float real, float imag) {
	initialValues<<< _gS3D, _bS3D , 0, _stream>>> ( a, slicePixels * _mc->n[2], 0.f, 0.f);
}
void CUDAFunctions::GetAtomicPotential(cufftComplex* V, int Z) {
	createAtomicPotential<<< _gS, _bS, 0, _stream>>> ( V, Z, _mc->n[0], _mc->n[1], _mc->d[0], _mc->d[1],_mc->sigma);
}
void CUDAFunctions::GetSincAtomicPotential(cufftComplex* V, int Z) {
	createAtomicPotential<<< _gS, _bS, 0, _stream >>> ( V, Z, _mc->n[0], _mc->n[1], _mc->d[0], _mc->d[1],_mc->sigma);
	divideBySinc<<< _gS, _bS, 0, _stream >>> ( V, _mc->n[0], _mc->n[1], PI);
}
void CUDAFunctions::GetAtomDeltaFunctions(cufftComplex* V, int Z, int slice, float* xyzPos_d, float* occupancy_d, int* znums_d) {
	int nAtom = _info->znums.size();
	putAtomDeltas<<< myGSize( nAtom ), myBSize( nAtom ), 0, _stream >>> ( V, nAtom, znums_d, Z, xyzPos_d, _mc->ImagPot,
			occupancy_d, slice, _mc->n[0], _mc->n[1], _mc->n[2], _mc->d[0], _mc->d[1], _mc->d[2]);
}
void CUDAFunctions::printPotArray(cufftComplex* V_d, int nx, int ny) {
	cufftComplex *V_host;
	V_host = (cufftComplex *) malloc(nx * ny * sizeof(cufftComplex));
	cuda_assert(cudaMemcpyAsync(V_host, V_d, nx * ny * sizeof(cufftComplex), cudaMemcpyDeviceToHost, afcu::getStream(af::getDevice())));
	cuda_assert(cudaStreamSynchronize(_stream));
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			float x = V_host[i * nx + j].x;
			float y = V_host[i * nx + j].y;
//			if (x > 1e-10 || y > 1e-10)
			cout << "(" << i << "," << j << ") = (" << x << ", " << y << ")" << endl;
		}
	}
	free(V_host);
}

void CUDAFunctions::printFloatArray(float_tt* f, int nx, int ny, int offset) {
	float_tt *f_host;
	f_host = (float_tt *) malloc(nx * ny * sizeof(float_tt));
	cuda_assert(cudaMemcpyAsync((void * )f_host, (void * )f, nx * ny * sizeof(float_tt), cudaMemcpyDeviceToHost, afcu::getStream(af::getDevice())));
	cuda_assert(cudaStreamSynchronize(_stream));
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
	cuda_assert(cudaMemcpyAsync(f_host, p, size * sizeof(int), cudaMemcpyDeviceToHost, afcu::getStream(af::getDevice())));
	cuda_assert(cudaStreamSynchronize(_stream));
	for (int i = 0; i < size; i++) {
		cout << f_host[i] << endl;
	}
	free(f_host);
}
void CUDAFunctions::initArrays() {
//	cuda_assert(cudaStreamSynchronize(_stream));
//	printf("wp 01\n");
//	cuda_assert(cudaMalloc(&xyzPos_d,_info->xyzPos.size()*sizeof(float)));
//	cuda_assert(cudaMalloc(&occupancy_d,_info->occupancy.size()*sizeof(float)));
//	cuda_assert(cudaMalloc(&znums_d,_info->znums.size()*sizeof(int)));
//	auto t1 = af::array(_info->xyzPos.size());
//	auto t2 = af::array(_info->occupancy.size());
//	auto t3 = af::array(_info->znums.size());
//	auto t1 = af::array(_info->xyzPos.size(), _info->xyzPos.data(),afHost);
//	auto t2 = af::array(_info->occupancy.size(), _info->occupancy.data(),afHost);
//	auto t3 = af::array(_info->znums.size(), _info->znums.data(),afHost);
//	t1 *= 1;
//	t2 *= 1;
//	t3 *= 1;
//	af::sync();
//	xyzPos_d = t1.device<float>();
//	occupancy_d = t2.device<float>();
//	znums_d = t3.device<int>();
////	af::sync();
//	printf("wp 02\n");
//	printf("znums_d:  %#08x\n", znums_d);
//	printf("occupancy_d:  %#08x\n", occupancy_d);
//	printf("xyzPos_d:  %#08x\n", xyzPos_d);
////	cuda_assert(cudaDeviceSynchronize());
//	printf("%g %g %g %g %g %g %g %g %g\n",_info->xyzPos[0],_info->xyzPos[1],_info->xyzPos[2],_info->xyzPos[3]
//			,_info->xyzPos[4],_info->xyzPos[5],_info->xyzPos[6],_info->xyzPos[7],_info->xyzPos[8]);
//	cuda_assert(cudaStreamSynchronize(_stream));
//	cuda_assert(cudaMemcpyAsync(xyzPos_d, _info->xyzPos.data(), _info->xyzPos.size() * sizeof(float_tt), cudaMemcpyHostToDevice, _stream));
//	cuda_assert(cudaMemcpyAsync(occupancy_d, _info->occupancy.data(), _info->occupancy.size() * sizeof(float_tt), cudaMemcpyHostToDevice, _stream));
//	cuda_assert(cudaMemcpyAsync(znums_d, _info->znums.data(), _info->znums.size() * sizeof(int), cudaMemcpyHostToDevice, _stream));
//	cuda_assert(cudaStreamSynchronize(_stream));
//	printFloatArray(xyzPos_d,_info->xyzPos.size()/100,100,0);
}
void CUDAFunctions::releaseArrays() {
//	cuda_assert(cudaFree(xyzPos_d));
//	cuda_assert(cudaFree(occupancy_d));
//	cuda_assert(cudaFree(znums_d));
}
void CUDAFunctions::initPotArrays(int slicePixels) {
	_V_elem = af::array(slicePixels, c32);
	_V_atom = af::array(slicePixels, c32);
	_v_accum = af::array(slicePixels, c32);
	_V_elem_ptr = (cufftComplex *) _V_elem.device<af::af_cfloat>();
	_V_atom_ptr = (cufftComplex *) _V_atom.device<af::af_cfloat>();
	_V_accum_ptr = (cufftComplex *) _v_accum.device<af::af_cfloat>();
	initArrays();
}

void CUDAFunctions::unlockArrays() {
	_V_elem.unlock();
	_V_atom.unlock();
	_v_accum.unlock();
}
void CUDAFunctions::limitBandwidth(cufftComplex* f) {
	zeroHighFreq<<< 2 * _gS, _bS, 0, _stream >>>(f,_mc->n[0], _mc->n[1]);
}
__global__ void zeroHighFreq(cufftComplex* f, int dim1, int dim2) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	float mindim = (float) dim1;

	if ((float) dim2 < mindim) {
		mindim = (float) dim2;
	}

	if (i < 2 * dim1 * dim2) {
		int i0 = i / 2;
		int i1, i2;
		dbCoord(i1, i2, i0, dim1);
		iwCoordIp(i1, dim1);
		iwCoordIp(i2, dim2);

		if (((float) (i1 * i1 + i2 * i2) * 9.f / (mindim * mindim)) > 1.f) {
			if ((i % 2) == 0) {
				f[i0].x = 0.f;
			}

			else {
				f[i0].y = 0.f;
			}
		}
	}
}
__global__ void initialValues(cuComplex* V, int size, float_tt initRe, float_tt initIm) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < size) {
		V[i].x = initRe;
		V[i].y = initIm;
	}
}

__global__ void putAtomDeltas(cufftComplex* V, int nAt, int *Z, int Z0, float_tt *xyz, float_tt imPot, float_tt *occ, int s, int nx, int ny,
		int nSlices, float_tt dx, float_tt dy, float_tt dz) //double linear interpolation
		{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
//	printf("Thread:(%d) 	i: %d		nAt:%d\n",threadIdx.x,i,nAt);
	if (i < nAt) {
//		printf("Thread:(%d) 	i: %d		Z[i]:%d\n",threadIdx.x,i,Z[i]);

		if (Z[i] == Z0) {

			const int m1 = nx;
			const int m2 = ny;
			float_tt x1, x2;
//			printf("Thread:(%d) xyz index: %d,%d \n",threadIdx.x, i * 3 + 0,i * 3 + 1);
			x1 = xyz[i * 3 + 0] / dx - 0.5f;
			x2 = xyz[i * 3 + 1] / dy - 0.5f;
			int i3 = (int) (roundf(xyz[i * 3 + 2] / dz));
//			printf("x1,m1:(%3.3f,%3.3f)		x2,m2:(%3.3f,%3.3f)	x3,slice:(%d,%d)\n",x1,m1,x2,m2,i3,s);
//			printf("%3.3f = (xyz[%04d]=%3.3f) / %3.3f - 0.5f; %3.3f = (xyz[%04d]=%3.3f) / %3.3f - 0.5f; %3.3f = (int) (roundf(xyz[%04d]=%3.3f / %3.3f - 0.5f));\n",
//			x1,i*3,xyz[i * 3 + 0],dx,
//			x2,i*3+1,xyz[i * 3 + 1],dy,
//			i3,i*3+2,xyz[i * 3 + 2],dz);

			if (((x1 > 1.f) && (x1 < ((float_tt) (m1 - 2)))) && ((x2 > 1.f) && (x2 < ((float_tt) (m2 - 2)))) && ((i3 > s - 1) && (i3 <= s))) {
//				if(Z0==79)
//				printf("s=%d	Z=%d (xyz)=(%3.3g,%3.3g,%3.3g)	(m1,m2,i3)=(%d,%d,%d)	(dx,dy,dz)=(%3.3g,%3.3g,%3.3g)\n",
//						s,Z0,xyz[i * 3 + 0],xyz[i * 3 + 1],xyz[i * 3 + 2],
//						m1,m2,i3,
//						dx,dy,dz);
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
//				if(j > m1*m2)
//					printf("Thread:(%d) V index: %d \n",threadIdx.x, j);
				atomicAdd(&(V[j].x), temp);
				atomicAdd(&(V[j].y), temp * imPot);

				i2 += sign(r2);
				sgCoord(j, i1, i2, m1);
				temp = (1 - fabsf(r1)) * fabsf(r2) * occ[i];
				//V[j].x += temp;
				//V[j].y += temp * imPot;
//				if(j > m1*m2)
//					printf("Thread:(%d) V index: %d \n",threadIdx.x, j);
				atomicAdd(&(V[j].x), temp);
				atomicAdd(&(V[j].y), temp * imPot);

				i1 += sign(r1);
				sgCoord(j, i1, i2, m1);
				temp = fabsf(r1) * fabsf(r2) * occ[i];
				//V[j].x += temp;
				//V[j].y += temp * imPot;
//				if(j > m1*m2)
//					printf("Thread:(%d) V index: %d \n",threadIdx.x, j);
				atomicAdd(&(V[j].x), temp);
				atomicAdd(&(V[j].y), temp * imPot);

				i2 -= sign(r2);
				sgCoord(j, i1, i2, m1);
				temp = fabsf(r1) * (1 - fabsf(r2)) * occ[i];
				//V[j].x += temp;
				//V[j].y += temp * imPot;
//				if(j > m1*m2)
//					printf("Thread:(%d) V index: %d \n",threadIdx.x, j);
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
		x /= (nx*ny);
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
struct potential2TransmissionFunctor {
	__host__ __device__
	thrust::complex<float> operator()(const thrust::complex<float>& x, const thrust::complex<float>& y) const {
		return thrust::complex<float>(expf(-x.imag()) * cosf(x.real()), expf(-x.imag()) * sinf(x.real()));
	}
};
__global__ void potential2Transmission(cufftComplex* t, cufftComplex* V, int size) {
	const int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < size) {
		float_tt Vx = V[i].x;
		float_tt Vy = V[i].y;
		t[i].x = expf(-Vy) * cosf(Vx);
		t[i].y = expf(-Vy) * sinf(Vx);
//		if(i%100)
//			printf("i=%-5d V=(%-10f,%-10f) t=(%-10f,%-10f)\n",i,Vx,Vy,t[i].x,t[i].y);
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
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, af::getDevice());

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
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, af::getDevice());

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
