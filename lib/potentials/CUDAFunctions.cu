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

CUDAFunctions::CUDAFunctions(superCellBoxPtr info, cModelConfPtr mc) {
	_info = info;
	_mc = mc;

	slicePixels = _mc->n[0] * _mc->n[1];
	gS2D = myGSize(slicePixels);
}

//void CUDAFunctions::GetPhaseGrating(cufftComplex* V_slice, int slice, std::map<int, cufftComplex*> & atomPot) {
//	int nAtom = _info->atoms.size();
//
//	SetComplex2D( V_slice, 0.f, 0.f);
//	SetComplex2D( _V_accum_ptr, 0.f, 0.f);
//	for (int& Z : _info->uniqueZ) {
//		SetComplex2D( _V_elem_ptr, 0.f, 0.f);
//		SetComplex2D( _V_atom_ptr, 0.f, 0.f);
//
//		GetAtomDeltaFunctions(_V_elem_ptr, Z, slice);
//
//		cuda_assert(cudaDeviceSynchronize());
//
//		_V_elem.unlock();
//		_v_accum.unlock();
//
//		_V_elem = af::fft(_V_elem);
//		_V_elem *= atomPot[Z];
//		_V_elem = af::ifft(_V_elem);
//
//		_V_elem *= slicePixels;
//		_v_accum = _V_elem + _v_accum;
//		af::sync();
//
//		PotentialToTransmission(V_slice, _V_accum_ptr);
//	}
//}
void CUDAFunctions::PotentialToTransmission(cufftComplex* pot, cufftComplex* trans){
    int af_id = af::getDevice();
    cudaStream_t af_stream = afcu::getStream(af_id);
	int slicePixels = _mc->n[0] * _mc->n[1];
	const int gS = myGSize(slicePixels);
	const int bS = myBSize(slicePixels);
	potential2Transmission<<< gS, bS, 0, af_stream >>> (pot, trans, slicePixels);
}
void CUDAFunctions::cmul(cufftComplex* a1, cufftComplex* a2){
    int af_id = af::getDevice();
    cudaStream_t af_stream = afcu::getStream(af_id);
	int slicePixels = _mc->n[0] * _mc->n[1];
	const int gS = myGSize(slicePixels);
	const int bS = myBSize(slicePixels);
	multiplyWithProjectedPotential_d<<< gS, bS, 0, af_stream >>> (a1,a2,_mc->n[0],_mc->n[1]);
}
void CUDAFunctions::SetComplex2D(cufftComplex* a, float real, float imag){
	int slicePixels = _mc->n[0] * _mc->n[1];
	const int gS = myGSize(slicePixels);
	const int bS = myBSize(slicePixels);
    int af_id = af::getDevice();
    cudaStream_t af_stream = afcu::getStream(af_id);
//    printf("pixels: %d	gs: %d,		bs: %d a: %#08x\n",slicePixels,gS,bS,a);
	initialValues<<< gS, bS , 0,  af_stream>>> ( a, slicePixels, 0.f, 0.f);
}
void CUDAFunctions::SetComplex3D(cufftComplex* a, float real, float imag){
	int slicePixels = _mc->n[0] * _mc->n[1];
	const int gS = myGSize(slicePixels * _mc->n[2]);
	const int bS = myBSize(slicePixels * _mc->n[2]);
    int af_id = af::getDevice();
    cudaStream_t af_stream = afcu::getStream(af_id);

	initialValues<<< gS, bS , 0,  af_stream>>> ( a, slicePixels * _mc->n[2], 0.f, 0.f);
}
void CUDAFunctions::GetAtomicPotential(cufftComplex* V, int Z) {
	int slicePixels = _mc->n[0] * _mc->n[1];
	const int bS = myBSize(slicePixels);
	const int gS2D = myGSize(slicePixels);
    int af_id = af::getDevice();
//    BOOST_LOG_TRIVIAL(info)<< format("device id: %d") % af_id;
    cudaStream_t af_stream = afcu::getStream(af_id);
	createAtomicPotential<<< gS2D, bS, 0,  af_stream>>> ( V, Z, _mc->n[0], _mc->n[1], _mc->d[0], _mc->d[1],_mc->sigma);
}
void CUDAFunctions::GetSincAtomicPotential(cufftComplex* V, int Z) {
	int slicePixels = _mc->n[0] * _mc->n[1];
	const int bS = myBSize(slicePixels);
	const int gS2D = myGSize(slicePixels);
    int af_id = af::getDevice();
    cudaStream_t af_stream = afcu::getStream(af_id);
//    printf("sigma: %g",_mc->sigma);
	createAtomicPotential<<< gS2D, bS, 0,  af_stream  >>> ( V, Z, _mc->n[0], _mc->n[1], _mc->d[0], _mc->d[1],_mc->sigma);
	divideBySinc<<< gS2D, bS, 0,  af_stream  >>> ( V, _mc->n[0], _mc->n[1], PI);
}
void CUDAFunctions::GetAtomDeltaFunctions(cufftComplex* V, int Z, int slice) {
	int nAtom = _info->znums.size();
    int af_id = af::getDevice();
    cudaStream_t af_stream = afcu::getStream(af_id);
	putAtomDeltas<<< myGSize( nAtom ), myBSize( nAtom ), 0,  af_stream  >>> ( V, nAtom, znums_d, Z, xyzPos_d, _mc->ImagPot,
			occupancy_d, slice, _mc->n[0], _mc->n[1], _mc->n[2], _mc->d[0], _mc->d[1], _mc->d[2]);
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
void CUDAFunctions::initArrays(){

	auto t1 = af::array(_info->xyzPos.size(),f32);
	auto t2 = af::array(_info->occupancy.size(),f32);
	auto t3 = af::array(_info->znums.size(),s32);

	xyzPos_d = t1.device<float>();
	occupancy_d = t2.device<float>();
	znums_d = t3.device<int>();

//	cuda_assert(cudaMalloc((void**)&xyzPos_d, _info->xyzPos.size()*sizeof(float_tt)));
//	cuda_assert(cudaMalloc((void**)&occupancy_d, _info->occupancy.size()*sizeof(float_tt)));
//	cuda_assert(cudaMalloc((void**)&znums_d, _info->znums.size()*sizeof(int)));

	printf("znums_d:  %#08x\n", znums_d);
	printf("occupancy_d:  %#08x\n", occupancy_d);
	printf("xyzPos_d:  %#08x\n", xyzPos_d);

	cuda_assert(cudaMemcpy(xyzPos_d, _info->xyzPos.data(), _info->xyzPos.size()*sizeof(float_tt), cudaMemcpyHostToDevice));
	cuda_assert(cudaMemcpy(occupancy_d, _info->occupancy.data(),  _info->occupancy.size()*sizeof(float_tt), cudaMemcpyHostToDevice));
	cuda_assert(cudaMemcpy(znums_d, _info->znums.data(), _info->znums.size()*sizeof(int), cudaMemcpyHostToDevice));

//	xyzPos = af::array(_info->xyzPos.size(),_info->xyzPos.data());
//	occupancy = af::array(_info->occupancy.size(),_info->occupancy.data());
//	znums = af::array(_info->znums.size(),_info->znums.data());
//	xyzPos *= 1;
//	occupancy *= 1;
//	znums *= 1;
//	BOOST_LOG_TRIVIAL(info)<< format("sizes: xyz %d occ %d znums %d atoms %d") % _info->xyzPos.size() %
//			_info->occupancy.size() % _info->znums.size() % _info->atoms.size();

//	xyzPos_d = xyzPos.device<float_tt>();
//	occupancy_d = occupancy.device<float_tt>();
//	znums_d = znums.device<int>();
//	af::sync();
//	BOOST_LOG_TRIVIAL(info)<< format("xyzPos_d: %d			size: %d") % xyzPos_d % xyzPos.dims(0);
//	BOOST_LOG_TRIVIAL(info)<< format("occupancy_d: %d		size: %d") % occupancy_d% occupancy.dims(0);
//	BOOST_LOG_TRIVIAL(info)<< format("znums_d: %d			size: %d") % znums_d% znums.dims(0);
}
void CUDAFunctions::releaseArrays(){
	cuda_assert ( cudaFree ( xyzPos_d ) );
	cuda_assert ( cudaFree ( occupancy_d ) );
	cuda_assert ( cudaFree ( znums_d ) );
}
void CUDAFunctions::initPotArrays(int slicePixels) {
	_V_elem = af::array(slicePixels, c32);
	_V_atom = af::array(slicePixels, c32);
	_v_accum = af::array(slicePixels, c32);
	af::sync();
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

__global__ void initialValues(cuComplex* V, int size, float_tt initRe, float_tt initIm) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i <  size) {
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
struct potential2TransmissionFunctor
{
  __host__ __device__
  thrust::complex<float> operator()(const thrust::complex<float>& x, const thrust::complex<float>& y) const
  {
    return thrust::complex<float>(expf(-x.imag())*cosf(x.real()),expf(-x.imag())*sinf(x.real()));
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
