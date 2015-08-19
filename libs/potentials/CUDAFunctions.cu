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

using namespace std;
namespace QSTEM {

CUDAFunctions::CUDAFunctions() {
	// TODO Auto-generated constructor stub

}

void  CUDAFunctions::phaseGrating(cufftComplex* V_d, int nAt, int nZ, float_tt*  xyz_d, float_tt imPot, int* Z_d, int* Zlist, float_tt* occ_d,  int s, int nx, int ny, int nSlices, float_tt dx, float_tt dy, float_tt dz ){
	cufftComplex alpha;
	int slicePixels = nx * ny;
	//printFloatArray(xyz_d, nAt, 3);

	af::sync();
	const int gS = myGSize(slicePixels);
	const int bS = myBSize(slicePixels * nSlices );
	const int gS2D = myGSize(slicePixels);
	alpha.x = 1.f;
	alpha.y = 0.f;


	initialValues <<< gS * 2, bS >>> ( V_d,  slicePixels, 0.f, 0.f);
	initialValues <<< gS * 2, bS >>> (V3_d,  slicePixels, 0.f, 0.f);
	for ( int j = 0; j < nZ; j++ )
	{
		initialValues <<< gS * 2, bS >>> ( V1_d, slicePixels, 0.f, 0.f );
		initialValues <<< gS * 2, bS >>> ( V2_d, slicePixels, 0.f, 0.f );
		squareAtoms_d <<< myGSize( nAt ), myBSize( nAt ) >>> ( V1_d, nAt, Z_d, Zlist[j], xyz_d, imPot, occ_d, s, nx, ny, nSlices, dx, dy, dz);
		projectedPotential_d <<< gS2D, bS >>> ( V2_d, Zlist[j], nx, ny, dx, dy);
		divideBySinc <<< gS2D, bS >>> ( V2_d, nx, ny, PI);
		cuda_assert(cudaDeviceSynchronize());

		V1.unlock();
		V2.unlock();
		V3.unlock();
		V1 = af::fft(V1);
		V1 *= V2;
		V1 = af::ifft(V1);
		V1 *= sqrtf(slicePixels);
		V3 = V1 + V3;
		af::sync();
		V1_d = (cufftComplex *)V1.device<af::af_cfloat>();
		V2_d = (cufftComplex *)V2.device<af::af_cfloat>();
		V3_d = (cufftComplex *)V3.device<af::af_cfloat>();

		potential2Transmission <<< gS, bS>>> (V_d, V3_d, slicePixels);
	}
}

void CUDAFunctions::printPotArray(cufftComplex* V_d, int nx, int ny){
	cufftComplex *V_host;
	V_host = (cufftComplex *)malloc(nx * ny * sizeof(cufftComplex));
	cuda_assert(cudaMemcpy(V_host, V_d, nx * ny * sizeof(cufftComplex), cudaMemcpyDeviceToHost));
	for (int i = 0; i < nx; i++){
		for (int j =0; j < ny; j++){
			float x = V_host[i*nx + j].x;
			float y = V_host[i*nx + j].y;
			if(x>1e-10 || y>1e-10)
				cout<<"("<<x<<", "<<y<<")"<<endl;
		}
	}
	free(V_host);
}

void CUDAFunctions::printFloatArray(float_tt* f, int nx, int ny, int offset){
	char *f_host;
	f_host = (char *)malloc(nx * ny * sizeof(float_tt));
	cuda_assert(cudaMemcpy((void *)f_host, (void *)f, nx * ny * sizeof(float_tt), cudaMemcpyDeviceToHost));
	f_host += offset;
	float_tt *ff = (float_tt *)f_host;
	for (int i = 0; i < nx; i++){
		for (int j =0; j < ny; j++){
			cout<<ff[i*ny + j]<<" ";
		}
		cout<<endl;
	}
	free(f_host);
}

void CUDAFunctions::printIntArray(int* p, int size){
	int *f_host;
	f_host = (int *)malloc(size * sizeof(int));
	cuda_assert(cudaMemcpy(f_host, p, size * sizeof(int), cudaMemcpyDeviceToHost));
	for (int i = 0; i < size; i++){
		cout<<f_host[i]<<endl;
	}
	free(f_host);
}
void CUDAFunctions::initPotArrays(int slicePixels){
	V1 = af::array(slicePixels, c32);
	V2 = af::array(slicePixels, c32);
	V3 = af::array(slicePixels, c32);
	af::sync();
	V1_d = (cufftComplex *)V1.device<af::af_cfloat>();
	V2_d = (cufftComplex *)V2.device<af::af_cfloat>();
	V3_d = (cufftComplex *)V3.device<af::af_cfloat>();
}
void CUDAFunctions::unlockArrays(){
	V1.unlock();
	V2.unlock();
	V3.unlock();
}
__global__ void initialValues ( cuComplex* V, int size, float_tt initRe, float_tt initIm )
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i < 2 * size )
    {
        if ( ( i % 2 ) == 0 )
        { V[i / 2].x = initRe; }

        else
        { V[i / 2].y = initIm; }
    }
}

__global__ void squareAtoms_d ( cufftComplex* V, int nAt, int *Z, int Z0, float_tt *xyz, float_tt imPot, float_tt *occ, int s, int nx, int ny, int nSlices, float_tt dx, float_tt dy, float_tt dz)//double linear interpolation
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if ( i < nAt )
	{
		if ( Z[i] == Z0 )
		{
			const int m1 = nx;
			const int m2 = ny;
			const int m3 = nSlices;
			float_tt x1, x2;
			x1 = xyz[i * 3 + 0] / dx - 0.5f;
			x2 = xyz[i * 3 + 1] / dy - 0.5f;
			int i3 = (int) ( roundf( xyz[i * 3 + 2] /dz - 0.5f ) );

			if ( ( ( x1 > 1.f ) && ( x1 < ( (float_tt) ( m1 - 2 ) ) ) )  &&  ( ( x2 > 1.f ) && ( x2 < ( (float_tt) ( m2 - 2 ) ) ) )  &&  ( ( i3 > s-1 ) && ( i3 <= s ) ) )
			{
				int i1 = (int) roundf( x1 );
				int i2 = (int) roundf( x2 );
				int j;
				float_tt r1 = x1 - ( (float_tt) i1 );
				float_tt r2 = x2 - ( (float_tt) i2 );
				float_tt temp;

				sgCoord(j, i1, i2, m1);
				temp = ( 1 - fabsf( r1 ) ) * ( 1 - fabsf( r2 ) ) * occ[i];
				//V[j].x += temp;
				//V[j].y += temp * imPot;
				atomicAdd( &( V[j].x ), temp );
				atomicAdd( &( V[j].y ), temp * imPot );

				i2 += mySignum_d( r2 );
				sgCoord(j, i1, i2, m1);
				temp = ( 1 - fabsf( r1 ) ) * fabsf( r2 ) * occ[i];
				//V[j].x += temp;
				//V[j].y += temp * imPot;
				atomicAdd( &( V[j].x ), temp );
				atomicAdd( &( V[j].y ), temp * imPot );

				i1 += mySignum_d( r1 );
				sgCoord(j, i1, i2, m1);
				temp = fabsf( r1 ) * fabsf( r2 ) * occ[i];
				//V[j].x += temp;
				//V[j].y += temp * imPot;
				atomicAdd( &( V[j].x ), temp );
				atomicAdd( &( V[j].y ), temp * imPot );

				i2 -= mySignum_d( r2 );
				sgCoord(j, i1, i2, m1);
				temp = fabsf( r1 ) * ( 1 - fabsf( r2 ) ) * occ[i];
				//V[j].x += temp;
				//V[j].y += temp * imPot;
				atomicAdd( &( V[j].x ), temp );
				atomicAdd( &( V[j].y ), temp * imPot );

			}
		}
	}
}

__global__ void divideBySinc ( cufftComplex* V, int nx, int ny, float_tt PI)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	const int m1 = nx;
	const int m2 = ny;

    if ( i <  m1 * m2 )
	{
		int i1, i2;
		dbCoord ( i1, i2, i, m1 );
        iwCoordIp ( i1, m1 );
        iwCoordIp ( i2, m2 );

		float_tt y = PI;
		float_tt x = ( (float_tt) i1 ) / ( (float_tt) m1 ) * y;
		x  = ( x + FLT_EPSILON ) / ( sinf( x ) + FLT_EPSILON );
		y *= ( (float_tt) i2 ) / ( (float_tt) m2 );
		x *= ( y + FLT_EPSILON ) / ( sinf( y ) + FLT_EPSILON );

		V[i].x *= x;
		V[i].y *= x;
	}
}
__global__ void multiplyWithProjectedPotential_d ( cufftComplex* V1, cufftComplex* V2, int nx, int ny)
{
	const int i = blockIdx.x * blockDim.x + threadIdx.x;
	const int m1 = nx;
	const int m2 = ny;

	if ( i < m1 * m2  )
	{
		float_tt V2x = V2[i].x;
		V1[i].x *= V2x;
		V1[i].y *= V2x;
	}
}

__global__ void potential2Transmission ( cufftComplex* t, cufftComplex* V, int size )
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i < size )
    {
        float_tt Vx = V[i].x;
        float_tt Vy = V[i].y;
        t[i].x = expf ( -Vy ) * cosf ( Vx );
        t[i].y = expf ( -Vy ) * sinf ( Vx );
    }
}

__device__ int mySignum_d( float x )
{
	int i;
	if ( x < 0.f )
	{    i = -1; }
	else
	{    i =  1; }

	return( i );
}


int CUDAFunctions::myGSize( int size )
{
    int dev = 0;
    cudaSetDevice(dev);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);

    const int maxGS = deviceProp.maxGridSize[0]/2 ; // HALF of max gridsize allowed by device, it is taken double elsewhere
    const int maxBS = deviceProp.maxThreadsDim[0]; // Maximum blocksize allowed by device.

    int bS = maxBS;
    int gS = size / bS + 1;

    if ( gS > maxGS )
    {    gS = maxGS; }

    if ( bS > maxBS )
    {    bS = maxBS; }

	if ( ( bS * gS ) < size )
    {    fprintf ( stderr, "    WARNING: Dimensions of the object too large for the GPU." ); }

	return gS;
}

int CUDAFunctions::myBSize( int size )
{
    int dev = 0;
    cudaSetDevice(dev);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);

    const int maxGS = deviceProp.maxGridSize[0]/2 ; // HALF of max gridsize allowed by device, it is taken double elsewhere
    const int maxBS = deviceProp.maxThreadsDim[0]; // Maximum blocksize allowed by device.

    int bS = maxBS;
    int gS = size / bS + 1;

    if ( gS > maxGS )
    {    gS = maxGS; }

    if ( bS > maxBS )
    {    bS = maxBS; }

	if ( ( bS * gS ) < size )
    {    fprintf ( stderr, "    WARNING: Dimensions of the object too large for the GPU." ); }

	return bS;
}

}
