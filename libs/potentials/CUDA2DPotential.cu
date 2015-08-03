/*
 * CUDA2DPotential.cpp
 *
 *  Created on: Jul 29, 2015
 *      Author: wenxuan
 */

#include "CUDA2DPotential.hpp"

namespace QSTEM {

CUDA2DPotential::CUDA2DPotential(const ConfigPtr& c,const PersistenceManagerPtr& persist): CPotential(c, persist) {
	// TODO Auto-generated constructor stub

}

CUDA2DPotential::~CUDA2DPotential() {
	// TODO Auto-generated destructor stub
}

void CUDA2DPotential::MakeSlices(superCellBoxPtr info){
	time_t time0, time1;
	cufftComplex *potential;
	int slicePixels, numAtoms, numAtUnique;
	float_tt imPot = _c->Wave.imPot;
	cufftHandle cufftPlanBatch;
	struct cublasContext *cublasHandle;
	cufftPlan2d (&cufftPlanBatch, _c->Model.nx, _c->Model.ny, CUFFT_C2C);
	slicePixels = _c->Model.nx * _c->Model.ny;
	numAtoms = info->atoms.size();
	numAtUnique = info->uniqueatoms.size();
	cudaMalloc ((void**) &potential, _c->Model.nSlices * slicePixels * sizeof(cufftComplex)) ;
	time(&time0);
	for (int islice = 0; islice < _c->Model.nSlices; islice++){
		phaseGrating(cublasHandle, potential[islice * slicePixels], numAtoms, numAtUnique, info->xyzPos, imPot, info->znums, info->uniqueatoms, info->occupancy, cufftPlanBatch, islice);
	}

	_t_af.write(potential, _c->Model.nSlices * slicePixels * sizeof(cufftComplex), afDevice);
	time(&time1);
	BOOST_LOG_TRIVIAL(info)<< format( "%g sec used for real space potential calculation (%g sec per atom)")
	% difftime(time1, time0)%( difftime(time1, time0) / info->atoms.size());

	if (_c->Output.SavePotential)
		_persist->SavePotential(_t_af);
	if (_c->Output.SaveProjectedPotential){
		if(!_persist->potSaved){
			_t_af.host(_t.data());
		}
		WriteProjectedPotential();
	}
}

void  CUDA2DPotential::phaseGrating(struct cublasContext *cublasHandle, cufftComplex* V_d, int nAt, int nZ, std::vector<float_tt> xyz_d, float_tt imPot, std::vector<int> Z_d, std::vector<int> Zlist, std::vector<float_tt> occ_d, cufftHandle cufftPlanBatch, int s ){
	cufftComplex *V1_d, *V2_d, alpha;
	int slicePixels = _c->Model.nx * _c->Model.ny;
	int gS = myGSize(slicePixels);
	int bS = myBSize(slicePixels * _c->Model.nSlices );
	int gS2D = myGSize(slicePixels);
	alpha.x = 1.f;
	alpha.y = 0.f;
	cudaMalloc((void**) &V1_d, slicePixels * sizeof (cufftComplex )); // V1_d m123
	cudaMalloc((void**) &V2_d, slicePixels * sizeof (cufftComplex ));
	initialValues <<< gS * 2, bS >>> ( V_d,  slicePixels, 0.f, 0.f );
	for ( int j = 0; j < nZ; j++ )
	{
		initialValues <<< gS * 2, bS >>> ( V1_d, slicePixels, 0.f, 0.f );
		squareAtoms_d <<< myGSize( nAt ), myBSize( nAt ) >>> ( V1_d, nAt, &Z_d[0], Zlist[j], &xyz_d[0], imPot, &occ_d[0], s, _c->Model.nx, _c->Model.ny, _c->Model.nSlices, _c->Model.dx, _c->Model.dy, _c->Model.dz);
		projectedPotential_d <<< gS2D, bS >>> ( V2_d, Zlist[j], _c->Model.nx, _c->Model.ny, _c->Model.dx, _c->Model.dy);
		divideBySinc <<< gS2D, bS >>> ( V2_d, _c->Model.nx, _c->Model.ny, PI);
		cufftExecC2C( cufftPlanBatch, V1_d, V1_d, CUFFT_FORWARD);

		multiplyWithProjectedPotential_d <<< gS, bS >>> ( V1_d, V2_d, _c->Model.nx, _c->Model.ny);
		cufftExecC2C( cufftPlanBatch, V1_d, V1_d, CUFFT_INVERSE );
		cublasCaxpy(cublasHandle, slicePixels, &alpha, V1_d, 1, V_d, 1);
	}
	cudaFree (V1_d);
	cudaFree (V2_d);

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
			x1 = xyz[i * 3 + 0] / dx + ( (float_tt) m1 ) * 0.5f - 0.5f;
			x2 = xyz[i * 3 + 1] / dy + ( (float_tt) m2 ) * 0.5f - 0.5f;
			int i3 = (int) ( roundf( xyz[i * 3 + 2] /dz + ( (float_tt) m3 ) * 0.5f - 0.5f ) );

			if ( ( ( x1 > 1.f ) && ( x1 < ( (float) ( m1 - 2 ) ) ) )  &&  ( ( x2 > 1.f ) && ( x2 < ( (float) ( m2 - 2 ) ) ) )  &&  ( ( i3 > s-1 ) && ( i3 <= s ) ) )
			{
				int i1 = (int) roundf( x1 );
				int i2 = (int) roundf( x2 );
				int j;
				float_tt r1 = x1 - ( (float_tt) i1 );
				float_tt r2 = x2 - ( (float_tt) i2 );
				float_tt temp;

				sgCoord(j, i1, i2, m1);
				temp = ( 1 - fabsf( r1 ) ) * ( 1 - fabsf( r2 ) ) * occ[i];
				atomicAdd( &( V[j].x ), temp );
				atomicAdd( &( V[j].y ), temp * imPot );

				i2 += mySignum_d( r2 );
				sgCoord(j, i1, i2, m1);
				temp = ( 1 - fabsf( r1 ) ) * fabsf( r2 ) * occ[i];
				atomicAdd( &( V[j].x ), temp );
				atomicAdd( &( V[j].y ), temp * imPot );

				i1 += mySignum_d( r1 );
				sgCoord(j, i1, i2, m1);
				temp = fabsf( r1 ) * fabsf( r2 ) * occ[i];
				atomicAdd( &( V[j].x ), temp );
				atomicAdd( &( V[j].y ), temp * imPot );

				i2 -= mySignum_d( r2 );
				sgCoord(j, i1, i2, m1);
				temp = fabsf( r1 ) * ( 1 - fabsf( r2 ) ) * occ[i];
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
		float V2x = V2[i].x;
		V1[i].x *= V2x;
		V1[i].y *= V2x;
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



void CUDA2DPotential::myGBSize( int* gbS, int size )
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

	gbS[0] = gS;
	gbS[1] = bS;
}

int CUDA2DPotential::myGSize( int size )
{
	int* gbS;
	gbS = ( int* ) malloc ( 2 * sizeof ( int ) );
	myGBSize( gbS, size );
	return( gbS[0] );
}

int CUDA2DPotential::myBSize( int size )
{
	int* gbS;
	gbS = ( int* ) malloc ( 2 * sizeof ( int ) );
	myGBSize( gbS, size );
	return( gbS[1] );
}


} /* namespace QSTEM */


