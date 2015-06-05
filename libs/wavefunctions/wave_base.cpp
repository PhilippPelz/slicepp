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

#include "wave_base.hpp"
using boost::format;

namespace QSTEM
{

void CreateWaveFunctionDataSets(int x, int y, std::vector<int> positions, std::string output_ext)
{
	//TODO user persistence CreateWaveFunctionDataSets
	//	CImageIO imageIO(x, y, "", output_ext);
	//	std::string potDataSetLabel = "Potential";
	//	std::string mulswavDataSetLabel = "mulswav";
	//	imageIO.CreateComplexDataSet(potDataSetLabel, positions);
	//	imageIO.CreateComplexDataSet(mulswavDataSetLabel, positions);
}
CBaseWave::CBaseWave(const ConfigPtr& c,const PersistenceManagerPtr& p) :
				_forward(fftwpp::fft2d(c->Model.nx,c->Model.ny,FFTW_FORWARD)),
				_backward(fftwpp::fft2d(c->Model.nx,c->Model.ny,FFTW_BACKWARD)),
				IWave(c,p)
{
	_nx = c->Model.nx;
	_ny = c->Model.ny;
	m_dx = c->Model.dx;
	m_dy = c->Model.dy;
	m_v0 = c->Beam.EnergykeV;
	m_realSpace = true;
	// TODO: where does this belong?
	//m_electronScale = m_beamCurrent*m_dwellTime*MILLISEC_PICOAMP;
	Initialize(".img", ".img");
}
void CBaseWave::InitializePropagators()
{
	_prop.resize(boost::extents[_nx][_ny]);
	std::fill(_prop.origin(), _prop.origin() + _prop.size(), complex_tt(0, 0));
	float_tt scale = _config->Model.dz * PI * GetWavelength();
	//t = exp(-i pi lam k^2 dz)
	for(int ixa = 0; ixa < _nx; ixa++)
		for(int iya = 0; iya < _ny; iya++){
			float_tt kx = m_kx2[ixa];
			float_tt ky = m_ky2[iya];
			float_tt s = scale*(m_kx2[ixa]+m_ky2[iya]);
			complex_tt tmp = complex_tt(cos(s),sin(s));
			_prop[ixa][iya] = tmp;
			BOOST_LOG_TRIVIAL(trace) << boost::format("p[%d][%d]= %g * exp(i %g) s=%g") % ixa % iya % abs(tmp) % arg(tmp) %s;
		}
	_persist->Save2DDataSet(_prop,"Propagator");
}

void CBaseWave::ShiftTo(float_tt x, float_tt y){

}

void CBaseWave::PropagateToNextSlice()
{
	float_tt wr, wi, tr, ti;
	float_tt scale, t;
	float_tt dzs = 0;

	_wave = fftwpp::fft2d::fftshift(_wave);
#pragma omp parallel for private(wr, wi, tr, ti)
	for(int i = 0; i < _nx; i++)
		for(int j = 0; j < _ny; j++) {
			try {
				if((m_kx2[i] + m_ky2[j]) < m_k2max){
					_wave[i][j] *= _prop[i][j];
				} else {
					_wave[i][j] = complex_tt(0,0);
				}
			} catch(const std::exception& e) {
				std::cerr << e.what();
			}
		} /* end for(ix..) */
	_wave = fftwpp::fft2d::ifftshift(_wave);
} /* end propagate */
void CBaseWave::Transmit(ComplexArray2DView t)
{
	double wr, wi, tr, ti;
	int nx, ny;
	for(int ix = 0; ix < _nx; ix++) {
		for(int iy = 0; iy < _ny; iy++) {
			complex_tt t1 = t[ix][iy];
			wr = _wave[ix][iy].real();
			wi = _wave[ix][iy].imag();
			tr = t1.real();
			ti = t1.imag();
			_wave[ix][iy] *= t1;
			//			BOOST_LOG_TRIVIAL(trace) << boost::format("w=(%g,%g) t=(%2.3f,%2.3f) w*t=(%g,%g)") % wr % wi % tr % ti %
			//					_wave[ix][iy].real() % _wave[ix][iy].imag();
		} /* end for(iy.. ix .) */
	}
} /* end transmit() */
/** Copy constructor - make sure arrays are deep-copied */
CBaseWave::CBaseWave(const CBaseWave &other): CBaseWave(other._config,other._persist)
{
	// TODO: make sure arrays are deep copied
	other.GetSizePixels(_nx, _ny);
	other.GetResolution(m_dx, m_dy);
	m_v0=other.GetVoltage();
	Initialize(".img", ".img");
}


CBaseWave::~CBaseWave()
{
}


void CBaseWave::Initialize(std::string input_ext, std::string output_ext)
{
	m_wavlen = Wavelength(m_v0);
}

void CBaseWave::InitializeKVectors()
{
	m_kx.resize(_nx);
	m_kx2.resize(_nx);
	m_ky.resize(_ny);
	m_ky2.resize(_ny);

	float_tt ax = m_dx*_nx;
	float_tt by = m_dy*_ny;

#pragma omp parallel for
	for(int ixa=0; ixa<_nx; ixa++)
	{
		//#pragma omp critical
		{
			float_tt t = (float_tt)(ixa-_nx/2)/ax;
			//			m_kx[ixa] = (ixa>_nx/2) ? (float_tt)(ixa-_nx)/ax : (float_tt)ixa/ax;
			m_kx[ixa] = t;
			m_kx2[ixa] = t*t;
			BOOST_LOG_TRIVIAL(trace) << boost::format("kx[%d]= %g;  kx2[%d]= %g") % ixa % t % ixa % m_kx2[ixa];
		}
	}
#pragma omp parallel for
	for(int iya=0; iya<_ny; iya++) {
		//#pragma omp critical
		{
			//			m_ky[iya] = (iya>_ny/2) ?	(float_tt)(iya-_ny)/by :(float_tt)iya/by;
			m_ky[iya] = (float_tt)(iya-_ny/2)/ax;
			m_ky2[iya] = m_ky[iya]*m_ky[iya];
		}
	}
	m_k2max = _nx/(2.0F*ax);
	if (_ny/(2.0F*by) < m_k2max ) m_k2max = _ny/(2.0F*by);
	m_k2max = 2.0/3.0 * m_k2max;
	m_k2max = m_k2max*m_k2max;
}
void  CBaseWave::GetExtents(int& nx, int& ny) const{
	nx = _nx;
	ny=_ny;
}
void CBaseWave::FormProbe(){
	_nx = _config->Model.nx;
	_ny = _config->Model.ny;
	_wave.resize(boost::extents[_nx][_ny]);
	_forward = fftwpp::fft2d(_nx,_ny,FFTW_FORWARD);
	_backward = fftwpp::fft2d(_nx,_ny,FFTW_BACKWARD);
	InitializeKVectors();
}
void CBaseWave::DisplayParams()
{
	BOOST_LOG_TRIVIAL(info) <<
			"*****************************  Wave  Parameters **************************************************";
	BOOST_LOG_TRIVIAL(info) <<
			"**************************************************************************************************";
	BOOST_LOG_TRIVIAL(info)<<format("* Real space res.:      %gA (=%gmrad)")%	(1.0/m_k2max)%(GetWavelength()*m_k2max*1000.0);
	BOOST_LOG_TRIVIAL(info)<<format("* Reciprocal space res: dkx=%g, dky=%g (1/A)")%	(1.0/(_nx*m_dx))%(1.0/(_ny*m_dy));

	BOOST_LOG_TRIVIAL(info)<<format("* Beams:                %d x %d ")%_nx%_ny;

	BOOST_LOG_TRIVIAL(info)<<format("* Acc. voltage:         %g (lambda=%gA)")%m_v0%(Wavelength(m_v0));

	if (k_fftMeasureFlag == FFTW_MEASURE)
		BOOST_LOG_TRIVIAL(info)<<format("* Probe array:          %d x %d pixels (optimized)")%_nx%_ny;
	else
		BOOST_LOG_TRIVIAL(info)<<format("* Probe array:          %d x %d pixels (estimated)")%_nx%_ny;
	BOOST_LOG_TRIVIAL(info)<<format("*                       %g x %gA")%(_nx*m_dx)%(_ny*m_dy);
}

/*
//TODO: where does this belong?
inline void CBaseWave::GetElectronScale(float_tt &electronScale)
{
  electronScale=m_electronScale;
}
 */
void CBaseWave::GetSizePixels(int &x, int &y) const
{
	x=_nx;
	y=_ny;
}

void CBaseWave::GetResolution(float_tt &x, float_tt &y) const 
{
	x=m_dx;
	y=m_dy;
}

void CBaseWave::GetPositionOffset(int &x, int &y) const
{
	x=m_detPosX;
	y=m_detPosY;
}

float_tt CBaseWave::GetK2(int ix, int iy) const
{
	return m_kx2[ix]+m_ky2[iy];
}

float_tt CBaseWave::GetIntegratedIntensity() const 
{
	float_tt intIntensity=0;
	//#pragma omp parallel for
	for (int i=0; i<_nx; i++)
		for(int j=0;j<_ny;j++)
		{
			//#pragma omp atomic
			intIntensity+=_wave[i][j].real()*_wave[i][j].real() + _wave[i][j].imag()*_wave[i][j].imag();
		}
	// TODO: divide by px or not?
	return intIntensity/(_nx*_ny);
}

void CBaseWave::ApplyTransferFunction(std::vector<complex_tt> &wave)
{
	// TODO: transfer function should be passed as a 1D vector that is half the size of the wavefunc.
	//       It should be applied by a radial lookup table (with interpolation?)
	//       Alternatively, is it easier to just use a 2D CTF?
	//       Whatever you do, use m_transferFunction as the storage for it.
	int px=GetTotalPixels();

	// multiply wave (in rec. space) with transfer function and write result to imagewave
	ToFourierSpace();
	for (int i=0;i<_nx;i++)
		for (int j=0;j<_ny;j++)
		{
			// here, we apply the CTF:
			// 20140110 - MCS - I think this is where Christoph wanted to apply the CTF - nothing is done ATM.

			// TODO: use these for calculating a radius (to get the CTF value from)
			//ix=i%m_nx;
			//iy=i/m_ny;

			//		wave[i][j] = m_wave[i][j];
		}
	ToRealSpace();
}


/*--------------------- wavelength() -----------------------------------*/
/*
	return the electron wavelength (in Angstroms)
	keep this is one place so I don't have to keep typing in these
	constants (that I can never remember anyhow)

	ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
		(The American Institute of Physics, New York) 1989
		page 4.

	kev = electron energy in keV

 */
float_tt CBaseWave::Wavelength(float_tt kev)
{
	double w;
	const double emass=510.99906; /* electron rest mass in keV */
	const double hc=12.3984244; /* Planck's const x speed of light*/

	/* electron wavelength in Angstroms */
	return hc/sqrt( kev * ( 2*emass + kev ) );
}  /* end wavelength() */

/*
void CBaseWave:WriteBeams(int absolute_slice) {
  static char fileAmpl[32];
  static char filePhase[32];
  static char fileBeam[32];
  static FILE *fp1 = NULL,*fpAmpl = NULL,*fpPhase=NULL;
  int ib;
  static std::vector<int> hbeam,kbeam;
  static float_tt zsum = 0.0f,scale;
  float_tt rPart,iPart,ampl,phase;
  static char systStr[64];
  // static int counter=0;

  if (!muls->lbeams)
    return;  	
}
 */

// FFT to Fourier space, but only if we're current in real space
void CBaseWave::ToFourierSpace()
{
	if (IsRealSpace())
	{
		m_realSpace = false;
		_forward.fft(_wave.data());

	}
}

// FFT back to realspace, but only if we're currently in Fourier space
void CBaseWave::ToRealSpace()
{
	if (!IsRealSpace())
	{
		m_realSpace = true;
		//		_wave = fftwpp::fft2d::ifftshift(_wave);
		_backward.fftNormalized(_wave.data());
	}
}

} //end namespace QSTEM
