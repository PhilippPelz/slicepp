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
#include <arrayfire.h>
#include <af/util.h>
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

	_prop(_nx, _ny, c32);
	_prop = af::constant(1, _nx, _ny);
	float_tt scale = _config->Model.dz * PI * GetWavelength();
	//t = exp(-i pi lam k^2 dz)
	af::array s, kx2D, ky2D;
// Tile the arrays to create 2D versions of the k vectors
	kx2D = af::tile(m_kx2, 1, _ny);
	ky2D = af::tile(m_ky2.T(), _nx);
	s = scale*(kx2D + ky2D);
	_prop = af::complex(af::cos(s), af::sin(s));

//	TODO:
	_prop = fftShift(_prop);
	_persist->Save2DDataSet(_prop,"Propagator");
	//_prop = fftShift(_prop);
}

void CBaseWave::ShiftTo(float_tt x, float_tt y){

}
af::array CBaseWave::fftShift(af::array wave){
	return af::shift(wave, wave.dims(0)/2, wave.dims(1)/2);
}

af::array CBaseWave::ifftShift(af::array wave){
	return af::shift(wave, (wave.dims(0)+1)/2, (wave.dims(1)+1)/2);
}

void CBaseWave::ApplyTransferFunction(){
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
void CBaseWave::PropagateToNextSlice()
{
	//_wave_af = fftShift(_wave_af);
	_wave_af = _condition*(_wave_af * _prop);
	//_wave_af = ifftShift(_wave_af);
} /* end propagate */

void CBaseWave::Transmit(af::array t_af)
{
	_wave_af(af::seq(_wave_af.dims(0)*_wave_af.dims(1))) *= t_af(af::seq(_wave_af.dims(0)*_wave_af.dims(1)));
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
	m_kx(_nx);
	m_kx2(_nx);
	m_ky(_ny);
	m_ky2(_ny);
	m_k2(_nx, _ny);

	float_tt ax = _config->Model.dx*_nx;
	float_tt by = _config->Model.dy*_ny;

	m_kx = (af::range(_nx) - _nx/2)/ax;
	m_kx2 = m_kx*m_kx;

	m_ky = (af::range(_ny) - _ny/2)/ax;
	m_ky2 = m_ky*m_ky;

	m_k2max = _nx/(2.0F*ax);
	if (_ny/(2.0F*by) < m_k2max ) m_k2max = _ny/(2.0F*by);
	m_k2max = 2.0/3.0 * m_k2max;
	m_k2max = m_k2max*m_k2max;

	GetK2();
	_k2max = af::constant(m_k2max, _nx, _ny);
	_condition = (m_k2 < _k2max);
	_condition = fftShift(_condition);
}
void  CBaseWave::GetExtents(int& nx, int& ny) const{
	nx = _nx;
	ny = _ny;
}
void CBaseWave::FormProbe(){
	_nx = _config->Model.nx;
	_ny = _config->Model.ny;
	_zero = af::constant(0, _nx, _ny);
	_wave_af = af::complex(_zero, _zero);
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
	BOOST_LOG_TRIVIAL(info)<<format("* Reciprocal space res: dkx=%g, dky=%g (1/A)")%	(1.0/(_nx*_config->Model.dx))%(1.0/(_ny*_config->Model.dy));

	BOOST_LOG_TRIVIAL(info)<<format("* Beams:                %d x %d ")%_nx%_ny;

	BOOST_LOG_TRIVIAL(info)<<format("* Acc. voltage:         %g (lambda=%gA)")%m_v0%(Wavelength(m_v0));

	if (k_fftMeasureFlag == FFTW_MEASURE)
		BOOST_LOG_TRIVIAL(info)<<format("* Probe array:          %d x %d pixels (optimized)")%_nx%_ny;
	else
		BOOST_LOG_TRIVIAL(info)<<format("* Probe array:          %d x %d pixels (estimated)")%_nx%_ny;
	BOOST_LOG_TRIVIAL(info)<<format("*                       %g x %gA")%(_nx*_config->Model.dx)%(_ny*_config->Model.dy);
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

void CBaseWave::GetK2()
{
	af::array kx2D(_nx, _ny), ky2D(_nx, _ny);
	kx2D = af::tile(m_kx2, 1, _ny);
	ky2D = af::tile(m_ky2.T(), _nx);
	m_k2 = kx2D + ky2D;
}

float_tt CBaseWave::GetIntegratedIntensity() const
{
	float_tt intensity;
	intensity = af::sum<float_tt>(af::sqrt(af::real(_wave_af)*af::real(_wave_af) + af::imag(_wave_af)*af::imag(_wave_af)));
	// TODO: divide by px or not?
	return (intensity)/(_nx*_ny);
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
		//_forward.fft(_wave.data());
		_wave_af = af::fft2(_wave_af);

	}
}

// FFT back to realspace, but only if we're currently in Fourier space
void CBaseWave::ToRealSpace()
{
	if (!IsRealSpace())
	{
		m_realSpace = true;
		//		_wave = fftwpp::fft2d::ifftshift(_wave);
		_wave_af = af::ifft2(_wave_af);
		//_backward.fftNormalized(_wave.data());
	}
}

} //end namespace QSTEM
