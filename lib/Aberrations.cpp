/*
 * Aberrations.cpp
 *
 *  Created on: Jan 13, 2016
 *      Author: philipp
 */

#include "Aberrations.hpp"

#include "boost/format.hpp"
#include <boost/log/trivial.hpp>
using boost::format;
namespace slicepp {

Aberrations::Aberrations(float_tt lambda,AberrationConfig c,string name) {
	_lambdaAngstrom = lambda;
	_c = c;
	_name = name;
}
void Aberrations::DisplayParams() {
	BOOST_LOG_TRIVIAL(info)<< _name << " aberrations:";
	BOOST_LOG_TRIVIAL(info)<< format("* C_3 (C_s):            %g mm") % (_c.Cs*1e-7);
	BOOST_LOG_TRIVIAL(info)<< format("* C_1 (Defocus):        %g nm") % (0.1*_c.Defocus ) ;
	BOOST_LOG_TRIVIAL(info)<< format("* Astigmatism:          %g nm, %g deg") % (0.1*_c.Astigmatism) % (RAD2DEG*_c.AstigmatismAngle);

	if (_c.a33 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_3,3:                %g nm, phi=%g deg") %(_c.a33*1e-1)%(_c.phi33*RAD2DEG);
	if (_c.a31 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_3,1:                %g nm, phi=%g deg") %(_c.a31*1e-1)%(_c.phi31*RAD2DEG);

	if (_c.a44 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_4,4:                %g um, phi=%g deg") %(_c.a44*1e-4)%(_c.phi44*RAD2DEG);
	if (_c.a42 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_4,2:                %g um, phi=%g deg") %(_c.a42*1e-4)%(_c.phi42*RAD2DEG);

	if (_c.a55 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_5,5:                %g um, phi=%g deg") %(_c.a55*1e-4)%(_c.phi55*RAD2DEG);
	if (_c.a53 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_5,3:                %g um, phi=%g deg") %(_c.a53*1e-4)%(_c.phi53*RAD2DEG);
	if (_c.a51 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_5,1:                %g um, phi=%g deg") %(_c.a51*1e-4)%(_c.phi51*RAD2DEG);

	if (_c.a66 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_6,6:                %g um, phi=%g deg") %(_c.a66*1e-7)%(_c.phi66*RAD2DEG);
	if (_c.a64 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_6,4:                %g um, phi=%g deg") %(_c.a64*1e-7)%(_c.phi64*RAD2DEG);
	if (_c.a62 > 0)
		BOOST_LOG_TRIVIAL(info)<< format("* a_6,2:                %g um, phi=%g deg") %(_c.a62*1e-7)%(_c.phi62*RAD2DEG);
	if (_c.C5 != 0)
		BOOST_LOG_TRIVIAL(info)<< format("* C_5:                  %g mm") % (_c.C5*1e-7);

	BOOST_LOG_TRIVIAL(info)<< format("* C_c:                  %g mm") %(_c.Cc*1e-7);
}
Aberrations::~Aberrations() {
	// TODO Auto-generated destructor stub
}

af::array Aberrations::getPhasePlate(af::array& k,af::array& kx,af::array& ky, PersistenceManagerPtr p) {
//	printf("lambda: %f\n",_lambdaAngstrom);
	auto phi = af::atan2(kx * _lambdaAngstrom, ky * _lambdaAngstrom);
	auto k0 = k*_lambdaAngstrom;
//	p->Save2DDataSet(phi,"phi");
	auto chi = af::constant(0, k.dims(), f32);
	chi += 0.5 * af::pow(k0, 2) * (_c.Defocus + _c.Cc * _c.dE_E + (_c.Astigmatism * af::cos(2.0 * (phi - _c.AstigmatismAngle))));
	if ((_c.a33 > 0) || (_c.a31 > 0))
		chi += af::pow(k0, 3) * (_c.a33 * af::cos(3.0 * (phi - _c.phi33)) + _c.a31 * af::cos(phi - _c.phi31)) / 3.0;
	if ((_c.a44 > 0) || (_c.a42 > 0) || (_c.Cs != 0))
		chi += af::pow(k0, 4) * (_c.a44 * af::cos(4.0 * (phi - _c.phi44)) + _c.a42 * af::cos(2.0 * (phi - _c.phi42)) + _c.Cs) / 4.0;
	if ((_c.a55 > 0) || (_c.a53 > 0) || (_c.a51 > 0))
		chi += af::pow(k0, 5)
				* (_c.a55 * af::cos(5.0 * (phi - _c.phi55)) + _c.a53 * af::cos(3.0 * (phi - _c.phi53)) + _c.a51 * af::cos(phi - _c.phi51)) / 5.0;
	if ((_c.a66 > 0) || (_c.a64 > 0) || (_c.a62 = 0) || (_c.C5 != 0))
		chi += af::pow(k0, 6)
				* (_c.a66 * af::cos(6.0 * (phi - _c.phi66)) + _c.a64 * af::cos(4.0 * (phi - _c.phi64)) + _c.a62 * af::cos(2.0 * (phi - _c.phi62))
						+ _c.C5) / 6.0;
//	p->Save2DDataSet(chi,"chi0");
	chi *= 2 * PI / _lambdaAngstrom;
//	p->Save2DDataSet(chi,"chi");
	return af::complex(af::cos(chi), -af::sin(chi));
}
} /* namespace slicepp */
