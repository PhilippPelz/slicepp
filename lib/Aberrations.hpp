/*
 * Aberrations.hpp
 *
 *  Created on: Jan 13, 2016
 *      Author: philipp
 */

#ifndef ABERRATIONS_HPP_
#define ABERRATIONS_HPP_
#include "stemtypes_fftw3.hpp"
#include <arrayfire.h>
#include "config_IO/configs.hpp"
#include "data_IO/PersistenceManager.hpp"
namespace slicepp {

class Aberrations {
public:
	Aberrations(float_tt lambda,AberrationConfig c,std::string name);
	virtual ~Aberrations();
	af::array getPhasePlate(af::array& k,af::array& kx,af::array& ky, PersistenceManagerPtr p);
	void DisplayParams();
protected:
	float_tt _lambdaAngstrom;
	AberrationConfig _c;
	std::string _name;
};

} /* namespace slicepp */

#endif /* ABERRATIONS_HPP_ */
