/*
 * ConfigReader.hpp
 *
 *  Created on: Dec 14, 2015
 *      Author: philipp
 */

#ifndef CONFIGREADER_HPP_
#define CONFIGREADER_HPP_

#include "configs.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace QSTEM {

typedef struct ConfigReader {
public:
	ConfigReader();
	ConfigPtr Read(boost::filesystem::path configPath);
} ConfigReader;

} /* namespace QSTEM */

#endif /* CONFIGREADER_HPP_ */
