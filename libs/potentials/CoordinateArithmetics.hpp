/*
 * CoordinateArithmetics.hpp
 *
 *  Created on: Aug 3, 2015
 *      Author: wenxuan
 */

#ifndef COORDINATEARITHMETICS_HPP_
#define COORDINATEARITHMETICS_HPP_
#define iwCoord(s,i,m) if ((i) > (m)/2) {(s) = (i) - (m);} else {(s) = (i);}

#define owCoord(t,i,m) ((t) = ((i) - ((m)/2))) //m >> 1;

#define iwCoordIp(i,m) if ((i) > (m)/2) {(i) -= (m);}

#define owCoordIp(i,m) ((i) -= ((m)/2))

#define dbCoord(i1, i2, j, m1) ((i1) = ((j) % (m1))); ((i2) = ((j) - (i1)) /(m1))

#define trCoord( i1, i2, i3, j, m1, m2 ) ( (i1) = ( ( (j) % ( (m1) * (m2) ) ) % (m1) ) ); ( (i2) = ( ( ( (j) % ( (m1) * (m2) ) ) - (i1) ) / (m1) ) ); ( (i3) = ( ( (j) - (i1) - ( (m1) * (i2) ) ) / ( (m1) * (m2) ) ) )

#define sgCoord(j, i1, i2, m1) ((j) = (((i2) * (m1)) + (i1)))

#define sgCoord3D(j, i1, i2, i3, m1, m2) ( (j) = ( ( (m1)*(m2)*(i3) ) + ( (m1)*(i2) ) + (i1) ) )

#endif /* COORDINATEARITHMETICS_HPP_ */
