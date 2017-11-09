//
// Created by toshiba on 24/10/2017.
//

#ifndef PASIVELOCATION_PASIVE_LOCATION_H
#define PASIVELOCATION_PASIVE_LOCATION_H

#define _USE_MATH_DEFINES
//#include "targetver.h"
#include <math.h>
//#include <cmath>
#include <complex>      // std::complex
#include "coordinates_system.h"
//#include "coordinates_system.h"
#include <iostream>
//#include <stdio.h>
//#include <tchar.h>

#define N 4
#define PI		3.14159265358979323846

class Pasive_Location {
private:
	long double alpha;
	long double beta;
	long double sigma;
	long double NoiseError;
	bearing* station_b;           // Bearing
	cartesian_normal* station_n;  // Cartesian_normal
	cartesian_special* station_s; // Cartesian_special
	long double var_0, var_1, var_2, var_3, var_4, var_5, var_6, var_7, var_8, var_9;
	
	void Turning_The_System();
	long double* Simulation(cartesian_special);
	cartesian_special* Solve(long double[3]);
public:
	Pasive_Location(bearing[N]);
	Pasive_Location(cartesian_normal[N]);
	
	void Change_Base(bearing[N]);
	void Change_Base(cartesian_normal[N]);
	
	bearing* Solve_Simulation(bearing, const uint8_t);
	cartesian_normal* Solve_Simulation(cartesian_normal, const uint8_t);
	cartesian_special* Solve_Simulation(cartesian_special);
	
	bearing* Solve_b(long double[3], const uint8_t);
	cartesian_normal* Solve_n(long double[3], const uint8_t);
	//cartesian_special* Solve(long double*, unsigned char);
	
	cartesian_special Bearing_To_Cartesian_special(bearing);
	cartesian_normal Bearing_To_Cartesian_normal(bearing);
	bearing Cartesian_special_To_Bearing(cartesian_special);
	bearing Cartesian_normal_To_Bearing(cartesian_normal);
	cartesian_special Cartesian_normal_To_Cartesian_special(cartesian_normal);
	cartesian_normal Cartesian_special_To_Cartesian_normal(cartesian_special);
};

#endif //PASIVELOCATION_PASIVE_LOCATION_H
