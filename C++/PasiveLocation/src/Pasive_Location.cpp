//
// Created by Vova Dumansky 24.10.2017.
//

#include "Pasive_Location.h"

Pasive_Location::Pasive_Location(bearing bear[N]) : station_b(bear) {
	//cartesian_normal base_n[N];
	station_n = (cartesian_normal*)malloc(sizeof(cartesian_normal) * N); //base_n;
	station_n[0] = Bearing_To_Cartesian_normal(station_b[0]);
	for (int i = 1; i < N; i++)
	{
		station_n[i] = Bearing_To_Cartesian_normal(station_b[i]) - this->station_n[0];
	}
	Turning_The_System();
}

Pasive_Location::Pasive_Location(cartesian_normal cart[N]) : station_n(cart) {
	//bearing base_b[N];
	//station_b = base_b;
	station_b = (bearing*)malloc(sizeof(bearing) * N);
	//station_n[N] = cart[N];
	for (int i = 0; i < N; i++)
	{
		station_b[i] = Cartesian_normal_To_Bearing(station_n[i]);
		if (i != 0)
			station_n[i] = station_n[i] - station_n[0];
	}
	Turning_The_System();
}

void Pasive_Location::Change_Base(bearing station_b[N]) {
	this->station_b = station_b;
	//station_n[0] = cartesian_normal();
	station_n[0].x = 0; station_n[0].y = 0; station_n[0].z = 0;
	for (int i = 0; i < N; i++)
	{
		station_n[i] = Bearing_To_Cartesian_normal(station_b[i]) - station_n[0];
	}
	Turning_The_System();
}

void Pasive_Location::Change_Base(cartesian_normal station_n[N]) {
	//this->station_n = station_n;
	this->station_n[0].x = 0; this->station_n[0].y = 0; this->station_n[0].z = 0;
	for (int i = 0; i < N; i++)
	{
		station_b[i] = Cartesian_normal_To_Bearing(station_n[i]);
		this->station_n[i] = station_n[i] - this->station_n[0];
	}
	Turning_The_System();
}

void Pasive_Location::Turning_The_System() {
	cartesian_normal v[N];
	//cartesian_special base_s[N];
	//station_s = base_s;
	station_s = (cartesian_special*)malloc(sizeof(cartesian_special) * N);
	
	this->alpha = atan(station_n[1].y / station_n[1].x);
	for (int i = 1; i < N; i++)
	{
		v[i].x = station_n[i].x * cos(alpha) + station_n[i].y * sin(alpha);
		v[i].y = -station_n[i].x * sin(alpha) + station_n[i].y * cos(alpha); // v[1].y = 0;
		//v[i].z = station_n[i].z;
	}
	
	this->beta = atan(station_n[1].z / v[1].x);
	for (int i = 1; i < N; i++)
	{
		station_s[i].x = v[i].x * cos(beta) + station_n[1].z * sin(beta);
		//station_s[i].y = v[i].y;
		//station_s[i].z = -v[i].x * sin(beta) + v[i].z * cos(beta); // station_s[1].z = 0;
		//v[i].z = station_n[i].z;
		v[i].z = -v[i].x * sin(beta) + station_n[1].z * cos(beta); // v[1].z = 0;
	}
	station_s[1].z = 0;
	station_s[1].y = 0;
	
	this->sigma = atan(v[2].z / v[2].y);
	for (int i = 2; i < N; i++)
	{
		station_s[i].y = v[i].y * cos(sigma) + v[i].z * sin(sigma);
		station_s[i].z = -v[i].y * sin(sigma) + v[i].z * cos(sigma); // v[2].z = 0;
	}
	station_s[2].z = 0;
	
	
	long double x12 = pow(station_s[1].x, 2);
	//long double x13 = pow(station_s[1].x, 3);
	
	long double x22 = pow(station_s[2].x, 2);
	long double x23 = pow(station_s[2].x, 3);
	//long double x24 = pow(station_s[2].x, 4);
	long double y22 = pow(station_s[2].y, 2);
	long double y23 = pow(station_s[2].y, 3);
	//long double y24 = pow(station_s[2].y, 4);
	
	long double x32 = pow(station_s[3].x, 2);
	//long double x33 = pow(station_s[3].x, 3);
	long double y32 = pow(station_s[3].y, 2);
	long double y33 = pow(station_s[3].y, 3); //5
	long double z32 = pow(station_s[3].z, 2);
	//long double z33 = pow(station_s[3].z, 3);
	
	
	var_1 = ( ( station_s[2].y*station_s[3].y - x32 + station_s[1].x*station_s[3].x - y32
	            + (x22*station_s[3].y - station_s[1].x*station_s[2].x*station_s[3].y)/(station_s[2].y) )/( z32 ) - 1 );
	var_2 = ( station_s[1].x*station_s[2].x*y32-x22*z32-x22*y32+station_s[1].x*station_s[2].x*z32+x32*station_s[2].y*station_s[3].y
	          - station_s[1].x*station_s[3].x*station_s[2].y*station_s[3].y-y22*y32-y22*z32+station_s[2].y*y33+station_s[2].y*station_s[3].y*z32 )/( y22*z32 );
	var_3 = ( x22*y32+x22*z32-2*station_s[2].x*station_s[3].x*station_s[2].y*station_s[3].y+x32*y22+y22*z32 )/( x12*y22*z32 );
	var_4 = ( station_s[2].x*y32-station_s[3].x*station_s[2].y*station_s[3].y+station_s[2].x*z32 )/( station_s[1].x*y22*z32 );
	var_5 = ( station_s[2].x*station_s[3].y-station_s[3].x*station_s[2].y )/( station_s[1].x*station_s[2].y*z32 );
	var_6 = ( y32 + z32 )/( y22*z32 );
	var_7 = ( ( x23*y32+x23*z32-x22*station_s[3].x*station_s[2].y*station_s[3].y-station_s[2].x*x32*station_s[2].y*station_s[3].y+station_s[2].x*y22*y32
	            + station_s[2].x*y22*z32-station_s[2].x*station_s[2].y*y33-station_s[2].x*station_s[2].y*station_s[3].y*z32+pow(station_s[3].x, 3)*y22
	            - station_s[3].x*y23*station_s[3].y+station_s[3].x*y22*y32+station_s[3].x*y22*z32 )/( station_s[1].x*y22*z32 )
	          - ( x22*y32+x22*z32-2*station_s[2].x*station_s[3].x*station_s[2].y*station_s[3].y+x32*y22 )/( y22*z32 ) - 1 );
	var_8 = ( station_s[3].y )/( station_s[2].y*z32 );
	var_9 = z32; // 1/var_9
	var_0 = ( x12*x22*y32+x12*x22*z32-2*x12*station_s[2].x*station_s[3].x*station_s[2].y*station_s[3].y+x12*x32*y22+x12*y22*z32
	          - 2*station_s[1].x*x23*y32-2*station_s[1].x*x23*z32+2*station_s[1].x*x22*station_s[3].x*station_s[2].y*station_s[3].y
	          + 2*station_s[1].x*station_s[2].x*x32*station_s[2].y*station_s[3].y-2*station_s[1].x*station_s[2].x*y22*y32
	          - 2*station_s[1].x*station_s[2].x*y22*z32+2*station_s[1].x*station_s[2].x*station_s[2].y*y33
	          + 2*station_s[1].x*station_s[2].x*station_s[2].y*station_s[3].y*z32-2*station_s[1].x*pow(station_s[3].x, 3)*y22
	          + 2*station_s[1].x*station_s[3].x*y23*station_s[3].y-2*station_s[1].x*station_s[3].x*y22*y32
	          - 2*station_s[1].x*station_s[3].x*y22*z32+pow(station_s[2].x, 4)*y32+pow(station_s[2].x, 4)*z32-2*x22*x32*station_s[2].y*station_s[3].y
	          + 2*x22*y22*y32+2*x22*y22*z32-2*x22*station_s[2].y*y33-2*x22*station_s[2].y*station_s[3].y*z32
	          + pow(station_s[3].x, 4)*y22-2*x32*y23*station_s[3].y+2*x32*y22*y32+2*x32*y22*z32+pow(station_s[2].y, 4)*y32
	          + pow(station_s[2].y, 4)*z32-2*y23*y33-2*y23*station_s[3].y*z32+y22*pow(station_s[3].y, 4)+2*y22*y32*z32+y22*pow(station_s[3].z, 4) ) / ( 4*y22*z32 );
	
}

cartesian_special Pasive_Location::Bearing_To_Cartesian_special(bearing bear) {
	return Cartesian_normal_To_Cartesian_special(
		Bearing_To_Cartesian_normal(bear) );
}

cartesian_normal Pasive_Location::Bearing_To_Cartesian_normal(bearing bear) {
	long double c;
	try
	{
		c = sqrt(pow(bear.d, 2) - pow(bear.z, 2));
		if (isnan(c)) {
			std::cout << "Problem using function.  Assigning a value of 0.\n";
			c = 0;
		}
	}
	catch (int i)//???? ?????????? ????? 123
	{
		std::cout << "Error #" << i << " - Problem using function.  Assigning a value of 0.\n";
		c = 0;
	}
	return cartesian_normal( c * cos(bear.ph), -c * sin(bear.ph), bear.z );
}

bearing Pasive_Location::Cartesian_special_To_Bearing(cartesian_special cart) {
	return Cartesian_normal_To_Bearing(
		Cartesian_special_To_Cartesian_normal(cart) );
}

bearing Pasive_Location::Cartesian_normal_To_Bearing(cartesian_normal cart) {
	return bearing( sqrt(pow(cart.x , 2) + pow(cart.y , 2) + pow(cart.z , 2)),
	                atan(-cart.y / cart.x) + PI * (cart.x < 0), cart.z );
}

cartesian_special Pasive_Location::Cartesian_normal_To_Cartesian_special(cartesian_normal cart) {
	cartesian_normal v = cart - station_n[0];
	
	cart.x = v.x * cos(alpha) + v.y * sin(alpha);
	cart.y = -v.x * sin(alpha) + v.y * cos(alpha);
	//cart.z = v.z;
	
	//v.x = cart.x * cos(beta) + v.z * sin(beta);
	//v.y = cart.y;
	cart.z = -cart.x * sin(beta) + v.z * cos(beta);
	
	return cartesian_special( cart.x * cos(beta) + v.z * sin(beta),
	                          cart.y * cos(sigma) + cart.z * sin(sigma),
	                          -cart.y * sin(sigma) + cart.z * cos(sigma) ); // (cartesian_special)cart
}

cartesian_normal Pasive_Location::Cartesian_special_To_Cartesian_normal(cartesian_special cart) {
	cartesian_special v;
	
	//cart.x;
	v.y = cart.y * cos(sigma) - cart.z * sin(sigma);
	v.z = cart.y * sin(sigma) + cart.z * cos(sigma);
	
	v.x = cart.x * cos(beta) - v.z * sin(beta);
	// = v.y;
	//cart.z = cart.x * sin(beta) + v.z * cos(beta);
	
	//v.x * cos(alpha) - v.y * sin(alpha) + station_n[0].x;
	//v.x * sin(alpha) + v.y * cos(alpha) + station_n[0].y;
	//cart.z + station_n[0].z;
	
	return cartesian_normal(v.x * cos(alpha) - v.y * sin(alpha) + station_n[0].x,
	                        v.x * sin(alpha) + v.y * cos(alpha) + station_n[0].y,
	                        cart.x * sin(beta) + v.z * cos(beta) + station_n[0].z );
}

bearing* Pasive_Location::Solve_Simulation (bearing goal, const uint8_t base = 0 ) {
	long double* dR = Simulation(Bearing_To_Cartesian_special(goal));
	return Solve_b(dR, base);
}

cartesian_normal* Pasive_Location::Solve_Simulation (cartesian_normal goal, const uint8_t base = 0) {
	long double *dR = Simulation(Cartesian_normal_To_Cartesian_special(goal));
	return Solve_n(dR, base);
}

cartesian_special* Pasive_Location::Solve_Simulation (cartesian_special goal) {
	long double* dR = Simulation(goal);
	return Solve(dR);
}

long double* Pasive_Location::Simulation(cartesian_special target) {
	long double* dR = (long double*)malloc(sizeof(long double) * N);
	dR[0] = sqrt( pow(target.x, 2) + pow(target.y, 2) + pow(target.z, 2) );
	for (int i = 1; i < N; ++i)
	{
		dR[i] = sqrt( pow(target.x - station_s[i].x, 2) + pow(target.y - station_s[i].y, 2) + pow(target.z - station_s[i].z, 2) ) - dR[0];
	}
	return &dR[1];
}

bearing* Pasive_Location::Solve_b(long double* d_radius, const uint8_t base = 0) {
	cartesian_normal* root_2 = Solve_n(d_radius, base);
	bearing* target = (bearing*)malloc(sizeof(bearing) * 2);
	target[0] = Cartesian_normal_To_Bearing(root_2[0]);
	target[1] = Cartesian_normal_To_Bearing(root_2[1]);
	return target;
}

cartesian_normal* Pasive_Location::Solve_n(long double* d_radius, const uint8_t base = 0) {
	cartesian_special* root_2 = Solve(d_radius);
	cartesian_normal* target = (cartesian_normal*)malloc(sizeof(cartesian_normal) * 2);
	target[0] = Cartesian_special_To_Cartesian_normal(root_2[0]);
	target[1] = Cartesian_special_To_Cartesian_normal(root_2[1]);
	if (base > 0 && base <= N) {
		target[0] = target[0] - station_n[base - 1];
		target[1] = target[1] - station_n[base - 1];
	}
	return target;
}

cartesian_special* Pasive_Location::Solve(long double d_radius[3]) {
	long double a = d_radius[2] * d_radius[2] / var_9
	                - 1
	                - d_radius[1] * d_radius[2] * var_8 * 2
	                + d_radius[1] * d_radius[1] * (var_6)
	                +d_radius[0] * d_radius[2] * (2 * var_5)
	                - d_radius[0] * d_radius[1] * (2 * var_4)
	                + d_radius[0] * d_radius[0] * var_3;
	
	long double b = d_radius[2] * var_1
	                + pow(d_radius[2], 3) / var_9
	                + d_radius[1] * var_2
	                - (d_radius[1] * pow(d_radius[2], 2) + pow(d_radius[1], 2)*d_radius[2])*var_8
	                + pow(d_radius[1], 3)*var_6
	                + d_radius[0] * var_7
	                + (d_radius[0] * pow(d_radius[2], 2) + pow(d_radius[0], 2)*d_radius[2])*var_5
	                - (d_radius[0] * pow(d_radius[1], 2) + pow(d_radius[0], 2)*d_radius[1])*var_4
	                + pow(d_radius[0], 3)*var_3;
	
	long double c = var_0
	                + pow(d_radius[2], 2)*var_1 / 2
	                + pow(d_radius[2], 4) / (var_9 * 4)
	                + pow(d_radius[1], 2)*var_2 / 2
	                - pow(d_radius[1], 2)*pow(d_radius[2], 2)*var_8 / 2
	                + pow(d_radius[1], 4)*var_6 / 4
	                + pow(d_radius[0], 2)*var_7 / 2
	                + pow(d_radius[0], 2)*pow(d_radius[2], 2)*var_5 / 2
	                - pow(d_radius[0], 2)*pow(d_radius[1], 2)*var_4 / 2
	                + pow(d_radius[0], 4)*var_3 / 4;
	
	long double d[3];
	d[2] = sqrt(pow(b, 2) - 4 * a*c);
	if (isnan(d[2])) {
		std::cout << "Problem using function.  Assigning a value of 0.\n";
		d[2] = 0;
	}
	d[0] = (-b + d[2]) / (2 * a);
	d[1] = (-b - d[2]) / (2 * a);
	
	/*
	long double x = -(  pow(d_radius[0], 2)-  pow(station_s[1].x, 2) +2*d[1]*d_radius[0] )/( 2*station_s[1].x );
	long double y = ( station_s[2].x* pow(d_radius[0], 2)+ 2*d[1]*station_s[2].x*d_radius[0] )/( 2*station_s[1].x*station_s[2].y ) -
	(  pow(d_radius[1], 2)-  pow(station_s[2].x, 2) - pow(station_s[2].y, 2) +station_s[1].x*station_s[2].x+2*d[1]*d_radius[1] )/( 2*station_s[2].y );
	long double z = (  pow(d_radius[0], 2)* station_s[3].x*station_s[2].y- pow(d_radius[0], 2)* station_s[2].x*station_s[3].y-2*d[1]*d_radius[0]*station_s[2].x*station_s[3].y+2*d[1]*d_radius[0]*station_s[3].x*station_s[2].y +
	 pow(d_radius[1], 2)* station_s[1].x*station_s[3].y+2*d[1]*d_radius[1]*station_s[1].x*station_s[3].y- pow(d_radius[2], 2)* station_s[1].x*station_s[2].y-2*d[1]*d_radius[2]*station_s[1].x*station_s[2].y +
	 pow(station_s[1].x, 2) *station_s[2].x*station_s[3].y- pow(station_s[1].x, 2) *station_s[3].x*station_s[2].y-station_s[1].x* pow(station_s[2].x, 2) *station_s[3].y+station_s[1].x* pow(station_s[3].x, 2) *station_s[2].y -
	station_s[1].x* pow(station_s[2].y, 2) *station_s[3].y+station_s[1].x*station_s[2].y* pow(station_s[3].y, 2) +station_s[1].x*station_s[2].y* pow(station_s[3].z, 2)  )/( 2*station_s[1].x*station_s[2].y*station_s[3].z );
	//*/
	//cartesian_special goal1[2] = { { 1,1,1 }, { 1,1,1 } };
	cartesian_special* goal = (cartesian_special*)malloc(sizeof(cartesian_special) * 2);
	goal[0] = {
		-(  pow(d_radius[0], 2)-  pow(station_s[1].x, 2) +2*d[0]*d_radius[0] )/( 2*station_s[1].x )
		,( station_s[2].x* pow(d_radius[0], 2)+ 2*d[0]*station_s[2].x*d_radius[0] )/( 2*station_s[1].x*station_s[2].y )
		 - (  pow(d_radius[1], 2)-  pow(station_s[2].x, 2) - pow(station_s[2].y, 2) +station_s[1].x*station_s[2].x
		      + 2*d[0]*d_radius[1] )/( 2*station_s[2].y )
		,(  pow(d_radius[0], 2)* station_s[3].x*station_s[2].y- pow(d_radius[0], 2)* station_s[2].x*station_s[3].y
		    - 2*d[0]*d_radius[0]*station_s[2].x*station_s[3].y+2*d[0]*d_radius[0]*station_s[3].x*station_s[2].y +
		    pow(d_radius[1], 2)* station_s[1].x*station_s[3].y+2*d[0]*d_radius[1]*station_s[1].x*station_s[3].y
		    - pow(d_radius[2], 2)* station_s[1].x*station_s[2].y-2*d[0]*d_radius[2]*station_s[1].x*station_s[2].y +
		    pow(station_s[1].x, 2) *station_s[2].x*station_s[3].y- pow(station_s[1].x, 2) * station_s[3].x*station_s[2].y
		    - station_s[1].x* pow(station_s[2].x, 2) *station_s[3].y+station_s[1].x* pow(station_s[3].x, 2) *station_s[2].y
		    - station_s[1].x* pow(station_s[2].y, 2) *station_s[3].y+station_s[1].x*station_s[2].y* pow(station_s[3].y, 2)
		    + station_s[1].x*station_s[2].y* pow(station_s[3].z, 2)  )/( 2*station_s[1].x*station_s[2].y*station_s[3].z )
		};
	goal[1] = {
		-(  pow(d_radius[0], 2)-  pow(station_s[1].x, 2) +2*d[1]*d_radius[0] )/( 2*station_s[1].x )
		,( station_s[2].x* pow(d_radius[0], 2)+ 2*d[1]*station_s[2].x*d_radius[0] )/( 2*station_s[1].x*station_s[2].y )
		 - (  pow(d_radius[1], 2)-  pow(station_s[2].x, 2) - pow(station_s[2].y, 2) +station_s[1].x*station_s[2].x
		      + 2*d[1]*d_radius[1] )/( 2*station_s[2].y )
		,(  pow(d_radius[0], 2)* station_s[3].x*station_s[2].y- pow(d_radius[0], 2)* station_s[2].x*station_s[3].y
		    - 2*d[1]*d_radius[0]*station_s[2].x*station_s[3].y+2*d[1]*d_radius[0]*station_s[3].x*station_s[2].y +
		    pow(d_radius[1], 2)* station_s[1].x*station_s[3].y+2*d[1]*d_radius[1]*station_s[1].x*station_s[3].y
		    - pow(d_radius[2], 2)* station_s[1].x*station_s[2].y-2*d[1]*d_radius[2]*station_s[1].x*station_s[2].y +
		    pow(station_s[1].x, 2) *station_s[2].x*station_s[3].y- pow(station_s[1].x, 2) * station_s[3].x*station_s[2].y
		    - station_s[1].x* pow(station_s[2].x, 2) *station_s[3].y+station_s[1].x* pow(station_s[3].x, 2) *station_s[2].y
		    - station_s[1].x* pow(station_s[2].y, 2) *station_s[3].y+station_s[1].x*station_s[2].y* pow(station_s[3].y, 2)
		    + station_s[1].x*station_s[2].y* pow(station_s[3].z, 2)  )/( 2*station_s[1].x*station_s[2].y*station_s[3].z )
	};
	
	//return &d[0];
	//cartesian_special goal[2] = {{1,1,1}, {1,2,2}} ;
	return goal;
}

//cartesian_special* Solve(long double*, unsigned char) {}