//
// Created by toshiba on 24/10/2017.
//

#ifndef PASIVELOCATION_COORDINATES_SYSTEM_H
#define PASIVELOCATION_COORDINATES_SYSTEM_H

struct bearing {
	long double d;
	long double ph;
	long double z;
	
	bearing(long double d, long double ph, long double z) : d(d), ph(ph), z(z) {}
	bearing() : d(0), ph(0), z(0) {}
	
	bearing operator=(const bearing v)
	{
		this->d = v.d, this->ph = v.ph, this->z = v.z;
		return *this;
	}
	
	bearing operator-(const bearing v)
	{
		this->d -= v.d, this->ph -= v.ph, this->z -= v.z;
		return *this;
	}
	
	bearing operator+(const bearing v)
	{
		this->d += v.d, this->ph += v.ph, this->z += v.z;
		return *this;
	}

//	bearing& operator = (const bearing &v);
//	bearing& operator - (const bearing &v);
//	bearing& operator + (const bearing &v);
};

struct cartesian_normal {
	long double x;
	long double y;
	long double z;
	
	cartesian_normal(long double x, long double y, long double z) : x(x), y(y), z(z) {}
	cartesian_normal() : x(0), y(0), z(0) {}
	
	cartesian_normal operator=(cartesian_normal v)
	{
		this->x = v.x, this->y = v.y, this->z = v.z;
		return *this;
	}
	
	cartesian_normal operator-(cartesian_normal v)
	{
		this->x -= v.x, this->y -= v.y, this->z -= v.z;
		return *this;
	}
	
	cartesian_normal operator+(cartesian_normal v)
	{
		this->x += v.x, this->y += v.y, this->z += v.z;
		return *this;
	}

//	cartesian_normal& operator = (const cartesian_normal &v);
//	cartesian_normal& operator - (const cartesian_normal &v);
//	cartesian_normal& operator + (const cartesian_normal &v);
}; //src/coordinates_system.cpp

struct cartesian_special {
	long double x;
	long double y;
	long double z;
	
	cartesian_special(long double x, long double y, long double z) : x(x), y(y), z(z) {}
	cartesian_special() : x(0), y(0), z(0) {}
	
	cartesian_special operator=(cartesian_special v)
	{
		this->x = v.x, this->y = v.y, this->z = v.z;
		return *this;
	}
	
	cartesian_special operator-(cartesian_special v)
	{
		this->x -= v.x, this->y -= v.y, this->z -= v.z;
		return *this;
	}
	
	cartesian_special operator+(cartesian_special v)
	{
		this->x += v.x, this->y += v.y, this->z += v.z;
		return *this;
	}

//	cartesian_special& operator = (const cartesian_special &v);
//	cartesian_special& operator - (const cartesian_special &v);
//	cartesian_special& operator + (const cartesian_special &v);
};

#endif //PASIVELOCATION_COORDINATES_SYSTEM_H