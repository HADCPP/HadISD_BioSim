#pragma once


#include <valarray>
#include <string>
#include <dlib/matrix.h>


typedef std::valarray<size_t> varraysize;
typedef std::valarray<float> varrayfloat;
typedef std::valarray<int> varrayInt;
typedef dlib::matrix<float> matrice;

typedef dlib::matrix<size_t > ind_matrice;
typedef std::vector<ind_matrice> vind_matrice;

struct _test
{
	bool duplicate = false;
	bool odd = false;
	bool frequent = false;
	bool diurnal = false;
	bool gap = true;
	bool records = true;
	bool streaks = true;
	bool climatological = true;
	bool spike = true;
	bool humidity = true;
	bool cloud = true;
	bool variance = true;
	bool winds = true;
};
typedef struct _test test;