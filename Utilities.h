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
