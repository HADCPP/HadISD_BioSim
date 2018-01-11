#pragma once


#include <dlib/matrix.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time.hpp>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime> 
#include <errno.h>
#include <exception>
#include <chrono>
#include <cmath>
#include <valarray>
#include <string>

typedef boost::numeric::ublas::vector<float> fvector;
typedef boost::numeric::ublas::vector<std::string> svector;
typedef boost::numeric::ublas::vector<int> ivector;

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
	bool cloud = false;
	bool variance = true;
	bool winds = true;
};
typedef struct _test test;

namespace
{
	const std::string LOG_OUTFILE_LOCS = "D:\\HadISD_BioSim\\Data\\Log\\";
	const std::string NETCDF_DATA_LOCS = "D:\\HadISD_BioSim\\Data\\NetCDF_files\\";
	const std::string CSV_OUTFILE_LOCS = "D:\\HadISD_BioSim\\Data\\Meteo_data\\East-Canada\\Quebec 2002-2015H\\";
	const std::string NEW_CSV_OUTFILE_LOCS = "D:\\HadISD_BioSim\\Data\\New_Meteo_data\\";
	const int INTMDI = -999;
	const float FLTMDI = -1e30;
	const int NBVAR = 69; // recuperer la taille des données dans le fichier netCDF
	std::vector<std::string> process_var{ "temperatures", "dewpoints", "windspeeds", "winddirs", "slp" };
	std::vector<std::string> carry_thru_vars;
	std::map<const char*, varraysize> FLAG_COL_DICT = { {"temperatures", { 0, 1, 4, 5, 8, 12, 16, 20, 24, 27, 41, 44, 54, 58 } },
	{"dewpoints", { 0, 2, 4, 6, 8, 9, 13, 17, 21, 25, 28, 30, 31, 32, 42, 45, 55, 59 } },
	{"slp", { 0, 3, 4, 7, 11, 15, 19, 23, 26, 29, 43, 46, 57, 60 } },
	{"windspeeds", { 0, 4, 10, 14, 18, 22, 47, 56, 61, 62, 63, 64, 65 } },
	{"winddirs", { 0, 4, 10, 14, 18, 22, 47, 48, 56, 61, 62, 63, 64, 65, 66, 67, 68 } } };
}

