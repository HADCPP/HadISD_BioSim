#pragma once
#include "station.h"
#include "Utilities.h"
#include "python_function.h"
#include "utils.h"

#include <vector>
#include <string>
#include <iostream>
#include <valarray>
#include <boost/math/constants/constants.hpp>


namespace
{
	const int OBS_PER_DAY = 4;
	const int DAILY_RANGE = 5;
	const int HOURS = 24;
	
	bool DYNAMIC_DIURNAL = true;
}
namespace INTERNAL_CHECKS
{
	const int INTMDI = -999;

	inline float __cdecl MyApplySinus(float n)
	{
		return std::sin(2 * boost::math::constants::pi<float>()*n);
	}
	
	//Return sine curve over 24 points
	inline std::valarray<float> dcc_make_sine()
	{
		std::valarray<float> val_sin = PYTHON_FUNCTION::arange<float>(float(HOURS));
		val_sin /= 24;
		return val_sin = val_sin.apply(MyApplySinus);
	}
	/*
	Check if >=3 quartiles of the day have data

	:param array day: 24 hours of data in array
	:returns: boolean
	*/
	bool dcc_quartile_check(CMaskedArray<float>& day);
	/*
			The diurnal cycle check.

		:param object station : the station object to be processed
		: param list variable_list : the variables to be processed
		: param list full_variable_list : the variables for flags to be applied to
		: param list flag_col : which column in the qc_flags array to work on
		: param file logfile : logfile to store outputs
		: returns :
		*/
	void dcc(CStation& station, std::vector<std::string> variable_list, std::vector<std::string>& full_variable_list, std::vector<int> flag_col,
		std::ofstream &logfile);
}