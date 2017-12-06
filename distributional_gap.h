#pragma once

#include "station.h"
#include "utils.h"
#include "python_function.h"
#include "Utilities.h"
#include "Statistic.h" 

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace
{
	 /// need to explain all of these
	 int OBS_LIMIT = 50;
	 int VALID_MONTHS = 5;
	 int SPREAD_LIMIT = 2;
	 int MONTH_LIMIT = 60;
	 int LARGE_LIMIT = 5; // # IQR
	 bool MEAN = false;
	 bool NON_GAP_FLAG = false;
	 double FREQUENCY_THRESHOLD = 0.1;
	 int GAP_SIZE = 2;
	 double	BIN_SIZE = 0.25;

}
namespace INTERNAL_CHECKS
{

	/*
	Walk the bins of the distribution to find a gap and return where it starts

	:param array hist: histogram values
	:param array bins: bin values
	:param flt threshold: limiting value
	:param int gap_size: gap size to record
	:returns:
	flt: gap_start
	*/
	int dgc_find_gap(varrayfloat& hist, varrayfloat& bins, float threshold, int gap_size = GAP_SIZE);
	varrayfloat dgc_all_obs(CStation& station, std::string variable, varrayfloat& flags, boost::gregorian::date start, boost::gregorian::date end, bool windspeeds = false, bool GH = false);
	float dgc_get_monthly_averages(CMaskedArray<float>&  data, int limit, double mdi, bool MEAN );
	float dgc_get_monthly_averages(varrayfloat&  Compress, int limit, double mdi, bool MEAN );
	/*
	'''
		Original Distributional Gap Check

		:param obj CStation: CStation object
		:param str variable: variable to act on
		:param array flags: flags array
		:param datetime start: data start
		:param datetime end: data end
		:param bool idl: run IDL equivalent routines for median
		:returns:
		flags - updated flag array
	'''

	*/
	void dgc_monthly(CStation& station, std::string variable, varrayfloat &flags, boost::gregorian::date start, boost::gregorian::date end);

	//'''Controller for two individual tests'''
	void dgc(CStation& station, std::vector<std::string> variable_list, std::vector<int> flag_col, boost::gregorian::date start,
		boost::gregorian::date  end, std::ofstream& logfile, bool GH = false);
}