#pragma once
#include "station.h"
#include "utils.h"
#include "python_function.h"
#include "Utilities.h"

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>

namespace
{

	int DEGREEBINS = 45;
	double PROB_THRESHOLD = 0.01;
}
namespace INTERNAL_CHECKS
{


	/*
	Select occurrences of wind speed and direction which are 
    logically inconsistent with measuring practices.

    :param Station station: station object
    :param array flag_col: which columns to use in QC flag array
    :param file logfile: logfile to output to
  */
	void logical_checks(CStation &station, std::vector<int> flag_col, std::ofstream& logfile);

	/*
		Checks for large differences in the year-to-year wind-rose shape.  
		Uses RMSE and fits Gaussian.  Finds gap in distribution to flag beyond

		:param  station: station object
		:param int flag_col: which column to store the flags in
		:param datetime start: start of data
		:param datetime end: end of data
		
	*/
	/*
	 create bins and bin centres from data 
		given bin width covers entire range'
	*/
	void wind_create_bins(CMaskedArray<float>&  indata, varrayfloat& binEdges, varrayfloat& bincenters);
	/*
	Adjust the normalisation of the fitted PDF to the histogram in natural units

	: param array indata : input data
	: param array bins : bins to use for histogram

	: returns : norm
	*/
	float get_histogram_norm(varrayfloat&  indata, varrayfloat& binEdges);
	/*
		 Checks for large differences in the year-to-year wind-rose shape.  
		Uses RMSE and fits Gaussian.  Finds gap in distribution to flag beyond

		:param MetStation station: station object
		:param int flag_col: which column to store the flags in
		:param datetime start: start of data
    
	*/
	
	void wind_rose_check(CStation &station, int flag_col, boost::gregorian::date start, boost::gregorian::date  end,
		std::ofstream& logfile);

	// Specific wind speed and direction checks
	void wdc(CStation& station, std::vector<int> flag_col, boost::gregorian::date start, boost::gregorian::date  end,
		std::ofstream& logfile);
}