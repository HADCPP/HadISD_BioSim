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
	float SSS_MONTH_FRACTION = float(0.2);
}
namespace INTERNAL_CHECKS
{
	/*
		Supersaturation check, on individual obs, and then if >20 % of month affected

		: param array T : temperatures
		: param array D : dewpoint temperatures
		: param array month_ranges : array of month start and end times
		: param datetime start : DATASTART(for plotting)
		: param file logfile : logfile to store outputs
		: param bool plots : do plots or not
		: param bool diagnostics : extra verbose output

		: returns : flags - locations where flags have been set
	*/
	void hcc_sss(varrayfloat& T, varrayfloat& D, std::vector<std::pair<int, int>> month_ranges, std::ofstream& logfile, varrayfloat& flags);
		
	/*
	Dew point Depression check.  If long string of DPD = 0, then flag
    
    :param array times: timestamps
    :param array T: temperatures
    :param array D: dewpoint temperatures
    :param array P: precipitation depth
    :param array C: cloud base ??
    :param array SX: past significant weather ??
    :param datetime start: DATASTART (for plotting)
    :param file logfile: logfile to store outputs
    :param bool plots: do plots or not
    :param bool diagnostics: extra verbose output
*/
	varrayfloat  hcc_dpd(varrayfloat& times, CMaskedArray<float>& T, CMaskedArray<float>& D, CMaskedArray<float>& P, CMaskedArray<float>& C, CMaskedArray<float>& SX, std::ofstream& logfile);
	/*
	Check each month to see if most T obs have corresponding D obs
    do in bins of 10C - if one bin has <50% of match, then remove month
    
    :param array T: temperatures
    :param array D: dewpoint temperatures
    :param array month_ranges: array of month start and end times
    :param file logfile: logfile to store outputs
        :returns: flags - locations where flags have been set
	*/

	varrayfloat hcc_cutoffs(CMaskedArray<float>& T, CMaskedArray<float>& D, std::vector<std::pair<int, int>> month_ranges, std::ofstream& logfile);
	

	/*
		 Run humidity cross checks on temperature and dewpoint temperature obs

		:param object station: station object
		:param array flag_col: which columns to fill in flag_array
		:param datetime start: start of dataset
		:param datetime end: end of dataset
		:param file logfile: logfile to store outputs
  

	*/
		void hcc(CStation& station, std::vector<int> flag_col, boost::gregorian::date start, boost::gregorian::date  end, std::ofstream& logfile);

}