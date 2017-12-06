#pragma once
#include "station.h"
#include "python_function.h"
#include "utils.h"
#include "Utilities.h"
#include "Statistic.h" 
#include "distributional_gap.h"

#include <boost/date_time/gregorian/gregorian.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <map>


namespace INTERNAL_CHECKS
{
	
	/*
		Get the weights for the low pass filter.
		:param array monthly_vqvs : monthly anomalies
		: param array monthly_subset : which values to take
		: param array filter_subset : which values to take
		: returns :
		weights
	*/

	float coc_get_weights(const CMaskedArray<float>& monthly_vqvs, varraysize monthly_subset, varraysize filter_subset);
		

	/*
			Run the low pass filter - get suitable ranges, get weights, and apply

			: param array normed_anomalies : input normalised anomalies
			: param array year_ids : look up array for which row is which year
			: param array monthly_vqvs : monthly average anomalies
			: param int years : number of years in data
			: returns :
			normed_anomalies - with weights applied
	*/
	void coc_low_pass_filter(std::valarray<CMaskedArray<float>>& normed_anomalies,const  std::vector<int>& year_ids,const CMaskedArray<float>& monthly_vqvs, int years);
	
	/*
	 Find which values need flagging, do so

    :param array month_ranges: start and end locations for each month
    :param array normed_anomalies: data to flag
    :param array flags: flag array
    :param array year_ids: year look-up array
    :param int threshold: critical value to start flagging
    :param int gap_start: location of gap if exists, else 0
    :param bool upper: upper or lower part of distribution
    :param bool plots: do plots


    :returns:
        flags - updated flag array
       */
	void coc_find_and_apply_flags(std::valarray<std::pair<int, int>>& month_ranges, std::valarray<CMaskedArray<float>>& normed_anomalies, varrayfloat& flags, const std::vector<int>& year_ids, float threshold, int gap_start, bool upper = true);
	void coc(CStation& station, std::vector<std::string> variable_list, std::vector<int> flag_col, boost::gregorian::date  start,
		boost::gregorian::date end, std::ofstream&  logfile, bool idl = false);
}