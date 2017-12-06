/*
Repeated Streaks Check (RSC)
*
*   Checks for replication of
*     1) checks for consecutive repeating values
*     2) checks if one year has more repeating strings than expected
*     3) checks for repeats at a given hour across a number of days
*     4) checks for repeats for whole days - all 24 hourly values
*
*
*   Some thresholds now determined dynamically
*/

#pragma once
#include "utils.h"
#include "station.h"
#include "python_function.h"
#include "Utilities.h"


#include <boost/date_time/gregorian/gregorian.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <array>


namespace INTERNAL_CHECKS
{
	/*Derive threshold number for strings / streaks of repeating values

		: param object st_var : station variable object
		: param datetime start : start of data
		: param datetime end : end of data
		: param float reporting : reporting accuracy
		: param bool diagnostics : do diagnostic output
		: param bool plots : do plots
		: param float old_threshold : old threshold to use as comparison
	*/
	inline float rsc_get_straight_string_threshold(CMetVar& st_var, boost::gregorian::date   start, boost::gregorian::date  end, float reporting = 0., float old_threshold = 0.);

	/*
	Check for strings/streaks of repeating values
    :param object st_var: station variable object
    :param int n_days: number of days to exceed
    :param int n_obs: number of observations to exceed
    :param datetime start: start of data
    :param datetime end: end of data    
    :param float reporting: reporting accuracy
    :param bool wind: whether there is wind data to account for - extra minimum value
    :param bool dynamic: calculate threshold of number of observations dynamically rather than using n_obs
	*/
	varrayfloat rsc_straight_strings(CMetVar& st_var, const varrayfloat& times, int n_obs, int n_days, boost::gregorian::date  start, boost::gregorian::date end, std::map<float, float> WIND_MIN_VALUE, bool wind = false, float reporting = 0., bool dynamic = true);
	

	/*
	Look through for repeats of full 24h of data
     :param array data: data to process
    :param int n_wday: number of whole days to exceed
    :param MetVar st_var: station variable - mainly for plotting
    */
	void rsc_whole_day_repeats(std::vector<CMaskedArray<float>>& data, int n_wday, varrayfloat& flags);
	
	/*
	Repeat of same value at given hour for >N days
	:param object st_var: station variable object
	:param array times: timestamps
	:param int n_hrs: number of hours to exceed
	:param int n_wday: number of whole days to exceed (passed on)

	*/
	void rsc_hourly_repeats(CMetVar& st_var, const varrayfloat& times, int n_hrs, int n_wdays, varrayfloat& flags, varrayfloat&  day_flags);

	/*
	Find years where have more strings than expected, but not long enough to set off test
    
    :param array all_filtered: data filtered by all flags set so far
    :param array value_starts: locations of start of strings/streaks of data
    :param array value_lengths: lengths of each streak
    :param array flags: array of flags to be set
    :param datetime start: start of data
    :param datetime end: end of data    
    :param bool diagnostics: do diagnostic output
    :param bool plots: do plots
	return flags
	*/

	void rsc_annual_string_expectance(const CMaskedArray<float>& all_filtered, std::vector<int>& value_starts, std::vector<int>& value_lengths, varrayfloat& flags, boost::gregorian::date  start, boost::gregorian::date end, CMetVar& st_var, const  varrayfloat&times);
	
	void rsc(CStation& station, std::vector<std::string> var_list, std::vector<std::vector<int>>  flag_col, boost::gregorian::date  start,
		boost::gregorian::date end, std::ofstream& logfile);
}