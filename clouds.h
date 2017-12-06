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


namespace INTERNAL_CHECKS
{

	/*
	Cloud observation code given as unobservable(== 9 or 10)

		:param obj station : station object
		: param list flag_col : flag columns to use
		: param file logfile : logfile to store outpu
		: param bool plots : to do any plots
		: param bool diagnostics : to do any extra diagnostic output
	*/
	void unobservable(CStation& station, std::vector<int> flag_col, std::ofstream& logfile);

	/*Total cloud cover less than maximum of low, mid and high

		: param obj station : station object
		: param list flag_col : flag columns to use
		: param file logfile : logfile to store output
		:
	*/
	void total_lt_max(CStation& station, std::vector<int> flag_col, std::ofstream& logfile);
	/*
		Call the logical cloud checks

	: param obj station : station object
	: param list flag_col : flag columns to use
	: param file logfile : logfile to store output

	*/
	/*
		Low cloud full, but values in mid or high
		:param obj station: station object
		:param list flag_col: flag columns to use
		:param file logfile: logfile to store output
   */
	void low_full(CStation& station, int flag_col, std::ofstream& logfile);
	/*
		Mid cloud full, but values in high
    
		:param obj station: station object
		:param list flag_col: flag columns to use
		:param file logfile: logfile to store output
	*/
	void mid_full(CStation& station, int flag_col, std::ofstream& logfile);

	/*
	If cloud base is 22000ft, then set to missing
    :param obj station: station object
	*/
	void fix_cloud_base(CStation& station);

	/*Non-sensical cloud value
    
    :param obj station: station object
    :param list flag_col: flag columns to use
    :param file logfile: logfile to store output*/
	void negative_cloud(CStation& station, int flag_col, std::ofstream& logfile);

	void ccc(CStation& station, std::vector<int> flag_col, std::ofstream& logfile);
		
}