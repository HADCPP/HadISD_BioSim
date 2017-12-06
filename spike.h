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

namespace INTERNAL_CHECKS
{

	/*
		Spike Check, looks for spikes up to 3 observations long, using thresholds
		calculated from the data itself.

		:param MetVar station: the station object
		:param list variable_list: list of observational variables to process
		:param list flag_col: the columns to set on the QC flag array
		:param datetime start: dataset start time
		:param datetime end: dataset end time
		:param file logfile: logfile to store outputs
		:param bool second: run for second time
	*/
	void sc(CStation &station, std::vector<std::string > variable_list, std::vector<int> flag_col, boost::gregorian::date  start,
		boost::gregorian::date end, std::ofstream&  logfile, bool second = false);

}