
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
#include <numeric>

namespace
{
	const float MAD_THRESHOLD = 4.;
	const int DATA_COUNT_THRESHOLD = 120;

}
namespace INTERNAL_CHECKS
{
	void evc(CStation& station, std::vector<std::string> variable_list, std::vector<int> flag_col, boost::gregorian::date start, boost::gregorian::date  end,
		std::ofstream& logfile, bool idl = false);
}