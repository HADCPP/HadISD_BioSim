#pragma once
#include "station.h"
#include "utils.h"
#include "utils.h"

#include <vector>
#include <string>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <iostream>

namespace INTERNAL_CHECKS
{
	void is_month_duplicated(const std::valarray<float>& source_data, const std::valarray<float>& target_data, const varraysize& valid, int sm,
		int tm, std::valarray<int>& duplicated);

	void duplication_test(const std::valarray<float>& source_data, const std::valarray<float>& target_data, const varraysize& valid, int sm, int tm,
		std::map<int, int>::const_iterator source_month, std::map<int, int>::const_iterator  target_month, std::vector<int> &duplicated, CStation& station, int flag_col);

	void dmc(CStation& station, std::vector<std::string>& variable_list, std::vector<std::string>& full_variable_list, int flag_col,
		boost::gregorian::date start,boost::gregorian::date end, std::ofstream &logfile);

	
}