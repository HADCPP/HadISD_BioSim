/** Known Records Check (KRC)
#
#   Check for exceedence of world records
*/

#pragma once
#include "station.h"
#include "Utilities.h"
#include "utils.h"
#include "python_function.h"

#include <vector>
#include <map>
#include	<string>
#include <iostream>


namespace INTERNAL_CHECKS
{

	/*Get the WMO region from the station id*/

	std::string krc_get_wmo_region(const std::string& wmoid);

	/*'''
		Set the flags to 1 for correct column if data exists

		: param array locs : locations for flags
		: param array flags : the flags
		: param int col : column to use
		'''
	*/
	void krc_set_flags(std::valarray<size_t>& locs, CStation& station, int col);
		
	/*
		Run the known records check for each variable in list

		:param object station: station to process
		:param list var_list: list of variables to process
		:param list flag_col: which columns to use for which variable
		:param file logfile: logfile to store output
	*/
	void krc(CStation& station, std::vector<std::string> variable_list, std::vector<int>flag_col,
		std::ofstream& logfile);
}