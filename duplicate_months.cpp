
#include "station.h"
#include <vector>
#include <string>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <map>
#include "utils.h"
#include <iostream>
#include "python_function.h"
using namespace std;
using namespace UTILS;
using namespace PYTHON_FUNCTION;


namespace INTERNAL_CHECKS
{
	/*
		Perform test whether month is duplicated
		: param array source_data : source data
		: param array target_data : target data
		: param array valid : valid values
		: param int sm : source month counter
		: param int tm : target month counter
		: param array duplicated : True / False array for duplicated months
	*/
	void is_month_duplicated(const valarray<float>& source_data, const valarray<float>& target_data, const varraysize& valid, int sm,
		int tm, valarray<int>& duplicated)
	{
		varrayfloat dum1 = source_data[valid];
		varrayfloat dum2 = target_data[valid];
		varraysize match = PYTHON_FUNCTION::npwhere(dum1,"=",dum2);
	
		if (match.size() == valid.size())
		{
			duplicated[sm] = 1;
			duplicated[sm + 1 + tm] = 1;
		}
	}

	/*
		Pass into test whether month is duplicated and set flags
		: param array source_data : source data
		: param array target_data : target data
		: param array valid : valid values
		: param int sm : source month counter
		: param int tm : target month counter
		: param array duplicated : True / False array for duplicated months
	*/

	void duplication_test(const valarray<float>& source_data, const valarray<float>& target_data, const varraysize& valid, int sm, int tm,
		std::vector<pair<int, int>> ::const_iterator source_month, std::vector<pair<int, int>> ::const_iterator  target_month, valarray<int>& duplicated, CStation& station, int flag_col)
	{
		
		is_month_duplicated(source_data, target_data, valid, sm, tm, duplicated);
		size_t taille_slice = source_month->second - source_month->first;
		varraysize indices0(taille_slice);
		
		for (size_t i = 0; i < taille_slice; i++)
		{
			indices0[i] = source_month->first + i;

		}
		
		taille_slice = target_month->second - target_month->first;
		varraysize indices1(taille_slice);
		for(size_t i = 0; i < taille_slice; i++)
		{
			indices1[i] = target_month->first + i;
		}

		if (duplicated[sm] == 1) //make flags if this month is duplicated somewhere
		{
			station.setQc_flags(float(1), indices0, flag_col);
			station.setQc_flags(float(1), indices1, flag_col);
		}
	}
	/*
	:param obj CStation stat: CStation object with suitable attributes
	:param list variable_list: list of netcdf variables to process
	:param list full_variable_list: the variables for flags to be applied to
	:param list flag_col: which column to set in flag array
	:param datetime start: data start
	:param datetime end: data end
	:param file logfile: logfile to store outputs

	*/
	void dmc(CStation& station, std::vector<std::string>& variable_list, std::vector<std::string>& full_variable_list,
		int flag_col, boost::gregorian::date start, boost::gregorian::date end, ofstream &logfile)
	{
		const  int MIN_DATA_REQUIRED = 20;
		//get array of MAP<INT,INT> start / end pairs
		std::vector<pair<int, int>> month_ranges;
		month_ranges = UTILS::month_starts_in_pairs(start, end);
		
		int v = 0;
		for (string variable: variable_list)
		{
			CMetVar & st_var = station.getMetvar(variable); //Recuperer la variable meteo de la CStation
			float missing_value = Cast<float>(st_var.getMdi());
			valarray<int> duplicated(0,month_ranges.size());
			
			int sm = 0;
			for (std::vector<pair<int, int>> ::const_iterator source_month = month_ranges.begin(); source_month != month_ranges.end(); ++source_month)
			{
				cout << "Month  " << (sm + 1) << "  of   " << month_ranges.size() << endl;
				int taille_slice = source_month->second - source_month->first;
				varraysize indices(taille_slice);
				for (size_t i = 0; i < taille_slice; i++)
				{
					indices[i] = source_month->first + i;
				}
				CMaskedArray<float> source_data=st_var.getData()[indices];
				source_data.masked(Cast<float>(st_var.getMdi()));
				if (duplicated[sm] == 0)   // don't repeat if already a duplicated
				{
					std::vector<pair<int, int>> ::const_iterator target_month = source_month;
					target_month++;

					int tm = 0;
					for (; target_month != month_ranges.end(); ++target_month)
					{
						
						int taille_slice = target_month->second - target_month->first ;
						varraysize indices(taille_slice);
						for (size_t i = 0; i < taille_slice; i++)
						{
							indices[i] = target_month->first + i;
						}
						CMaskedArray<float> target_data = st_var.getData()[indices];
						target_data.masked(Cast<float>(st_var.getMdi()));
						//match the data periods
						size_t overlap = std::min(source_data.size(), target_data.size());
						CMaskedArray<float> s_data = source_data.m_data[slice(0, overlap, 1)];  
						CMaskedArray<float> t_data = target_data.m_data[slice(0, overlap, 1)]; 
						s_data.masked( missing_value );
						t_data.masked(missing_value);
												
						varraysize s_valid = npwhere<float>(s_data.compressed(), "!", st_var.getFdi());
						varraysize t_valid = npwhere<float>(t_data.compressed(), "!", st_var.getFdi());
						//if enough of an overlap
						if (s_valid.size() >= MIN_DATA_REQUIRED && t_valid.size() >= MIN_DATA_REQUIRED)
						{
							if (s_valid.size() < t_valid.size())
								duplication_test(source_data.m_data, target_data.m_data, s_valid, sm, tm, source_month, target_month, duplicated, station, flag_col);
							else //swap the list of valid points 
								duplication_test(source_data.m_data, target_data.m_data, t_valid, sm, tm, source_month, target_month, duplicated, station, flag_col);
						}
						tm++;
					}//target month
				}
				sm++;
			}//souce month
			varraysize flag_locs = npwhere(station.getQc_flags()[flag_col], "!", float(0));
			UTILS::print_flagged_obs_number(logfile, "Duplicate Month", variable, flag_locs.size());
			//copy flags into attribute
			st_var.setFlags(flag_locs, float(1));
		}

		UTILS::apply_flags_all_variables(station, variable_list, 0, logfile, "Duplicate Months check");
		append_history(station, "Duplicate Months Check");

	}//End Duplicate Check

}