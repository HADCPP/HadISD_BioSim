#include "utils.h"
#include "station.h"
#include<vector>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <map>
#include <string>
#include "python_function.h"
#include <valarray> 
#include <boost/numeric/ublas/matrix.hpp>
#include <algorithm>
#include<fstream>
#include <ctime>

using namespace std;
using namespace boost;
using namespace PYTHON_FUNCTION;

namespace UTILS
{
	/*
		Returns locations of month starts(using hours as index)
	*/
	


	void month_starts(boost::gregorian::date start, boost::gregorian::date end,std::vector<int>& month_locs)
	{
		
		boost::gregorian::date Date = start;
		while (Date < end)
		{
			boost::gregorian::date_duration difference = Date - start;
			month_locs.push_back(difference.days() * 24);
				//increment counter
			if (Date.month() < 9)
			{
				/*cout << Date.month() << "";
				cout << to_string(Date.year()) + "0" + to_string((Date.month() + 1)) + "01" << endl;*/
				Date = boost::gregorian::date_from_iso_string(to_string(Date.year()) + "0" + to_string((Date.month() + 1)) + "01");
			}
			else if (Date.month() >= 9 && Date.month() < 12)
			{
				/*cout << Date.month() << "";
				cout << to_string(Date.year()) + to_string((Date.month() + 1)) + "01" << endl;*/
				Date = boost::gregorian::date_from_iso_string(to_string(Date.year()) + to_string((Date.month() + 1)) + "01");
			}
			else
			{
				/*cout << Date.month() << "";
				cout << to_string(Date.year()+1) + "0101" << endl;*/
				Date = boost::gregorian::date_from_iso_string(to_string((Date.year() + 1)) + "0101");
			}
					
		}
		
	}
	/*
		Create array of month start / end pairs
			: param datetime start : start of data
			: param datetime end : end of data
			: returns : month_ranges : Nx2 array
	*/
	inline std::vector<std::pair<int,int>> month_starts_in_pairs(boost::gregorian::date start, boost::gregorian::date end)
	{
		//set up the arrays of month start locations
		std::vector<int> m_starts;
		month_starts(start, end,m_starts);

		std::vector<std::pair<int, int>> month_ranges;
			
		for (int i =0 ; i < m_starts.size()-1; i++)
		{
			month_ranges.push_back(make_pair(m_starts[i],m_starts[i + 1]));
		}
		boost::gregorian::date_duration difference = end - start;
			
		month_ranges.push_back(make_pair(m_starts[m_starts.size()-1],difference.days() * 24));

		return month_ranges;
	}

	valarray<bool> create_fulltimes(CStation& station, std::vector<string>& var_list, boost::gregorian::date start, boost::gregorian::date end, std::vector<string>& opt_var_list, bool do_input_station_id, bool do_qc_flags, bool do_flagged_obs)
	{
		//expand the time axis of the variables
		boost::gregorian::date_duration DaysBetween = end - start;
		valarray < float > fulltimes = PYTHON_FUNCTION::arange<float>(DaysBetween.days() * 24,0);
				
		//adjust if input netCDF file has different start date to desired
		string time_units = station.getTime_units();
		//Extraire la date de la chaîne
		time_units = time_units.substr(time_units.find("since "));  
		time_units = time_units.substr(time_units.find(" "));
		time_units.erase(std::remove(time_units.begin(), time_units.end(), ' '), time_units.end());
		boost::gregorian::date  netcdf_start = boost::gregorian::date_from_iso_string(time_units);
		boost::gregorian::date_duration offset = start- netcdf_start;
		
		fulltimes += offset.days()*24;
		

		valarray<bool> match, match_reverse;

		match = PYTHON_FUNCTION::in1D<float, float>(fulltimes, station.getMetvar("time").getData());
		match_reverse = PYTHON_FUNCTION::in1D<float, float>(station.getMetvar("time").getData(), fulltimes);

		//if optional/carry through variables given, then set to extract these too
		std::vector<string> full_var_list = var_list;
		
		if (opt_var_list.size() != 0)
		{
			copy(opt_var_list.begin(), opt_var_list.end(), std::back_inserter(full_var_list));
		}
	
		//if (do_input_station_id)  full_var_list.push_back("input_station_id");

		for (string variable : full_var_list)
		{
			CMetVar& st_var = station.getMetvar(variable);
			//use masked arrays for ease of filtering later
			
			CMaskedArray<float> news(Cast<float>(st_var.getMdi()), fulltimes.size());
			news.masked(Cast<float>(st_var.getMdi()));
			
			if (match_reverse.size()!= 0)	
			{
				CMaskedArray<float> dummy=st_var.getData()[match_reverse];
				dummy.masked(Cast<float>(st_var.getMdi()));
				news.fill(match,dummy);
			}
			//but re-mask those filled timestamps which have missing data
						
			station.getMetvar(variable).setData(news);
						
			if (find(var_list.begin(), var_list.end(), variable) != var_list.end() && do_flagged_obs == true)
			{
				//flagged values
				news.m_data =Cast<float>(st_var.getMdi());
				if (st_var.getFlagged_obs().size() != 0)
				{
					varrayfloat dummy = st_var.getFlagged_obs().m_data[match_reverse];
					news.fill(match,dummy);
				}
				st_var.setFlagged_obs(news);

				//flags - for filtering
				news.m_data = Cast<float>(st_var.getMdi());
				if (st_var.getFlags().size() != 0)
				{

					varrayfloat dummy = st_var.getFlags().m_data[match_reverse];
					news.fill(match,dummy);
				}
				st_var.setFlags(news);
				
			}	
		}
		//do the QC flags, using try / except
		if (do_qc_flags == true && station.getQc_flags().size()!=0)
		{
			std::valarray<varrayfloat> qc_var = station.getQc_flags();
			std::valarray<varrayfloat> news(fulltimes.size());
			varrayfloat col(float(0),69);
			for (size_t i = 0; i< news.size(); ++i)
			{
				news[i] = col;
			}
			std::valarray<varrayfloat> dummy = qc_var[match];
			news[match] = dummy;
			station.setQc_flags(news);
		}

		//working in fulltimes throughout and filter by missing
		if (offset.days() != 0)	fulltimes = PYTHON_FUNCTION::arange<float>(DaysBetween.days() * 24, 0);

		station.getMetvar("time").setData(fulltimes);
		
		return match;
	}

		

	valarray<std::pair<float, float>> monthly_reporting_statistics(CMetVar& st_var, boost::gregorian::date start, boost::gregorian::date end)
	{
		std::vector<pair<int, int>> monthly_ranges;
		monthly_ranges= month_starts_in_pairs(start, end);
		valarray<std::pair<float, float>> reporting_stats(monthly_ranges.size());

		for (size_t i = 0; i < reporting_stats.size(); ++i)
		{
			reporting_stats[i] = make_pair(float(-1), float(-1));
		}
		int m = 0;
		for (pair<int,int> month:monthly_ranges)
		{
			size_t taille = month.second - month.first + 1;
			varraysize indices(taille);
			for (size_t i = 0; i < taille; i++)
			{
				indices[i] = month.first + i;
			}
			CMaskedArray<float> dummy = st_var.getAllData()[indices];
			reporting_stats[m] = make_pair(reporting_frequency(dummy), reporting_accuracy(dummy));
		}
		return reporting_stats;
	}
	 void print_flagged_obs_number(std::ofstream& logfile, string test, string variable, int nflags)
	{
		logfile << test << "     Check Flags  :     " << variable << "    :    " << nflags << endl;
	}
	/*
	Apply these flags to all variables
		: param object CStation : the CStation object to be processed
		: param list all_variables : the variables where the flags are to be applied
		: param file logfile : logfile to store outputs
		: returns :
		*/
	 void apply_flags_all_variables(CStation& station, const std::vector<string>& full_variable_list, int flag_col, ofstream& logfile, string test)
	{
		varraysize flag_locs = npwhere(station.getQc_flags()[flag_col], "!", float(0));

		for (string var : full_variable_list)
		{
			CMetVar& st_var = station.getMetvar(var);
			//copy flags into attribute
			st_var.setFlags(flag_locs, float(1));
			logfile << "Applying  " << test<< "flags to  "<<var<< endl;
		}
	}
	/*
	Applying windspeed flags to wind directions synergistically
    Called after every test which assess windspeeds

    :param object station: the station object to be processed
	*/
	void apply_windspeed_flags_to_winddir(CStation& station)
	{
		CMetVar & windspeeds = station.getMetvar("windspeeds");
		CMetVar & winddirs = station.getMetvar("winddirs");

		winddirs.setFlags(windspeeds.getFlags());
	}

	 void append_history(CStation& station, string text)
	{
		time_t _time;
		struct tm timeInfo;
		char format[32];
		time(&_time);
		localtime_s(&timeInfo, &_time);
		strftime(format, 32, "%Y-%m-%d %H-%M", &timeInfo);
		station.setHistory(station.getHistory() + text + format);
	}
	/*
		Return the data masked by the flags
		Retourner un array où les valeurs qui correspondent aux indices 
		où st_var.flags==1 sont masquées
	*/

	CMaskedArray<float> apply_filter_flags(CMetVar& st_var)
	{
		return  PYTHON_FUNCTION::ma_masked_where<float,float>(st_var.getFlags().m_data,float(1),st_var.getData(),Cast<float>(st_var.getMdi()));
	}

	inline float reporting_accuracy(CMaskedArray<float>& indata, bool winddir)
	{
		varrayfloat good_values = indata.compressed();
		float resolution = -1;
		if (winddir)
		{
			resolution = 1;
			if (good_values.size() > 0)
			{
				varrayfloat binEdges = PYTHON_FUNCTION::arange<float>(362, 0);
				varrayfloat hist = PYTHON_FUNCTION::histogram(good_values, binEdges);
				//normalise
				hist = hist / hist.sum();

				varrayfloat hist1 = hist[PYTHON_FUNCTION::arange<size_t>(360 + 90, 90, 90)];
				if (hist1.sum() >= 0.6) resolution = 90;
				hist1 = hist[PYTHON_FUNCTION::arange<size_t>(360 + 45, 45, 45)];
				if (hist1.sum() >= 0.6) resolution = 45;
				hist1 = hist[PYTHON_FUNCTION::Arange(360 + 22.5, 22.5, 22.5)];
				if (hist1.sum() >= 0.6) resolution = 10;

				cout << "Wind dir resolution =" << resolution << " degrees" << endl;
			}
		}
		else
		{
			if (good_values.size() > 0)
			{
				varrayfloat remainders = std::abs(good_values) - std::abs(good_values.apply(UTILS::MyApplyRoundFunc));
				varrayfloat binEdges = PYTHON_FUNCTION::arange<float>(1.05, -0.05, 0.1);
				varrayfloat hist = PYTHON_FUNCTION::histogram(good_values, binEdges);

				hist = hist / hist.sum();
				if (hist[0] >= 0.3)
					if (hist[5] >= 0.15)
						resolution = 0.5;
					else
						resolution = 1.0;
				else
					resolution = 0.1;

			}
		}
		return resolution;
	}

	inline float reporting_accuracy(varrayfloat& good_values, bool winddir)
	{
		
		float resolution = -1;
		if (winddir)
		{
			resolution = 1;
			if (good_values.size() > 0)
			{
				varrayfloat binEdges = PYTHON_FUNCTION::arange<float>(362, 0);
				varrayfloat hist = PYTHON_FUNCTION::histogram(good_values, binEdges);
				//normalise
				hist = hist / hist.sum();

				varrayfloat hist1 = hist[PYTHON_FUNCTION::arange<size_t>(360 + 90, 90, 90)];
				if (hist1.sum() >= 0.6) resolution = 90;
				hist1 = hist[PYTHON_FUNCTION::arange<size_t>(360 + 45, 45, 45)];
				if (hist1.sum() >= 0.6) resolution = 45;
				hist1 = hist[PYTHON_FUNCTION::Arange(360 + 22.5, 22.5, 22.5)];
				if (hist1.sum() >= 0.6) resolution = 10;

				cout << "Wind dir resolution =" << resolution << " degrees" << endl;
			}
		}
		else
		{
			if (good_values.size() > 0)
			{
				varrayfloat remainders = std::abs(good_values) - std::abs(good_values.apply(UTILS::MyApplyRoundFunc));
				varrayfloat binEdges = PYTHON_FUNCTION::arange<float>(1.05, -0.05, 0.1);
				varrayfloat hist = PYTHON_FUNCTION::histogram(good_values, binEdges);

				hist = hist / hist.sum();
				if (hist[0] >= 0.3)
					if (hist[5] >= 0.15)
						resolution = 0.5;
					else
						resolution = 1.0;
				else
					resolution = 0.1;

			}
		}
		return resolution;
	}
	float reporting_frequency(CMaskedArray<float>& indata)
	{
		varraysize masked_locs = npwhere(indata.m_mask, "=", false);
		float frequency = float(-1);

		if (masked_locs.size()>0)
		{
			varraysize  difference_series = npDiff(masked_locs);
			varrayfloat binEdges = PYTHON_FUNCTION::arange<float>(25, 1);
			varrayfloat hist = PYTHON_FUNCTION::histogram(difference_series, binEdges,true);
			if (hist[0] >= 0.5) frequency = float(1);
			else if (hist[1] >= 0.5) frequency = float(2);
			else if (hist[2] >= 0.5) frequency = float(3);
			else if (hist[3] >= 0.5) frequency = float(4);
			else if (hist[5] >= 0.5) frequency = float(6);
			else frequency = float(24);
		}
		return frequency;
	}
	
	/* create bins and bin centres from data 
		given bin width covers entire range */
	
	void create_bins(varrayfloat& indata, float binwidth, varrayfloat& bins, varrayfloat& bincenters)
	{
		//set up the bins
		int bmins = int(floor(indata.min()));
		int bmax = int(ceil(indata.max()));
		bins = PYTHON_FUNCTION::arange<float>(bmax + (3. * binwidth), bmins - binwidth, binwidth);
		bincenters.resize(bins.size());
		for (int i = 0; i < bins.size() - 1; i++)
				bincenters[i] = 0.5*(bins[i] + bins[i+1]);

	}
	/* create bins and bin centres from data given bin width covers entire range */

	void create_bins(std::valarray<CMaskedArray<float>>& indata, float binwidth, varrayfloat& bins, varrayfloat& bincenters)
	{
		
		//Search min and min of the matrice
		varrayfloat Min(indata.size());
		varrayfloat Max(indata.size());

		for (size_t i = 0; i < indata.size(); i++)
		{
			Min[i] = indata[i].compressed().min();
			Max[i] = indata[i].compressed().max();
		}

		int bmins = int(floor(Min.min()));
		int bmax = int(ceil(Max.max()));
		//set up the bins
		bins = PYTHON_FUNCTION::arange<float>(bmax + (3. * binwidth), bmins - binwidth, binwidth);
		bincenters.resize(bins.size()-1);
		for (int i = 0; i < bins.size() - 1; i++)
			bincenters[i] = 0.5*(bins[i] + bins[i + 1]);

	}
	//Sum up a single month across all years(e.g.all Januaries)
	//return this_month, year_ids, datacount

	std::valarray<int>  concatenate_months(std::valarray<std::pair<int, int>>& month_ranges, CMaskedArray<float>& data, std::vector<CMaskedArray<float>>& this_month,
		std::vector<int>& year_ids, float missing_value, bool hours)
	{

		valarray<int> datacount(month_ranges.size());
		for (size_t y = 0; y < month_ranges.size();y++)
		{

			
			CMaskedArray<float> this_year = data(month_ranges[y].first, month_ranges[y].second);
			this_year.masked(data.m_fill_value);
			datacount[y] = this_year.compressed().size();

			if (y == 0)
			{
				//store so can access each hour of day separately
				if (hours)
					this_month = C_reshape(this_year, 24);
				else
					this_month.push_back(this_year);
				for (int i = 0; i < this_month.size(); ++i)
					year_ids.push_back(y);
			}
			else
			{
				if (hours)
				{
					std::vector<CMaskedArray<float>> t_years = C_reshape(this_year, 24);
					std::copy(t_years.begin(), t_years.end(), std::back_inserter(this_month));
				}
				else
					this_month.push_back(this_year);
				for (int i = 0; i < int(this_year.size()/24); ++i)
					year_ids.push_back(y);
			}
			
		}
		return datacount;
	}

	//''' Calculate the percentile of data '''
	inline float percentiles(varrayfloat& data, float percent, bool idl)
	{
		varrayfloat sorted_data(data);
		std::sort(std::begin(sorted_data), std::end(sorted_data));
		int n_data;
		float percentile;
		if (idl)
		{
			n_data = data.size() - 1;
			percentile = sorted_data[int(ceil(n_data * percent))]; // matches IDL formulation
		}
		else
		{
			n_data = data.size();
			percentile = sorted_data[int(n_data * percent)];
		}

		return percentile;
	}
		
	void winsorize(varrayfloat& data, float percent, bool idl )
	{
		for (float pct : { percent, 1 - percent })
		{
			float percentile;
			varraysize locs(data.size());
			if (pct < 0.5)
			{
				percentile = percentiles(data, pct, idl );
				locs = npwhere(data, "<", percentile);
			}
			else
			{
				percentile = percentiles(data, pct, idl);
				locs = npwhere(data, ">", percentile);
			}
			data[locs] = percentile;

		}
	}
	/*Calculate the IQR of the data*/
	float IQR(varrayfloat data, double percentile )
	{
		varrayfloat sorted_data(data);
		std::sort(std::begin(sorted_data), std::end(sorted_data));
		int n_data = sorted_data.size();
		int quartile = int(std::round(percentile*n_data));
		return sorted_data[n_data - quartile] - sorted_data[quartile];
	}

	void reshapeMonth(std::vector<std::valarray<std::pair<int, int>>>& month_ranges_years,std::map<int, int>&  month_ranges)
	{
		int taille = int(month_ranges.size() / 12);
		std::valarray<pair<int, int>> month(taille);
		int index = 0;
		int compteur = 0;
		int iteration = 1;
		map<int, int>::iterator month_it = month_ranges.begin();
		for (int i = 0; i < month_ranges.size(); i += 12)
		{
			if (iteration <= taille && i < month_ranges.size() && compteur < 12)
			{
				month[index++] = make_pair(month_it->first, month_it->second);
				iteration++;
				if (i + 12 < month_ranges.size()) std::advance(month_it, 12);
			}
			else
			{
				if (compteur == 12) break;
				month_ranges_years.push_back(month);
				compteur++;
				month.resize(taille);
				index = 0;
				i = month_ranges_years.size();
				month_it = month_ranges.begin();
				std::advance(month_it, i);
				month[index++] = make_pair(month_it->first, month_it->second);
				std::advance(month_it, 12);
				iteration = 2;
			}
			if (i + 12 >= month_ranges.size() && compteur < 12)
			{
				i = month_ranges_years.size();
				month_it = month_ranges.begin();
				std::advance(month_it, i + 1);
			}
		}
	}
	void reshapeYear(std::vector<std::vector<std::pair<int, int>>>& month_ranges_years, std::map<int, int>&  month_ranges)
	{
		int iteration = 1;
		std::vector<std::pair<int, int>> month;
		
		for (map<int, int>::iterator month_it = month_ranges.begin(); month_it != month_ranges.end(); month_it++)
		{
			if (iteration <= 12)
			{
				month.push_back(make_pair(month_it->first, month_it->second));
				iteration++;
			}
			else
			{
				month_ranges_years.push_back(month);
				
				month.clear();
				month.push_back(make_pair(month_it->first, month_it->second));
				iteration = 2;
			}

		}
		month_ranges_years.push_back(month);
	}

	float mean_absolute_deviation(const varrayfloat& data, bool median )
	{
		//''' Calculate the MAD of the data '''
		float mad=0;
		if (median)
		{
			float med = idl_median(data);
			for (size_t i = 0; i < data.size();i++)
				mad += abs(data[i]-med);

			mad /= data.size();
		}

		else
		{

			float mean = data.sum()/data.size();
			for (size_t i = 0; i < data.size(); i++)
				mad += abs(data[i] - mean);

			mad /= data.size();
		}
		
		return mad;// mean_absolute_deviation
	}

}
