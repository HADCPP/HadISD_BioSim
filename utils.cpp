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
using namespace boost::gregorian;
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

	valarray<bool> create_fulltimes(CStation& station,  std::vector<string> var_list, boost::gregorian::date start, boost::gregorian::date end, std::vector<string> opt_var_list, bool do_input_station_id, bool do_qc_flags, bool do_flagged_obs)
	{
		//expand the time axis of the variables
		boost::gregorian::date_duration DaysBetween = end - start;
		valarray <float> fulltimes = PYTHON_FUNCTION::arange<float>(DaysBetween.days() * 24,0);
				
		//adjust if input  file has different start date to desired
		string time_units = station.getMetvar("time").getUnits();
		//Extraire la date de la chaîne
		/*time_units = time_units.substr(time_units.find("since"));  
		time_units = time_units.substr(time_units.find(" "));
		time_units.erase(std::remove(time_units.begin(), time_units.end(), ' '), time_units.end());
		boost::gregorian::date  netcdf_start = boost::gregorian::date_from_iso_string(time_units);
		boost::gregorian::date_duration offset = start- netcdf_start;*/
		
		//fulltimes += offset.days()*24;
		

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
			
			CMaskedArray<float> news(st_var.getMdi(), fulltimes.size());
			news.masked(st_var.getMdi());
			
			if (match_reverse.size()!= 0)	
			{
				CMaskedArray<float> dummy=st_var.getData()[match_reverse];
				dummy.masked(st_var.getMdi());
				news.fill(match,dummy);
			}
			//but re-mask those filled timestamps which have missing data
						
			station.getMetvar(variable).setData(news);
						
			if (find(var_list.begin(), var_list.end(), variable) != var_list.end() && do_flagged_obs == true)
			{
				//flagged values
				news.m_data =st_var.getMdi();
				if (st_var.getFlagged_obs().size() != 0)
				{
					varrayfloat dummy = st_var.getFlagged_obs().m_data[match_reverse];
					news.fill(match,dummy);
				}
				st_var.setFlagged_obs(news);

				//flags - for filtering
				news.m_data = st_var.getMdi();
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
		fulltimes = PYTHON_FUNCTION::arange<float>(DaysBetween.days() * 24, 0);

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
			logfile << "  Applying  " << test<< "   flags to     "<<var<< endl;
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
	
	void create_bins(const varrayfloat& indata, float binwidth, varrayfloat& bins, varrayfloat& bincenters)
	{
		//set up the bins
		int bmins = int(floor(indata.min()));
		float max = indata.max();
		int bmax = int(ceil(indata.max()));
		bins = PYTHON_FUNCTION::arange<float>(bmax + (3. * binwidth), bmins - binwidth, binwidth);
		bincenters.resize(bins.size());
		for (int i = 0; i < bins.size() - 1; i++)
				bincenters[i] = 0.5*(bins[i] + bins[i+1]);
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

	//Get the distance between two points long Earth's surface
	pair<int, int>  get_dist_and_bearing(pair<float, float> coord1, pair<float, float> coord2)
	{
		float lat1 = coord1.first;
		float lon1 = coord1.second;
		float lat2 = coord2.first;
		float lon2 = coord2.second; 

		double R = 6371.229; // km

		float phi1 = get_phi(lat1);
		float phi2 = get_phi(lat2);

		float theta1 = deg2rad(lon1);
		float theta2 = deg2rad(lon2);

		float cosinus = (std::sin(phi1) * std::sin(phi2) * std::cos(theta1 - theta2) + std::cos(phi1) * std::cos(phi2));
		float arc = std::acos(cosinus);

		lat1 = deg2rad(lat1);
		lon1 = deg2rad(lon1);
		lat2 = deg2rad(lat2);
		lon2 = deg2rad(lon2);
			
		float bearing = rad2deg(std::atan2(std::sin(lon2 - lon1)*std::cos(lat2), std::cos(lat1)*std::sin(lat2) - std::sin(lat1)*std::cos(lat2)*std::cos(lon2 - lon1)));

		if (bearing < 0)  bearing += 360.;

		float distance = arc*R;

		return make_pair(int(distance), int(bearing));
			
	}

}

namespace DATA_READING
{
	void read_station_metadata(const std::string file, std::vector<CStation>& station_info)
	{


		char_separator<char> sep(",");
		std::ifstream input;
		std::stringstream sst;
		sst << file;
		std::string chaine = sst.str();
		try
		{
			input.open(chaine.c_str());
		}
		catch (std::exception& e)
		{
			cout << e.what() << endl;
		}
		if (!input) exit(1);
		std::string ligne = "";
		getline(input, ligne);
		/*while (!input.eof())
		{
		getline(input, ligne);
		int i = 0;
		std::string data[22] = { "0" };
		tokenizer<char_separator<char>> tokens(ligne, sep);
		for (tokenizer<char_separator<char>>::const_iterator t = tokens.begin(); t != tokens.end();t++)
		{

		(t->c_str() != "") ? data[i] = t->c_str() : data[i] = "";
		i++;
		}

		CStation station = CStation::CStation(data[0], data[1], atof(data[2].c_str()), atof(data[3].c_str()),
		atof(data[4].c_str()), data[12].c_str());
		station_info.push_back(station);
		}*/


		char  delimiter = ',';
		while (!input.eof())
		{
			getline(input, ligne);
			std::stringstream iss(ligne);
			std::string token;
			int i = 0;
			std::string data[22] = { "0" };

			while (getline(iss, token, delimiter))
			{
				data[i] = token.c_str();
				i++;
			}
			CStation stat = CStation(data[0], data[1], atof(data[2].c_str()), atof(data[3].c_str()), atof(data[4].c_str()), data[12].c_str());
			station_info.push_back(stat);

		}
	}
	// look if the csv file exist for the station
	std::string find_CSV_file(const CStation& station)
	{

		std::string nom_file = (station).getName() + " [" + (station).getId() + "]";
		std::string file_name = nom_file + ".csv";
		std::string fichier = CSV_OUTFILE_LOCS + file_name;
		
		boost::filesystem::path p{ fichier };

		if (!boost::filesystem::exists(p))
		{
			
			std::string nom = (station).getName();
			auto pos = std::string::npos;
			while ((pos = nom.find('é')) != std::string::npos || (pos = nom.find('è')) != std::string::npos || (pos = nom.find('ê')) != std::string::npos)
			{
				nom.replace(pos, 1, "e");
			}
			while ((pos = nom.find('à')) != std::string::npos)
			{
				nom.replace(pos, 1, "a");
			}
			while ((pos = nom.find('\'')) != std::string::npos)
			{
				nom.replace(pos, 1, " ");
			}
			file_name = nom + " [" + (station).getId() + "]" + ".csv";
			fichier = CSV_OUTFILE_LOCS + file_name;
			boost::filesystem::path p1{ fichier };
			if (boost::filesystem::exists(p1))
			{
				
				cout << p1.filename() << endl;
				cout << "\n Reading data from CSV files \n" << endl;
				cout << "Reading data from  " << p1.parent_path() << endl;

				return file_name;
			}
			else exit(1);
		}

		else
		{
			cout << p.filename() << endl;
			cout << "\n Reading data from CSV files \n" << endl;
			cout << "Reading data from  " << p.parent_path() << endl;
			return file_name;
		}
	}
	//Read each station variables data and put them in the appropriate data structures.
	bool readData(CStation& station, boost::gregorian::date  DATESTART, boost::gregorian::date DATEEND, test& internal_tests)
	{

		string fichier = find_CSV_file(station);
		fichier = CSV_OUTFILE_LOCS + fichier;
		cout << "Station Identifier:  " << (station).getId() << endl;

		std::string namefile = station.getId();

		//duree entre DATASTART et DATAEND
		date end = DATEEND;// + years(1);
		boost::gregorian::date_duration DaysBetween = end - DATESTART;
		int HoursBetween = int(DaysBetween.days()) * 24;
		std::vector<int> TimeStamps;
		PYTHON_FUNCTION::linspace<int>(TimeStamps, 0, HoursBetween - 1, HoursBetween);
		std::vector<int> ValidYears;
		PYTHON_FUNCTION::linspace<int>(ValidYears, int(DATESTART.year()), int(DATEEND.year())-1, int(DATEEND.year() - DATESTART.year()));
		std::string dubiousfile = LOG_OUTFILE_LOCS + "dubious_data_files.txt";
		boost::gregorian::date  dbg_sttime = boost::gregorian::day_clock::local_day();

		//Création et initialisation des tableaux où mettre les variables meteo

		std::valarray<float> temperatures(INTMDI, HoursBetween);
		std::valarray<int> temperature_flags(INTMDI, HoursBetween);
		std::valarray<float> dewpoints(INTMDI, HoursBetween);
		std::valarray<int> dewpoint_flags(INTMDI, HoursBetween);
		std::valarray<float> windspeeds(INTMDI, HoursBetween);
		std::valarray<int> windspeeds_flags(INTMDI, HoursBetween);
		std::valarray<float> humidity(INTMDI, HoursBetween);
		std::valarray<int> humiditys_flags(INTMDI, HoursBetween);
		std::valarray<float> winddirs(INTMDI, HoursBetween);
		std::valarray<int> winddirs_flags(INTMDI, HoursBetween);
		std::valarray<float>  slp(INTMDI, HoursBetween);
		std::valarray<int> slp_flags(INTMDI, HoursBetween);
		std::valarray<float> precip(INTMDI, HoursBetween);
		std::valarray<int> precip_flags(INTMDI, HoursBetween);

		//if extra : Ajouter des variables meteo supplémentaires

		boost::gregorian::date  dbg_lasttime = boost::gregorian::day_clock::local_day();

		//ouvrir le fichier csv qui contient des données
		int  last_obs_time = 0;
		std::ifstream input;
		std::stringstream sst;
		sst << fichier;
		std::string chaine = sst.str();

		boost::filesystem::exists(chaine.c_str());
		try
		{
			input.open(chaine.c_str());
		}
		catch (std::exception e)
		{
			std::cout << e.what() << endl;
		}
		if (!input) return true;
		std::string ligne = "";
		std::string token;
		char  delimiter = ',';
		map<std::string, int> Headings;
		// Lecture de l'entête du fichier et récupération des entêtes des principales variables (si elles existent)
		// Indices des variables;
		int temp, dew, hum, windS, windD, pres, precipitation;
		getline(input, ligne);
		char_separator<char> sep(",");
		tokenizer<char_separator<char>> tokens(ligne, sep);
		int i = 0;
		for (tokenizer<char_separator<char>>::const_iterator t = tokens.begin(); t != tokens.end(); t++)
		{
			Headings[t->c_str()] = i++;
		}
		bool sortie = true;

		if (Headings.find("Tair") != Headings.end())
		{
			temp = Headings["Tair"];
			internal_tests.climatological = true;
			internal_tests.diurnal = true;
			internal_tests.duplicate = true;
			internal_tests.frequent = true;
			internal_tests.gap = true;
			internal_tests.humidity = true;
			internal_tests.odd = true;
			internal_tests.records = true;
			internal_tests.spike = true;
			internal_tests.streaks = true;
			internal_tests.variance = true;

		}
		if (Headings.find("Tdew") != Headings.end())
		{
			dew = Headings["Tdew"];

			internal_tests.climatological = true;
			internal_tests.frequent = true;
			internal_tests.gap = true;
			internal_tests.humidity = true;
			internal_tests.odd = true;
			internal_tests.records = true;
			internal_tests.spike = true;
			internal_tests.streaks = true;
			internal_tests.variance = true;

		}
		if (Headings.find("Pres") != Headings.end())
		{
			pres = Headings["Pres"];


			internal_tests.frequent = true;
			internal_tests.gap = true;
			internal_tests.odd = true;
			internal_tests.records = true;
			internal_tests.spike = true;
			internal_tests.streaks = true;
			internal_tests.variance = true;

		}
		if (Headings.find("WndS") != Headings.end())
		{
			windS = Headings["WndS"];


			internal_tests.frequent = true;
			internal_tests.odd = true;
			internal_tests.records = true;
			internal_tests.spike = true;
			internal_tests.streaks = true;
			internal_tests.variance = true;
			internal_tests.winds = true;
		}
		if (Headings.find("WndD") != Headings.end())
		{
			windD = Headings["WndD"];

			internal_tests.streaks = true;
			internal_tests.winds = true;
		}
		if (Headings.find("Prcp") != Headings.end())  precipitation = Headings["Prcp"];
		if (Headings.find("RelH") != Headings.end())  hum = Headings["RelH"];
		int compt = 0; //compteur de la boucle while
		int year_since; // Date de début des données dans le fichier csv
		while (!input.eof())
		{
			getline(input, ligne);
			std::stringstream iss(ligne);
			std::string year, month, day;
			int hour;
			if (ligne == "") break;

			vector<std::string> data;

			while (getline(iss, token, delimiter))
			{
				data.push_back(token.c_str());
			}
			year = data[0];

			// tester si la date de l'observation appartient à l'intervalle d'années désiré
			int annee = std::stoi(year);
			if (compt == 0) year_since = annee;
			if (!(std::find(ValidYears.begin(), ValidYears.end(), annee) != ValidYears.end())) continue;

			month = (data[1].size() == 1) ? "0" + data[1] : data[1];
			day = (data[2].size() == 1) ? "0" + data[2] : data[2];
			hour = atoi(data[3].c_str());
			boost::gregorian::date  dt_time;    // Date au format  aaaammdd
			try
			{
				dt_time = boost::gregorian::date_from_iso_string(year + month + day);
			}
			catch (boost::gregorian::bad_day_of_month)
			{
				std::cout << "error with boost::gregorian";
			}
			//duree entre DATESTART et la date de l'observation

			boost::gregorian::date_duration obs_date = dt_time - DATESTART;

			int  obs_time = obs_date.days() * 24 + hour;

			int time_loc = obs_time;
			//test if this time_loc out of bounds:

			if (time_loc != HoursBetween)
			{
				(Headings.find("Tair") != Headings.end()) ? temperatures[time_loc] = atof(data[temp].c_str()) : temperatures[time_loc] = -999;
				temperature_flags[time_loc] = INTMDI;
				(Headings.find("Tdew") != Headings.end()) ? dewpoints[time_loc] = atof(data[dew].c_str()) : dewpoints[time_loc] = -999;
				dewpoint_flags[time_loc] = INTMDI;
				(Headings.find("RelH") != Headings.end()) ? humidity[time_loc] = atof(data[hum].c_str()) : humidity[time_loc] = -999;
				humiditys_flags[time_loc] = INTMDI;
				(Headings.find("WndS") != Headings.end()) ? windspeeds[time_loc] = atof(data[windS].c_str()) : windspeeds[time_loc] = -999;
				windspeeds_flags[time_loc] = INTMDI;
				(Headings.find("WndD") != Headings.end()) ? winddirs[time_loc] = atof(data[windD].c_str()) : winddirs[time_loc] = -999;
				winddirs_flags[time_loc] = INTMDI;
				(Headings.find("Pres") != Headings.end()) ? slp[time_loc] = atof(data[pres].c_str()) : slp[time_loc] = -999;
				slp_flags[time_loc] = INTMDI;
				(Headings.find("Prcp") != Headings.end()) ? precip[time_loc] = atof(data[precipitation].c_str()) : precip[time_loc] = -999;
				precip_flags[time_loc] = INTMDI;
			}
		}


		vector<std::string> full_var_list = process_var;
		if (carry_thru_vars.size() != 0)
		{
			copy(carry_thru_vars.begin(), carry_thru_vars.end(), std::back_inserter(full_var_list));
		}
		full_var_list.push_back("time");

		for (std::string variable : full_var_list)
		{
			varrayfloat dummy(HoursBetween);
			// temperature  min et max
			if (variable == "temperatures")
			{
				dummy = temperatures[temperatures != float(INTMDI)];
				float max_temp;
				(Headings.find("Tair") != Headings.end() && dummy.size()>0) ? max_temp = dummy.max() : max_temp = -999;
				float min_temp;
				(Headings.find("Tair") != Headings.end() && dummy.size()>0) ? min_temp = dummy.min() : min_temp = -999;

				CMetVar this_var(variable, "Dry bulb air temperature at screen height (~2m)");
				this_var.setMdi(INTMDI);
				this_var.setUnits("degree_Celsius");
				this_var.setValidMin(min_temp);
				this_var.setValidMax(max_temp);
				this_var.setCoordinates("latitude longitude elevation");
				this_var.setStandard_name("surface_temperature");
				this_var.setFdi(-2.e30);

				this_var.setData(temperatures);
				varrayfloat flagged_obs(INTMDI, HoursBetween);
				this_var.setFlagged_obs(flagged_obs);
				flagged_obs.resize(HoursBetween,0);
				this_var.setFlags(flagged_obs);
				//Ajouter la variable méteo à la liste des variables de la station stat
				station.setMetVar(this_var, variable);
			}
			if (variable == "dewpoints")
			{
				dummy = dewpoints[dewpoints != float(INTMDI)];
				float max_dew;
				(Headings.find("Tdew") != Headings.end() && dummy.size()>0) ? max_dew = dummy.max() : max_dew = -999;
				float min_dew;
				(Headings.find("Tdew") != Headings.end() && dummy.size()>0) ? min_dew = dummy.min() : min_dew = -999;

				CMetVar this_var(variable, "Dew point temperature at screen height (~2m)");
				this_var.setMdi(INTMDI);
				this_var.setUnits("degree_Celsius");
				this_var.setCellmethods("latitude: longitude: time: point (nearest to reporting hour)");
				this_var.setValidMin(min_dew);
				this_var.setValidMax(max_dew);
				this_var.setCoordinates("latitude longitude elevation");
				this_var.setStandard_name("dew_point_temperature");
				this_var.setFdi(-2.e30);

				this_var.setData(dewpoints);
				varrayfloat flagged_obs(INTMDI, HoursBetween);
				this_var.setFlagged_obs(flagged_obs);
				flagged_obs.resize(HoursBetween, 0);
				this_var.setFlags(flagged_obs);
				//Ajouter la variable méteo à la liste des variables de la station stat
				station.setMetVar(this_var, variable);
			}
			if (variable == "windspeeds")
			{
				dummy = windspeeds[windspeeds != float(INTMDI)];
				float max_ws;
				(Headings.find("WndS") != Headings.end() && dummy.size()>0) ? max_ws = dummy.max() : max_ws = -999;
				float min_ws;
				(Headings.find("WndS") != Headings.end() && dummy.size()>0) ? min_ws = dummy.min() : min_ws = -999;

				CMetVar this_var(variable, "Wind speed at mast height (~10m)");
				this_var.setMdi(INTMDI);
				this_var.setUnits("meters per second");
				this_var.setCellmethods("latitude: longitude: time: point (nearest to reporting hour)");
				this_var.setValidMin(min_ws);
				this_var.setValidMax(max_ws);
				this_var.setCoordinates("latitude longitude elevation");
				this_var.setStandard_name("wind_speed");
				this_var.setFdi(-2.e30);

				this_var.setData(windspeeds);
				varrayfloat flagged_obs(INTMDI, HoursBetween);
				this_var.setFlagged_obs(flagged_obs);
				flagged_obs.resize(HoursBetween, 0);
				this_var.setFlags(flagged_obs);
				//Ajouter la variable méteo à la liste des variables de la station stat
				station.setMetVar(this_var, variable);
			}
			if (variable == "winddirs")
			{

				CMetVar this_var(variable, "Wind Direction at mast height (~10m)");
				this_var.setMdi(INTMDI);
				this_var.setUnits("degree");
				this_var.setCellmethods("latitude: longitude: time: point (nearest to reporting hour)");
				this_var.setValidMin(0);
				this_var.setValidMax(360);
				this_var.setCoordinates("latitude longitude elevation");
				this_var.setStandard_name("wind_from_direction");
				this_var.setFdi(-2.e30);

				this_var.setData(winddirs);
				varrayfloat flagged_obs(INTMDI, HoursBetween);
				this_var.setFlagged_obs(flagged_obs);
				flagged_obs.resize(HoursBetween, 0);
				this_var.setFlags(flagged_obs);
				//Ajouter la variable méteo à la liste des variables de la station stat
				station.setMetVar(this_var, variable);
			}

			if (variable == "slp")
			{
				dummy = slp[slp != float(INTMDI)];
				float smax;
				float smin;
				(Headings.find("pres") != Headings.end() && dummy.size()>0) ? smax = dummy.max() : smax = -999;
				(Headings.find("pres") != Headings.end() && dummy.size()>0) ? smin = dummy.min() : smin = -999;

				CMetVar this_var(variable, "Reported Sea Level Pressure at screen height (~2m)");
				this_var.setMdi(INTMDI);
				this_var.setUnits("hPa");
				this_var.setCellmethods("latitude: longitude: time: point (nearest to reporting hour)");
				this_var.setValidMin(smin);
				this_var.setValidMax(smax);
				this_var.setCoordinates("latitude longitude elevation");
				this_var.setStandard_name("air_pressure_at_sea_level");
				this_var.setFdi(-2.e30);

				this_var.setData(slp);
				varrayfloat flagged_obs(INTMDI, HoursBetween);
				this_var.setFlagged_obs(flagged_obs);
				flagged_obs.resize(HoursBetween, 0);
				this_var.setFlags(flagged_obs);
				//Ajouter la variable méteo à la liste des variables de la station stat
				station.setMetVar(this_var, variable);
			}
			if (variable == "time")
			{
				CMetVar this_var(variable, "time_of_measurement");
				this_var.setUnits("hours since " + to_string(year_since));
				this_var.setCalendar("gregorian | start_year : " + to_string(DATESTART.year()) + " end_month : " + to_string(DATEEND.year()));
				this_var.setValidMin(0);
				this_var.setStandard_name("time");

				dummy.resize(HoursBetween);
				std::iota(std::begin(dummy), std::end(dummy), 0);

				this_var.setData(dummy);
				
				//Ajouter la variable méteo à la liste des variables de la station stat
				station.setMetVar(this_var, variable);

			}
		}

		std::cout << " Reading data to netcdf file" << endl;

		std::cout << "Done station " << station.getName() << endl;
		std::cout << boost::gregorian::day_clock::local_day() - dbg_sttime << " s" << endl;
		std::cout << boost::gregorian::day_clock::local_day() << endl;

		return sortie;
	}

	void writeCSV(CStation& station, boost::gregorian::date  DATESTART, boost::gregorian::date DATEEND)
	{

		//duree entre DATASTART et DATAEND
		date end = DATEEND + years(1);
		boost::gregorian::date_duration DaysBetween = end - DATESTART;
		int HoursBetween = int(DaysBetween.days()) * 24;
		std::vector<int> TimeStamps;
		PYTHON_FUNCTION::linspace<int>(TimeStamps, 0, HoursBetween - 1, HoursBetween);
		std::vector<int> ValidYears;
		PYTHON_FUNCTION::linspace<int>(ValidYears, int(DATESTART.year()), int(DATEEND.year()), int(DATEEND.year() - DATESTART.year() + 1));

		//ouvrir le fichier csv initial en mode lecture qui contient des données
		string csv_file = DATA_READING::find_CSV_file(station);
		// new file for output
		string new_csv_file = NEW_CSV_OUTFILE_LOCS + csv_file;

		csv_file = CSV_OUTFILE_LOCS + csv_file;

		std::ifstream input;
		std::stringstream sst;
		sst << csv_file;
		std::string chaine = sst.str();
		boost::filesystem::exists(chaine.c_str());
		try
		{
			input.open(chaine.c_str());
		}
		catch (std::exception e)
		{
			std::cout << e.what() << endl;
		}
		// ouvrir le fichier de sortie en mode écriture

		std::ofstream output;
		output.open(new_csv_file, std::ofstream::out);

		std::string ligne = "";
		std::string token;
		char  delimiter = ',';
		map<std::string, int> Headings;

		// Lecture de l'entête du fichier et récupération des entêtes des principales variables (si elles existent)
		// Indices des variables;

		int temp, dew, hum, windS, windD, pres, precipitation;
		getline(input, ligne);
		char_separator<char> sep(",");
		tokenizer<char_separator<char>> tokens(ligne, sep);
		int i = 0;
		for (tokenizer<char_separator<char>>::const_iterator t = tokens.begin(); t != tokens.end(); t++)
		{
			Headings[t->c_str()] = i++;
		}

		if (Headings.find("Tair") != Headings.end())  temp = Headings["Tair"];
		if (Headings.find("Tdew") != Headings.end())  dew = Headings["Tdew"];
		if (Headings.find("Pres") != Headings.end())  pres = Headings["Pres"];
		if (Headings.find("WndS") != Headings.end())  windS = Headings["WndS"];
		if (Headings.find("WndD") != Headings.end())  windD = Headings["WndD"];
		if (Headings.find("Prcp") != Headings.end())  precipitation = Headings["Prcp"];
		if (Headings.find("RelH") != Headings.end())  hum = Headings["RelH"];
		
		// Ecriture de l'entête dans Output, ajout d'une colonne pour le flag

		ligne = ligne + ",Flag";
		output << ligne << endl;;

		while (!input.eof())  // parcours du fichier à copier
		{
			getline(input, ligne);
			std::stringstream iss(ligne);
			std::string year, month, day;
			int hour;
			if (ligne == "") break;

			vector<std::string> data;

			while (getline(iss, token, delimiter))
			{
				data.push_back(token.c_str());
			}
			year = data[0];
			data.push_back("0"); //Valeur par défaut s'il n'y a aucun flag. 
			// tester si la date de l'observation appartient à l'intervalle d'années désiré
			int annee = std::stoi(year);

			if (!(std::find(ValidYears.begin(), ValidYears.end(), annee) != ValidYears.end()))
			{
				ligne += "," + data[Headings.size()];
				output << ligne << endl;
				continue;
			}

			month = (data[1].size() == 1) ? "0" + data[1] : data[1];
			day = (data[2].size() == 1) ? "0" + data[2] : data[2];
			hour = atoi(data[3].c_str());
			boost::gregorian::date  dt_time;    // Date au format  aaaammdd
			try
			{
				dt_time = boost::gregorian::date_from_iso_string(year + month + day);
			}
			catch (boost::gregorian::bad_day_of_month)
			{
				std::cout << "error with boost::gregorian";
			}
			//duree entre DATESTART et la date de l'observation

			boost::gregorian::date_duration obs_date = dt_time - DATESTART;

			int  obs_time = obs_date.days() * 24 + hour;  // indice de l'observation dans le vecteur des flags associé à la station

			int time_loc = obs_time;
			//test if this time_loc out of bounds:
			
			if (time_loc != HoursBetween)
			{
				bool flag_libelle = false;

				// Eliminer l'observation qui a été flaggée.
				if (Headings.find("Tair") != Headings.end())
				{
					for (int ligne : FLAG_COL_DICT["temperatures"])
					{
						flag_libelle = false;

						if (station.getQc_flags()[ligne][time_loc] != 0)
						{
							data[Headings["Tair"]] = "-999";
							if (flag_libelle == false)
							{
								(data[Headings.size()] == "0") ? data[Headings.size()] = "temp" : data[Headings.size()] = data[Headings.size()] + "+ temp ";
								flag_libelle = true;
							}
						}
					}
				}
				if (Headings.find("Tdew") != Headings.end())
				{
					for (int ligne : FLAG_COL_DICT["dewpoints"])
					{
						flag_libelle = false;

						if (station.getQc_flags()[ligne][time_loc] != 0)
						{
							data[Headings["Tdew"]] = "-999";
							if (flag_libelle == false)
							{
								(data[Headings.size()] == "0") ? data[Headings.size()] = "dew" : data[Headings.size()] = data[Headings.size()] + "+ dew ";
								flag_libelle = true;
							}
						}
					}
				}
				if (Headings.find("WndS") != Headings.end())
				{
					for (int ligne : FLAG_COL_DICT["windspeeds"])
					{
						flag_libelle = false;

						if (station.getQc_flags()[ligne][time_loc] != 0)
						{
							data[Headings["WndS"]] = "-999";
							if (flag_libelle == false)
							{
								(data[Headings.size()] == "0") ? data[Headings.size()] = "winds" : data[Headings.size()] = data[Headings.size()] + "+ winds ";
								flag_libelle = true;
							}
						}
					}
				}
				if (Headings.find("Pres") != Headings.end())
				{
					for (int ligne : FLAG_COL_DICT["slp"])
					{
						flag_libelle = false;

						if (station.getQc_flags()[ligne][time_loc] != 0)
						{
							data[Headings["Pres"]] = "-999";
							if (flag_libelle == false)
							{
								(data[Headings.size()] == "0") ? data[Headings.size()] = "Pres" : data[Headings.size()] = data[Headings.size()] + "+ Pres ";
								flag_libelle = true;
							}
						}
					}
				}
				if (Headings.find("WndD") != Headings.end())
				{
					for (int ligne : FLAG_COL_DICT["winddirs"])
					{
						flag_libelle = false;

						if (station.getQc_flags()[ligne][time_loc] != 0)
						{
							data[Headings["WndD"]] = "-999";
							if (flag_libelle == false)
							{
								(data[Headings.size()] == "0") ? data[Headings.size()] = "WndD" : data[Headings.size()] = data[Headings.size()] + "+ WndD ";
								flag_libelle = true;
							}
						}
					}
				}

                ///////////////////////////////////  Copier la ligne dans le nouveau fichier       ////////////////
				ligne = "";
				for (int i = 0; i < data.size()-1; i++)
				{
					ligne += data[i] + ",";
				}
				ligne += data[data.size() - 1];

				output << ligne << endl;
			}

		}
		output.close();
	}
}