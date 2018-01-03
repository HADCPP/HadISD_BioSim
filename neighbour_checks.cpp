#include "neighbour_checks.h"

using namespace std;
using namespace UTILS;
using namespace boost;
using namespace boost::numeric::ublas;
using namespace boost::gregorian;
using namespace PYTHON_FUNCTION;


namespace NEIGHBOUR_CHECKS
{
	void get_distances_angles(CStation station, const std::vector<CStation>& station_info, ivector& distances, ivector& angles)
	{
		//Return two big symmetrical arrays of the station separations and bearings
		
		pair<float, float> coord1 = make_pair(station.getLat(), station.getLon());
		int st = 0;
		for (CStation stat : station_info)
		{
			if (stat != station)
			{
				pair<int, int> result = get_dist_and_bearing(coord1, make_pair(stat.getLat(), stat.getLon()));
				distances(st) = result.first;
				angles(st) = result.second;
				st++;
			}
			else
			{
				distances(st) = 0;
				angles(st) = 0;
				st++;
			}
		} 

		
	}

	void get_all_neighbours(int station_loc,float st_elev, ivector& distances, ivector& bearings, fvector& elevations, ivector& neighbours, ivector& neighbour_quadrants, float sep_limit, float elev_limit, int max_neighbours)
	{
		varrayInt quadrants(float(0), bearings.size());
		int a = 0;
		for (int angle : {270, 180, 90, 0})
		{
			varraysize locs = npwhereAnd(bearings, angle, angle + 90, ">=", "<");
			
			quadrants[locs] = a + 1;
		}
		// start by sorting by distance
		varraysize ordering = arg_sort(bearings);
		
		ivector neighbour_locations;
		

		for (int index : ordering)
		{
			//if within range of distance and elevation
			if (index == station_loc) continue;
			int i = 0;
			if (distances[index] <= sep_limit)
			{
				if (std::abs(elevations[index] - st_elev) <= elev_limit)

				{
					neighbour_locations[i] = index;
					neighbour_quadrants[i] = quadrants[index];
					i++;
				}
			}
		}
		int stop;
		(neighbour_locations.size() <= max_neighbours) ? stop = neighbour_locations.size() : stop = max_neighbours;

		vector_range<ivector> vr(neighbour_locations, boost::numeric::ublas::range(0, stop));
		neighbours = vr;
		vector_range<ivector> vr1(neighbour_quadrants, boost::numeric::ublas::range(0, stop));
		neighbour_quadrants = vr1;
	}
	
	//Format the time series to get allow for sensible correlations
	CMaskedArray<float> hourly_daily_anomalies(CMaskedArray<float> timeseries, int obs_per_day)
	{
		std::vector<CMaskedArray<float>> time_series = C_reshape(timeseries, 24);

		varrayfloat _daily_mean( time_series.size());
		varrayfloat not_mask_count(time_series.size());
		for (int i = 0; i < _daily_mean.size(); i++)
		{
			_daily_mean[i] = time_series[i].ma_mean();

			not_mask_count[i] = time_series[i].compressed().size();
		}
		//completeness check - only for correlations
		CMaskedArray<float> daily_mean = ma_masked_where(not_mask_count, "<", float(obs_per_day), _daily_mean,false);
		daily_mean.masked(INTMDI, true);

		std::vector<CMaskedArray<float>> hourly_anomalies; //removed annual cycle

		for (int i = 0; i < time_series.size(); i++)
		{
			CMaskedArray<float> dummy = time_series[i] - daily_mean[i];

			hourly_anomalies.push_back(dummy);
		}

		CMaskedArray<float> hourly = Shape(hourly_anomalies);

		std::vector<CMaskedArray<float>> hourly_anomalies_axis0 = L_reshape(hourly, 24);
		CMaskedArray<float> hourly_mean;
		for (int i = 0; i < 24; i++)
		{
			hourly_mean.m_data[i] = hourly_anomalies_axis0[i].ma_mean();
		}
		hourly_mean.masked(float(INTMDI));

		std::vector<CMaskedArray<float>> _anomalies = hourly_anomalies; //removed diurnal cycle

		for (int i = 0; i < _anomalies.size(); i++)
		{
			_anomalies[i] = _anomalies[i]- hourly_mean;
		}

		return Shape(_anomalies);
	}

	bool readNeighbourData(const string fichier, boost::gregorian::date DATESTART, boost::gregorian::date  DATEEND, CStation& station)
	{
		string namefile = station.getId();

		//duree entre DATASTART et DATAEND
		date end = DATEEND + years(1);
		boost::gregorian::date_duration DaysBetween = end - DATESTART;
		int HoursBetween = int(DaysBetween.days()) * 24;
		std::vector<int> TimeStamps;
		PYTHON_FUNCTION::linspace<int>(TimeStamps, 0, HoursBetween - 1, HoursBetween);
		std::vector<int> ValidYears;
		PYTHON_FUNCTION::linspace<int>(ValidYears, int(DATESTART.year()), int(DATEEND.year()), int(DATEEND.year() - DATESTART.year() + 1));
		string dubiousfile = LOG_OUTFILE_LOCS + "dubious_data_files.txt";
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
		ifstream input;
		stringstream sst;
		sst << fichier;
		string chaine = sst.str();

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
		string ligne = "";
		string token;
		char  delimiter = ',';
		map<string, int> Headings;
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
		bool sortie = false;
		
			
			if (Headings["Tair"])  temp = Headings["Tair"];
			if (Headings["Tdew"])  dew = Headings["Tdew"];
			if (Headings["Pres"])  pres = Headings["Pres"];
			if (Headings["WndS"])  windS = Headings["WndS"];
			if (Headings["WndD"])  windD = Headings["WndD"];
			if (Headings["Prcp"])  precipitation = Headings["Prcp"];
			if (Headings["RelH"])  hum = Headings["RelH"];
			

		

		while (!input.eof())
		{
			getline(input, ligne);
			istringstream iss(ligne);
			string year, month, day;
			int hour;
			if (ligne == "") break;

			std::vector<string> data;

			while (getline(iss, token, delimiter))
			{
				data.push_back(token.c_str());
			}
			year = data[0];

			// tester si la date de l'observation appartient à l'intervalle d'années désiré
			int annee = std::stoi(year);

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
				(Headings["Tair"]) ? temperatures[time_loc] = atof(data[temp].c_str()) : temperatures[time_loc] = -999;
				temperature_flags[time_loc] = INTMDI;
				(Headings["Tdew"]) ? dewpoints[time_loc] = atof(data[dew].c_str()) : dewpoints[time_loc] = -999;
				dewpoint_flags[time_loc] = INTMDI;
				(Headings["RelH"]) ? humidity[time_loc] = atof(data[hum].c_str()) : humidity[time_loc] = -999;
				humiditys_flags[time_loc] = INTMDI;
				(Headings["WndS"]) ? windspeeds[time_loc] = atof(data[windS].c_str()) : windspeeds[time_loc] = -999;
				windspeeds_flags[time_loc] = INTMDI;
				(Headings["WndD"]) ? winddirs[time_loc] = atof(data[windD].c_str()) : winddirs[time_loc] = -999;
				winddirs_flags[time_loc] = INTMDI;
				(Headings["Pres"]) ? slp[time_loc] = atof(data[pres].c_str()) : slp[time_loc] = -999;
				slp_flags[time_loc] = INTMDI;
				(Headings["Prcp"]) ? precip[time_loc] = atof(data[precipitation].c_str()) : precip[time_loc] = -999;
				precip_flags[time_loc] = INTMDI;
			}
		}


		std::vector<string> full_var_list = process_var;
		if (carry_thru_vars.size() != 0)
		{
			copy(carry_thru_vars.begin(), carry_thru_vars.end(), std::back_inserter(full_var_list));
		}
		full_var_list.push_back("time");

		for (string variable : full_var_list)
		{
			varrayfloat dummy(HoursBetween);
			// temperature  min et max
			if (variable == "temperatures")
			{
				dummy = temperatures[temperatures != float(INTMDI)];
				float max_temp;
				(Headings["Tair"]) ? max_temp = dummy.max() : max_temp = -999;
				float min_temp;
				(Headings["Tair"]) ? min_temp = dummy.min() : min_temp = -999;

				CMetVar this_var(variable, "Dry bulb air temperature at screen height (~2m)");
				this_var.setMdi(INTMDI);
				this_var.setUnits("degree_Celsius");
				this_var.setValidMin(min_temp);
				this_var.setValidMax(max_temp);
				this_var.setCoordinates("latitude longitude elevation");
				this_var.setStandard_name("surface_temperature");
				this_var.setFdi(-2.e30);

				this_var.setData(temperatures);

				//Ajouter la variable méteo à la liste des variables de la station stat
				station.setMetVar(this_var, variable);
			}
			if (variable == "dewpoints")
			{
				dummy = dewpoints[dewpoints != float(INTMDI)];
				float max_dew;
				(Headings["Tdew"]) ? max_dew = dummy.max() : max_dew = -999;
				float min_dew;
				(Headings["Tdew"]) ? min_dew = dummy.min() : min_dew = -999;

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

				//Ajouter la variable méteo à la liste des variables de la station stat
				station.setMetVar(this_var, variable);
			}
			if (variable == "windspeeds")
			{
				dummy = windspeeds[windspeeds != float(INTMDI)];
				float max_ws;
				(Headings["WndS"]) ? max_ws = dummy.max() : max_ws = -999;
				float min_ws;
				(Headings["WndS"]) ? min_ws = dummy.min() : min_ws = -999;

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

				//Ajouter la variable méteo à la liste des variables de la station stat
				station.setMetVar(this_var, variable);
			}

			if (variable == "slp")
			{
				dummy = slp[slp != float(INTMDI)];
				float smax;
				float smin;
				(Headings["pres"]) ? smax = dummy.max() : smax = -999;
				(Headings["pres"]) ? smin = dummy.min() : smin = -999;

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

				//Ajouter la variable méteo à la liste des variables de la station stat
				station.setMetVar(this_var, variable);
			}
			if (variable == "time")
			{
				CMetVar this_var(variable, "time_of_measurement");
				
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
	//Pearson product-moment correlation coefficients

	float corrcoef(CMaskedArray<float> X, CMaskedArray<float> Y)
	{
		float Xmean = X.ma_mean();
		float Ymean = Y.ma_mean();

		float Xt, Yt, Sxx=0, Syy=0, Sxy=0;

		for (int i = 0; i < X.size(); i++)
		{
			if (X.m_mask[i] == true || Y.m_mask[i] == true) continue;
			else
			{
				Xt = X[i] - Xmean;
				Yt = Y[i] - Ymean;
				Sxx += Xt*Xt;
				Syy += Yt*Yt;
				Sxy += Xt*Yt;
			}
		}
		if (Sxx == 0 || Syy == 0) return float(INTMDI);

		return Sxy / std::sqrtf(Sxx*Sxy);
	}

	/*From the list of nearby stations select the ones which will be good neighours for the test.
	Select on basis of correlation, overlap of data points and bearing (quadrants)
	*/
	void  select_neighbours(CStation station, string variable, const std::vector<CStation>& station_info, ivector neighbours,
		ivector neighbour_distances, ivector neighbour_quadrants, boost::gregorian::date start, boost::gregorian::date  end, std::ofstream& logfile)
	{

		//set up storage arrays
		varrayfloat n_correlations(float(0), neighbours.size());
		varrayfloat	n_distances(float(0), neighbours.size());
		varrayfloat	n_quadrants(float(0), neighbours.size());
		varrayfloat	n_overlaps(float(0), neighbours.size());
		varrayfloat	combined_score(float(0), neighbours.size());

		//get station data

		CMaskedArray<float> st_var_data = station.getMetvar(variable).getAllData();
		CMaskedArray<float> st_anomalies = hourly_daily_anomalies(st_var_data);

		//go through initial list and extract correlations and overlaps
		int nn = 0;
		for (int nn_loc : neighbours)
		{
			CStation neigh = station_info[nn_loc];
			/////////////////////////////////////////Lire les données de la station neigh   /////////////////////////////////////////////////////////////////////
			test internal_tests2;

			if (!DATA_READING::readData(neigh, start, end)) continue;
						
			//////////////////////////////////////////////////////////////  FIN LECTURE DATA ///////////////////////////////////////////////////

			///////Internal test  for neigh/////////////////////
			ofstream logfile;
			stringstream sst;
			sst << LOG_OUTFILE_LOCS << (neigh).getId() << ".neighbour";
			try{ logfile.open(sst.str().c_str()); }
			catch (std::exception e)
			{
				std::cout << e.what() << endl;
			}
			INTERNAL_CHECKS::internal_checks(neigh, internal_tests2, start, end, logfile);
			
			valarray<bool> dummy = create_fulltimes(neigh, {variable }, start, end, { "" });

			//get the correlations of data to this neighbour

			CMetVar neigh_var = neigh.getMetvar(variable);

			CMaskedArray<float> neigh_anomalies = hourly_daily_anomalies(neigh_var.getAllData());

			float correlation = corrcoef(neigh_anomalies, st_anomalies);

			float overlap = (npwhereOr(neigh_var.getAllData().m_mask, "=", false, st_var_data.m_mask, "=", false).size()) / (st_var_data.compressed().size());

			if (correlation != float(INTMDI))
			{
				n_correlations[nn] = correlation;
				n_overlaps[nn] = overlap;
				combined_score[nn] = correlation + overlap;
				n_quadrants[nn] = neighbour_quadrants[nn_loc];
			}

			dummy.free();
			neigh_anomalies.free();





		}
	}
	
	void neighbour_checks(CStation& station, const std::vector<CStation>& station_info, boost::gregorian::date start, boost::gregorian::date  end, std::ofstream&  logfile)
	{
		//calculate distances and angles 
		
		
		ivector distances(station_info.size() - 1);
		ivector angles(station_info.size() - 1);


		cout << " calculating distances and bearings matrix" << endl;

		get_distances_angles(station, station_info, distances, angles);

		fvector neighbour_elevations(station_info.size() );
		svector neighbour_ids(station_info.size() );

		int station_loc;
		for (int i = 0; i < station_info.size();i++)
		{
			if (station_info[i] == station) station_loc = i;
				neighbour_elevations(i) = station_info[i].getElev();
				neighbour_ids(i) = station_info[i].getId();
		}

		//process each neighbour

		
		
		cout << " Neighbour Check for station  " << station.getName() << endl;

		logfile << "Neighbour Check  \n" ;

		//return all neighbours up to a limit from the distance and elevation offsets (500km and 300m respectively)

		ivector neighbours, neighbour_quadrants;

		get_all_neighbours(station_loc,station.getElev(), distances, angles, neighbour_elevations, neighbours, neighbour_quadrants,500,300,20);

		if (neighbours.size() == 0)
		{
			logfile << " No neighbour found for the station \n    " << station.getName();
			exit(1);
		}
		logfile << " Neighbours \n  "  ;
		for (int n : neighbours)
		{
			logfile << "| " << neighbour_ids[n] << "  ----   ";
		}
		//if sufficient neighbours
		if (neighbours.size() > 3)
		{
			map<const char*, int>::const_iterator iflag = FLAG_OUTLIER_DICT.cbegin();

			for (; iflag != FLAG_OUTLIER_DICT.cend(); iflag++)
			{
				string variable = iflag->first;
				int col = iflag->second;
				CMetVar st_var = station.getMetvar(variable);

				logfile << "Length of  " << variable << "  record:  " << st_var.getAllData().compressed().size(); 

				if (st_var.getAllData().compressed().size() > 0)
				{

				}
			}
		}
		
	}
}