#include "neighbour_checks.h"

using namespace std;
using namespace UTILS;
using namespace boost;
using namespace boost::numeric::ublas;
using namespace boost::gregorian;
using namespace PYTHON_FUNCTION;


namespace NEIGHBOUR_CHECKS
{
	void get_distances_angles(CStation& station, const std::vector<CStation>& station_info, ivector& distances, ivector& angles)
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
	CMaskedArray<float> hourly_daily_anomalies(CMaskedArray<float>& timeseries, int obs_per_day)
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

	
	//Pearson product-moment correlation coefficients

	float corrcoef(CMaskedArray<float>& X, CMaskedArray<float>& Y)
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
	vector<size_t>  select_neighbours(CStation& station, string variable, const std::vector<CStation>& station_info, ivector neighbours,
		ivector neighbour_distances, ivector neighbour_quadrants, boost::gregorian::date start, boost::gregorian::date  end, std::ofstream& logfile)
	{

		//set up storage arrays
		varrayfloat n_correlations(float(0), neighbours.size());
		varrayfloat	n_distances(float(0), neighbours.size());
		varrayInt	n_quadrants(float(0), neighbours.size());
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
			try { logfile.open(sst.str().c_str()); }
			catch (std::exception e)
			{
				std::cout << e.what() << endl;
			}
			INTERNAL_CHECKS::internal_checks(neigh, internal_tests2, start, end, logfile);

			valarray<bool> dummy = create_fulltimes(neigh, { variable }, start, end, { "" });

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
			nn++;
			dummy.free();
			neigh_anomalies.free();
		}
		varraysize sort_order = arg_dsort(combined_score);

		// sort out the quadrants
		valarray<size_t> locs = neighbours[sort_order];
		valarray<int> quadrants = n_quadrants[sort_order];
		valarray<size_t> locs1 = locs[quadrants == 1];
		valarray<size_t> locs2 = locs[quadrants == 2];
		valarray<size_t> locs3 = locs[quadrants == 3];
		valarray<size_t> locs4 = locs[quadrants == 4];

		vector<size_t> final_locs;

		if (locs1.size() > 2)
		{
			final_locs.push_back(locs1[0]);
			final_locs.push_back(locs1[1]);
		}
		else if (locs1.size() == 1) final_locs.push_back(locs1[0]);

		if (locs2.size() > 2)
		{
			final_locs.push_back(locs2[0]);
				final_locs.push_back(locs2[1]);
		}
		else if (locs2.size() == 2) final_locs.push_back(locs2[0]);
		if (locs3.size() > 2)
		{
			final_locs.push_back(locs3[0]);
			final_locs.push_back(locs3[1]);
		}
		else if (locs3.size() == 3) final_locs.push_back(locs3[0]);
		if (locs4.size() > 2)
		{
			final_locs.push_back(locs4[0]);
			final_locs.push_back(locs4[1]);
		}
		else if (locs4.size() == 1) final_locs.push_back(locs4[0]);

		//and add the rest in order of combined score

		for(int index : locs)
		{
			if (!(std::find(final_locs.begin(), final_locs.end(), index) != final_locs.end()))
			{
				final_locs.push_back(index);
			}
			if (final_locs.size() == N_NEIGHBOURS) break;
		}
		
		return final_locs;
	}
	
	void detect(CStation& station, CStation&  neighbour, string  variable, varrayfloat flags, varrayfloat neighbour_count, boost::gregorian::date start, end, int distance = 0)
	{
		map<const char*, vector<int>> FILTERING_FLAG_COL = { {"temperatures",{0,1,4,5,8,12,16,20,27,41,44,58}},
		{"dewpoints" ,{[0,2,4,6,8,9,13,17,21,28,30,31,32,42,45,59}},
		{"slp" , {0,3,4,7,11,15,19,23,29,43,46,60}},
		{"windspeeds" ,{ [0,4,10,14,18,22,56,62,63,64}} };

		CMetVar st_var = station.getMetvar(variable);
		CMetVar neigh_var = neighbour.getMetvar(variable);

		//filter by flags - not all (no Climatological [24,25], or Odd cluster [54,55,56,57]), T record check not in D, 

		valarray<float> total_flags(st_var.getAllData().size());
		for (int i = 0; i < total_flags.size(); i++)
		{
			float sum = 0;
			for (int line : FILTERING_FLAG_COL[variable])
			{
				sum += station.getQc_flags[line][i];
			}
			total_flags[i] = sum;
		}
		CMaskedArray<float> st_filtered = ma_masked_where<float, float>(total_flags	, float(1),st_var.getData(), float(INTMDI));
		CMaskedArray<float> neigh_filtered = ma_masked_where<float, float>(total_flags, float(1), neigh_var.getData(), float(INTMDI));
	
		//match the observation times

		varraysize match = npwhereAnd(st_filtered.m_data, "!", float(INTMDI), neigh_filtered.m_data,"!", float(INTMDI));
		
		vector<pair<int, int>> month_starts_in_pairs(start, end);
		std::vector<std::valarray<pair<int, int>>> month_ranges_years = L_reshape3(month_ranges, 12);

		
		if (match.size() >= 100)
		{
			neighour_count[match] += 1;
			
			CMaskedArray<float> differences(INTMDI, st_filtered.size());
			differences.m_mask = true;
			varrayfloat dummy = st_filtered.m_data[match];
			differences.m_data[match] = dummy;
			dummy= neigh_filtered.m_data[match];
			differences.m_data[match] -= dummy;
			differences.m_mask[match] = false;

			varrayfloat all_iqrs(float(0), differences.size());

			//get monthly IQR values

			for (int month = 0; month < 12; month++)
			{
				float iqr;

				vector<int> year_ids;
				varrayInt datacount(month_ranges.size());
				vector<CMaskedArray<float>> this_month;
				datacount = concatenate_months(month_ranges_years[month], differences, this_month, year_ids, mdi, true);
				year_ids.clear();
				datacount.free();
				dummy = CompressedMatrice(this_month);
				if (dummy.size()> 4)
				{
					iqr = IQR(dummy);
					if (iqr <= 2) iqr = 2;
					else iqr = 2;
				}
				//and copy back into the array
				for (pair<int, int> year : month_ranges_years[month])
				{
					all_iqr[std::slice(year.first, year.second - year.first, 1)] = iqr;
				}
				
			}

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

					vector<size_t> final_neighbours = select_neighbours(station,variable, station_info, neighbours,
						neighbour_distances, neighbour_quadrants, start, end, logfile));
					//now read in final set of neighbours and process
					varrayfloat neigh_flags(float(0), st_var.getAlldata().size());// count up how many neighbours think this obs is bad
					varrayfloat neigh_count = (float(0), st_var.getAlldata().size()); //number of neighbours at each time stamp
					varrayfloat	dpd_flags = (float(0), st_var.getAlldata().size()); // number of neighbours at each time stamp
					varrayfloat	reporting_accuracies = (float(0), st_var.getAlldata().size()); // reporting accuracy of each neighbour
				
					vector<CMaskedArray<float>> all_data;
					int nn = 0;
					varrayfloat reporting_accuracies(final_neighbours.size());
					for (int nn_loc : final_neighbours)
					{
						CStation neigh = station_info[nn_loc];
						varraysize dummy= create_fulltimes(neigh, { variable }, start, end, { "" },false,true,true);
						//lecture des données de neigh
						all_data[nn] = apply_filter_flags(neigh.getMetvar(variable));

						//detect
						reporting_accuracies[nn] = reporting_accuracies(neigh.getMetvar(variable).getAllData());
					}
					//gone through all neighbours
					//if at least 2 / 3 of neighbours have flagged this point(and at least 3 neighbours)

					varraysize some_flags = npwhere(neigh_flags, ">", 0);
					varraysize outlier_locs;//= npwhereAnd(neigh_count[some_flags],">=",float(3),)

					///flag where < 3 neighbours

					varraysize locs = npwhere(neigh_count[some_flags], "<", 3);
					//station.qc_flags[some_flags[locs], col] = -1

				}
			}
		}
		
	}
}