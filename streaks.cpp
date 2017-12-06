#include "streaks.h"



using namespace std;
using namespace PYTHON_FUNCTION;
using namespace UTILS;

namespace INTERNAL_CHECKS
{


	inline float rsc_get_straight_string_threshold(CMetVar& st_var, boost::gregorian::date   start, boost::gregorian::date   end, float reporting , float old_threshold)
	{
		
		CMaskedArray<float> all_filtered = apply_filter_flags(st_var);
		//find and count the length of all repeating strings
		float prev_value = Cast<float>(st_var.getMdi());
		vector<int> this_string;
		vector<int> string_length;
		//run through all obs
		int o = 0;
		for (float obs : all_filtered.m_data)
		{
			if (all_filtered.m_mask[o] == false)
			{
				if (obs != prev_value)
				{
					//if different values to before
					string_length.push_back(this_string.size());
					this_string.clear();
					this_string.push_back(o);
				}
				else
				{
					// if same value as before, note and continue
					this_string.push_back(o);
				}
				prev_value = obs;
			}
			o++;
		}
		float threshold = get_critical_values(string_length, 1, 1, old_threshold);
		return  threshold;
		 
	}
	
	void rsc_annual_string_expectance(CMaskedArray<float>& all_filtered, std::vector<int>& value_starts, std::vector<int>& value_lengths, varrayfloat& flags, boost::gregorian::date  start, boost::gregorian::date end, CMetVar& st_var, const  varrayfloat& times)
	{
		vector<int> month_starts;
		UTILS::month_starts(start, end,month_starts);  //  FONCTION RESHAPE
		vector<valarray<int>> m_month_starts = L_reshape(month_starts, 12);
		valarray<float> year_proportions(Cast<float>(st_var.getMdi()), m_month_starts.size());
		CMaskedArray<float> year(Cast<float>(st_var.getMdi()), all_filtered.size());
		// churn through each year in turn
		for (int y = 0; y < m_month_starts.size(); y++)
		{
			if (y != m_month_starts.size() - 1)
			{
				/*varraysize indices(m_month_starts[y + 1][0] - m_month_starts[y][0]);
				for (size_t i = 0; i < indices.size(); i++)
				{
					indices[i] = m_month_starts[y][0] + i;
				}*/
				year = all_filtered(m_month_starts[y][0], m_month_starts[y+1][0]);
			}
			else
			{
				year = all_filtered(m_month_starts[y][0], all_filtered.size());
			}

			if (year.compressed().size() >= 200)
			{
				varraysize string_starts(value_starts.size());
				//if there are strings(streaks of same value) in this year
				if (y != m_month_starts.size() - 1)
					string_starts = npwhereAnd(value_starts, m_month_starts[y][0], m_month_starts[y + 1][0], ">=", "<");
				else
					string_starts = npwhere(value_starts, ">=", m_month_starts[y][0]);

				year_proportions[y] = 0;
				if (string_starts.size() >= 1)
				{
					//work out the proportion of the amount of data
					year_proportions[y] = Somme(value_lengths)/ float(year.compressed().size());
				}
			}
		}
		//if enough dirty years
		varraysize good_years = npwhere(year_proportions, "!", Cast<float>(st_var.getMdi()));
		
		if (good_years.size() >= 10)
		{
			varrayfloat dummy = year_proportions[good_years];
			float median = idl_median(dummy);
			dummy.free();
			if (median < 0.005) median = 0.005;
			//find the number which have proportions > 5 x median
			varraysize bad_years = npwhere(year_proportions, ">", 5 * median);
			if (bad_years.size()>1)
			{
				varraysize locs(value_starts.size());
				for (size_t bad : bad_years)
				{
					//and flag
					if (bad == m_month_starts.size() - 1)
					{
						//if last year, just select all
						locs = npwhere(value_starts, ">=", m_month_starts[bad][0]);
					}
					else
						locs = npwhereAnd(value_starts, m_month_starts[bad][0], m_month_starts[bad + 1][0], ">=", "<=");
					for (size_t loc : locs)
					{
						//need to account for missing values here
						valarray<bool> dummy = all_filtered.m_mask[std::slice(value_starts[loc], all_filtered.m_mask.size() - value_starts[loc], 1)];
						varraysize goods = npwhere(dummy, "=", false);
						goods= goods[std::slice(0, value_lengths[loc], 1)];
						goods += value_starts[loc];
						flags[goods] = 1;
					}
				}
			}
		}			

	}

	varrayfloat rsc_straight_strings(CMetVar& st_var,const varrayfloat& times, int n_obs, int n_days, boost::gregorian::date  start, boost::gregorian::date end, std::map<float, float> WIND_MIN_VALUE, bool wind, float reporting, bool dynamic)
	{
		float threshold;
		
		CMaskedArray<float> all_filtered = apply_filter_flags(st_var);
		if (st_var.getName() == "winddirs")
		{
			//remove calm periods for this check
			CMetVar wd_st_var = st_var;
			varraysize calms = npwhere(st_var.getData(), "=", float(0));  //True calms have direction set to 0, northerlies to 360
			wd_st_var.setData(calms,Cast<float>( wd_st_var.getMdi()));

			if (dynamic)
			{
				threshold = rsc_get_straight_string_threshold(wd_st_var, start, end, reporting, float(n_obs));
				if (threshold < n_obs)  n_obs = int(threshold);
			}
			all_filtered = apply_filter_flags(wd_st_var); // calms have been removed
		}
		else
		{
			if (dynamic)
			{
				threshold = rsc_get_straight_string_threshold(st_var, start, end, reporting, float(n_obs));
				if (threshold < n_obs)  n_obs = int(threshold);
			}
		}
		varrayfloat flags(0., all_filtered.size());
		//Look for continuous straight strings

		float prev_value = Cast<float>(st_var.getMdi());
		vector<int> string_points;
		//storage for excess over years
		vector<int> values_starts;
		vector<int> values_lengths;
		

		int o = 0;
		for (float obs : all_filtered.m_data)
		{
			if (all_filtered.m_mask[o] == false)
			{
				if (obs != prev_value)
				{
					if (st_var.getName() != "winddirs" || prev_value != 0)
					{//if different value to before, which is long enough (and large enough for Wind)
						if (string_points.size() >= 10)
						{
							if (wind == false || (wind == true && prev_value > WIND_MIN_VALUE[reporting]))
							{
								//note start and length for the annual excess test
								values_starts.push_back(string_points[0]);
								values_lengths.push_back(string_points.size());

								int time_diff = times[string_points[string_points.size() - 1]] - times[string_points[0]];
								
								//if length above threshold and spread over sufficient time frame, flag

								if (string_points.size() >= n_obs || time_diff >= n_days * 24) //measuring time in hours
								{
									valarray<bool> dummy(false, flags.size());
									for (int val: string_points)
									{
										dummy[val] = true;
									}
									flags[dummy] = 1;
								}
							}
						}
					}
					string_points.clear();
					string_points.push_back(o);
				}
				else
				{
					//if same value as before, note and continue
					string_points.push_back(o);
				}
				prev_value = obs;
			}
			o++;
		}
		rsc_annual_string_expectance(all_filtered,values_starts, values_lengths, flags, start, end, st_var, times);

		return flags;
	}

	void rsc_whole_day_repeats(vector<CMaskedArray<float>>& data, int n_wday, varrayfloat& day_flags)
	{
		
		valarray<varrayfloat> flags(data[0].m_data, data.size());
		vector<int > ndays;
		int nday;
		varrayfloat prev_day(data[0].size());
		for (int day = 0; day < data[0].size(); ++day)
		{
			if (day == 0)
			{
				prev_day = data[day].m_data;
				nday = 1;
				continue;
			}
			else
			{
				varraysize matches = npwhere(prev_day, data[day].m_data, "=");
				//if this day matches previous one(not IDL wording, but matches logic)
				if (matches.size() == 24)
				{
					if (nday >= n_wday)
					{
						if (data[day].compressed().size()>0)
						{
							//if above threshold and not because all values are missing, flag
							for (int i = day - nday; i < day; ++i)
								flags[i] = 1.;
						}
					}
				}
				else
				{
					ndays.push_back(nday);
					prev_day = data[day].m_data;
					nday = 1;
				}
				nday++;
			}
		}

		for (size_t i = 0; i < flags.size(); ++i)
			for (size_t j = 0; j < 24; j++)
			{
				day_flags[i + j] = flags[i][j];
			}

	}

	void rsc_hourly_repeats(CMetVar& st_var, const varrayfloat& times, int n_hrs, int n_wday, varrayfloat& flags, varrayfloat&  day_flags)
	{
		
		CMaskedArray<float> hourly_data = apply_filter_flags(st_var);
		vector<CMaskedArray<float>> h_hourly_data = C_reshape(hourly_data, 24);
		vector<valarray<float>> hourly_times = C_reshape(times, 24);
		int count = 1;
	
		for (size_t hour = 0; hour < 24; hour++)
		{
			float match_values = -999.;
			vector<int> len_matches;
			vector<float> match_times;
			
			for (int day = 0; day < h_hourly_data.size(); day++)
			{
				
				cout << "Iteration " << count++<< " / " << 24 * h_hourly_data.size() << endl;
				//for each day at each given hour
				if (h_hourly_data[day].m_mask[hour] == false)
				{
					if (h_hourly_data[day][hour] != match_values)
					{
						//if different value, check if string/streak above threshold
						if (match_times.size() > n_hrs)
						{
							varraysize bad = npwhere(match_times, "=", times);
							flags[bad] = 1;
						}
						len_matches.push_back(match_times.size());
						match_values = h_hourly_data[day][hour];
						match_times.clear();
						match_times.push_back(hourly_times[day][hour]);
					}
					else
					{
						//if same value
						match_times.push_back(hourly_times[day][hour]);
					}
				}	
			}
				
		}
		rsc_whole_day_repeats(h_hourly_data, n_wday, day_flags);
			
	}
	void rsc(CStation& station, std::vector<std::string> var_list, std::vector<vector<int>>  flag_col, boost::gregorian::date  start,
		boost::gregorian::date end, std::ofstream& logfile)
	{
		/*
		threshold values for low, mid and high resolution for each of the variables tested
		# consecutive repeats (ignoring missing data); spread over N days; same at each hour of day for N days; full day repeats
		*/

		array<int,4> tab1 {{40, 14, 25, 10 }};
		array<int,4> tab2 {{ 30, 10, 20, 7 }};
		array<int,4> tab3 {{ 24, 7, 15, 5 }};
		array<int,4> tab4 {{ 80, 14, 25, 10 }};
		array<int,4> tab5 {{ 60, 10, 20, 7 }};
		array<int,4> tab6 {{ 48, 7, 15, 2 }};
		array<int,4> tab7 {{ 120, 28, 25, 10 }};
		array<int,4> tab8 {{ 100, 21, 20, 7 }};
		array<int,4> tab9 {{ 72, 14, 15, 5 }};
		array<int,4> tab10{{ 120, 28, 28, 10 } };
		array<int,4> tab11 {{ 96, 28, 28, 10 }};
		array<int,4> tab12 {{ 72, 21, 21, 7 }};
		array<int,4> tab13	{{ 48, 14, 14, 7 }};
		array<int, 4> tab14{{ 24, 7, 14, 5 } };

		map<float ,array<int,4> > T, D, S, WS, WD;
		T[1] = tab1;
		T[0.5] = tab2;
		T[0.1] = tab3;
		D[1] =	tab4;
		D[0.5] = tab5;
		D[0.1] = tab6;
		S[1] = tab7;
		S[0.5] = tab8;
		S[0.1] = tab9;
		WS[1] = tab1;
		WS[0.5] = tab2;
		WS[0.1] = tab3;
		WD[90] = tab10;
		WD[45] = tab11;
		WD[22] = tab12;
		WD[10] = tab13;
		WD[1] =	tab14;

		map<string, map<float, array<int, 4> >> limits_dict;

		limits_dict["temperatures"] = T;
		limits_dict["dewpoints"] = D;
		limits_dict["slp"] = S;
		limits_dict["windspeeds"] = WS;
		limits_dict["winddirs"] = WD;

		map<float, float> WIND_MIN_VALUE;
		WIND_MIN_VALUE[1] = 0.5;
		WIND_MIN_VALUE[0.5] = 1.0;
		WIND_MIN_VALUE[0.1] = 0.5;

		cout << "Streaks" << endl;
		varrayfloat times = station.getMetvar("time").getData();
		int v = 0;
		for (string variable : var_list)
		{
			CMetVar& st_var = station.getMetvar(variable);
			if (apply_filter_flags(st_var).compressed().size() > 0)
			{
				bool wind = false;
				if (variable == "windspeeds") wind = true;
				bool winddir = false;
				if (variable == "winddirs") winddir = true;

				float reporting_resolution = reporting_accuracy(apply_filter_flags(st_var), winddir);
				array<int,4> limits = limits_dict[variable][reporting_resolution];

				// need to apply flags to st_var.flags each time for filtering
								
				station.setQc_flags(rsc_straight_strings(st_var, times, limits[0], limits[1], start, end,WIND_MIN_VALUE, wind, reporting_resolution, true), flag_col[v][0]);

				varrayfloat flags(st_var.getData().size());
				varrayfloat day_flags(st_var.getData().size());

				rsc_hourly_repeats(st_var, times, limits[2], limits[3], flags,day_flags);

				station.setQc_flags(flags, flag_col[v][1]);
				station.setQc_flags(day_flags, flag_col[v][2]);

				for (int streak_type = 0; streak_type < 3; ++streak_type)
				{
					varraysize flag_locs = npwhere(station.getQc_flags(flag_col[v][streak_type]), "!", float(0));

					print_flagged_obs_number(logfile, "Streak Check", variable, flag_locs.size());

					//copy flags into attribute
					st_var.setFlags(flag_locs, 1.);
				}
			}
			v++;
		}
		append_history(station, "Streak Check");
	}
}

