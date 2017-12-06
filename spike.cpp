#include "spike.h"

using namespace std;
using namespace UTILS;
using namespace boost;
using namespace PYTHON_FUNCTION;
using namespace boost::numeric::ublas;


namespace INTERNAL_CHECKS
{
	void sc(CStation &station, std::vector<std::string> variable_list, std::vector<int> flag_col, boost::gregorian::date  start,
		boost::gregorian::date end, std::ofstream&  logfile, bool second)
	{

		cout << "refactor" << endl;

		int v = 0;
		for (string variable : variable_list)
		{

			varrayfloat flags = station.getQc_flags(flag_col[v]);
			int prev_flag_number = 0;
			if (second)
			{
				// count currently existing flags :
				varrayfloat dummy = flags[flags != float(0)];
				prev_flag_number = dummy.size();
			}
			CMetVar& st_var = station.getMetvar(variable);
			CMaskedArray<float>  all_filtered = apply_filter_flags(st_var);
			float reporting_resolution = reporting_accuracy(all_filtered);

			//to match IDL system - should never be called as would mean no data
			if (reporting_resolution == float(-1))  reporting_resolution = 1;

			std::vector<pair<int, int>> month_ranges = month_starts_in_pairs(start, end);
			// array de taille 12*86*2 où chaque ligne correspond à un mois
			std::vector<std::valarray<pair<int, int>>> month_ranges_years = L_reshape3(month_ranges, 12);

			varraysize good = npwhere(all_filtered.m_mask, "=", false);

			CMaskedArray<float> full_time_diffs(0, all_filtered.size());

			full_time_diffs.m_mask = all_filtered.m_mask;
			full_time_diffs.m_masked_indices = all_filtered.m_masked_indices;
			varrayfloat time_data = station.getMetvar("time").getData();

			varrayfloat dummy(full_time_diffs.size());
			varraysize indices(good.size());
			for (size_t i = 0; i < indices.size()-1; i++)
			{
				indices[i] = good[i+1];
			}
			dummy = time_data[good];
			
			full_time_diffs.fill(good, dummy[std::slice(1, good.size() - 1, 1)]);
			time_data.resize(dummy.size(), 0);
			time_data -= dummy[std::slice(0, good.size() - 1, 1)];
			full_time_diffs.add(good, time_data);
			full_time_diffs.m_data[good[good.size() - 1]] = full_time_diffs.m_data[good[0]];
			

			//develop critical values using clean values
			cout << "sort the differencing if values were flagged rather than missing" << endl;

			CMaskedArray<float> full_filtered_diffs(float(0), all_filtered.size());
			full_filtered_diffs.m_mask = all_filtered.m_mask;
			full_filtered_diffs.m_masked_indices = all_filtered.m_masked_indices;

			full_filtered_diffs.fill(good,all_filtered.compressed()[std::slice(1, all_filtered.compressed().size() - 1, 1)]);
			dummy.resize(good.size(), 0);
			dummy -= all_filtered.compressed()[std::slice(0, all_filtered.compressed().size() - 1, 1)];
			full_filtered_diffs.add(good, dummy);
			full_filtered_diffs.m_data[good[good.size() - 1]] = full_filtered_diffs.m_data[good[0]];
			
			//test all values
			varraysize good_to_uncompress = npwhere(st_var.getAllData().m_mask, "=", false);

			CMaskedArray<float> full_value_diffs(float(0), st_var.getData().size());
			full_value_diffs.m_mask = st_var.getAllData().m_mask;
			full_value_diffs.m_masked_indices = st_var.getAllData().m_masked_indices;

			full_value_diffs.fill(good_to_uncompress, st_var.getAllData().compressed()[std::slice(1, st_var.getAllData().compressed().size() - 1, 1)]);

			dummy.resize(good_to_uncompress.size(), 0);

			dummy -= st_var.getAllData().compressed()[std::slice(0, st_var.getAllData().compressed().size() - 1, 1)];
			full_value_diffs.add(good_to_uncompress, dummy);
			full_value_diffs.m_data[good_to_uncompress[good_to_uncompress.size() - 1]] = full_value_diffs.m_data[good_to_uncompress[0]];

			
			// convert to compressed time to match IDL

			varrayfloat value_diffs = full_value_diffs.compressed();
			varrayfloat time_diffs = full_time_diffs.compressed();
			varrayfloat filtered_diffs = full_filtered_diffs.compressed();

			dummy = flags[good_to_uncompress];
			flags = dummy;
			matrix<float> critical_values(9, 12);
			/*  Initialisationd de critical_values*/
			for (size_t i = 0; i < critical_values.size1(); ++i)
			{
				for (size_t j = 0; j < critical_values.size2(); ++j)
					critical_values(i, j) = Cast<float>(st_var.getMdi());
			}

			//link observation to calendar month

			varrayfloat month_locs(float(0), full_time_diffs.size());

			CMaskedArray<float> this_month_time_diff(full_time_diffs.size());

			CMaskedArray<float> this_month_filtered_diff(full_filtered_diffs.size());
			for (size_t month = 0; month < 12; month++)
			{
				for (size_t year = 0; year < int(month_ranges.size() / 12); year++)
				{
					if (year == 0)
					{
						this_month_time_diff = full_time_diffs(month_ranges_years[month][year].first,
							month_ranges_years[month][year].second);
						this_month_filtered_diff = full_filtered_diffs(month_ranges_years[month][year].first,
							month_ranges_years[month][year].second);
					}
					else
					{
						Concatenate(this_month_time_diff, full_time_diffs(month_ranges_years[month][year].first,
							month_ranges_years[month][year].second));
						Concatenate(this_month_filtered_diff, full_filtered_diffs(month_ranges_years[month][year].first,
							month_ranges_years[month][year].second));
					}
					month_locs[std::slice(month_ranges_years[month][year].first,
						month_ranges_years[month][year].second - month_ranges_years[month][year].first, 1)] = month;

				}
				for (size_t delta = 1; delta < 9; delta++)
				{
					float iqr;
					varraysize locs = ma_masked_where(this_month_time_diff, "=", float(delta));
					if (locs.size() >= 100)
					{
						dummy = this_month_filtered_diff.m_data[locs];
						iqr = IQR(dummy);
						if (iqr == float(0) && delta == 1)
							critical_values(delta - 1, month) = float(6);
						else if (iqr == float(0))
							critical_values(delta - 1, month) = Cast<float>(st_var.getMdi());
						else
						{
							critical_values(delta - 1, month) = 6 * iqr;
						}
						
						float threshold = get_critical_values<float>(dummy, float(0), 0.5, critical_values(delta - 1, month));
						if (threshold < critical_values(delta - 1, month)) critical_values(delta - 1, month) = threshold;

						cout << critical_values(delta - 1, month) << "   " << iqr << "  " << 6 * iqr << endl;
					}
				}
			}
			 dummy= month_locs[good_to_uncompress];
			 month_locs = dummy;
			//not less than 5x reporting accuracy
			for (size_t i = 0; i < critical_values.size1(); ++i)
			{
				for (size_t j = 0; j < critical_values.size2(); ++j)
				{
					if (critical_values(i, j) != Cast<float>(st_var.getMdi()))
					{

						if (critical_values(i, j) <= 5 * reporting_resolution)
						{
							critical_values(i, j) = 5 * reporting_resolution;
						}

					}
				}
			}
			
			//check hourly against 2 hourly, if <2/3 the increase to avoid crazy rejection rate
			for(int month = 0; month < 12; month++)
			{
				if (critical_values(0, month) != Cast<float>(st_var.getMdi()) && critical_values(1, month) != Cast<float>(st_var.getMdi()))
				{
					if (critical_values(0, month) / critical_values(1, month) <= 0.66)
					{
						critical_values(0, month) = 0.66*critical_values(1, month);
					}
				}
			}
			
			//get time differences for unfiltered data
			full_time_diffs.resize(st_var.getData().size());

			full_time_diffs.m_mask=st_var.getAllData().m_mask;
			full_time_diffs.m_masked_indices = st_var.getAllData().m_masked_indices;
			time_data.resize(station.getMetvar("time").getData().size(),0);
			time_data = station.getMetvar("time").getData();
			

			full_time_diffs.fillData(good_to_uncompress, time_data[std::slice(1, time_data.size() - 1, 1)]);

			dummy.resize(good_to_uncompress.size(), 0);

			dummy -= time_data[std::slice(0, time_data.size() - 1, 1)];
			full_time_diffs.addData(good_to_uncompress, dummy);
			full_time_diffs.m_data[good_to_uncompress[good_to_uncompress.size() - 1]] = full_time_diffs.m_data[good_to_uncompress[0]];
			
			time_diffs = full_time_diffs.compressed();

			// go through each difference, identify which month it is in if passes spike thresholds 

			// spikes at the beginning or ends of sections

			for (size_t t = 0; t < time_diffs.size(); t++)
			{
				float median_diff;
				int j;
				(t < 1) j = time_diffs.size() - 1 : j = t - 1;

				if (std::abs(time_diffs[j]) > 240 && std::abs(time_diffs[t]) < 3)
				{
					//10 days before but short gap thereafter

					varraysize indice = good_to_uncompress[std::slice(t + 1, good_to_uncompress.size()-(t+1), 1)];

					CMaskedArray<float> next_values = st_var.getAllData()[indice];

					good = npwhere(next_values.m_mask, "=", false);

					indice = good;
					if (good.size() > 10)
					{
						//indice.resize(good.size() - 10);
						indice = good[std::slice(0, 10, 1)];
					}
					float next_median = idl_median(next_values[indice].compressed());

					float next_diff = std::abs(value_diffs[t]);  //out of spike

					median_diff = std::abs(next_median - st_var.getAllData().m_data[good_to_uncompress[t]]);  //are the remaining onees

					if (critical_values(time_diffs[t] - 1, month_locs[t]) != Cast<float>(st_var.getMdi()))
					{
						//jump from spike > critical but average after < critical / 2
						if (	std::abs(median_diff) < (critical_values(time_diffs[t] - 1, month_locs[t])/2.) &&
								std::abs(next_diff) > critical_values(time_diffs[t] - 1, month_locs[t]))
						{
								flags[t] = 1;

						}
					}
				}
				else if (std::abs(time_diffs[j]) < 3 && std::abs(time_diffs[t]) > 240)
				{
					//10 days after but short gap before

					varraysize indice = good_to_uncompress[std::slice(0, j, 1)];

					CMaskedArray<float> prev_values = st_var.getAllData()[indice];

					good = npwhere(prev_values.m_mask, "=", false);
					indice=good;
					if (good.size() >= 10)
					{
						//indice.resize(good.size() - 10);
						indice = good[std::slice(good.size() - 10, good.size(), 1)];
					}

					CMaskedArray<float> prevalues = prev_values[indice];

					float prev_median = idl_median(prev_values[indice].compressed());

					float prev_diff = std::abs(value_diffs[j]); 

					median_diff = std::abs(prev_median - st_var.getAllData().m_data[good_to_uncompress[t]]);

					if (critical_values(time_diffs[j] - 1, month_locs[t]) != Cast<float>(st_var.getMdi()))
					{
						//jump from spike > critical but average before< critical / 2
						if (std::abs(median_diff) < critical_values(time_diffs[j] - 1, month_locs[t])/2. &&
							std::abs(prev_diff) > critical_values(time_diffs[t] - 1, month_locs[t]))
						{
							flags[t] = 1;
						}
					}
				}
			}			
			for (size_t t = 0; t < time_diffs.size(); ++t)
			{
				for (int spk_len : {1, 2, 3})
				{
					if (t >= spk_len && t < time_diffs.size() - spk_len)
					{
						//check if time differences are appropriate, for multi-point spikes
						if ((std::abs(time_diffs[t - spk_len]) <= spk_len * 3) &&
							(std::abs(time_diffs[t]) <= spk_len * 3) &&
							(time_diffs[t - spk_len - 1] - 1 < spk_len * 3) &&
							(time_diffs[t + 1] - 1 < spk_len * 3) &&
							((spk_len == 1) ||
							((spk_len == 2) && (std::abs(time_diffs[t - spk_len + 1]) <= spk_len * 3)) ||
							((spk_len == 3) && (std::abs(time_diffs[t - spk_len + 1]) <= spk_len * 3) && (std::abs(time_diffs[t - spk_len + 2]) <= spk_len * 3))))
						{
							//check if differences are valid   
							if ((value_diffs[t - spk_len] != Cast<float>(st_var.getMdi())) &&
								(value_diffs[t - spk_len] != st_var.getFdi()) &&
								(critical_values(time_diffs[t - spk_len] - 1, month_locs[t]) != Cast<float>(st_var.getMdi())))
							{
								//if exceed critical values
								if (std::abs(value_diffs[t - spk_len]) >= critical_values(time_diffs[t - spk_len] - 1, month_locs[t]))
								{
									//are signs of two differences different
									if (std::copysignf(1, value_diffs[t]) != copysignf(1, value_diffs[t - spk_len]))
									{
										//are within spike differences small
										if ((spk_len == 1) ||
											((spk_len == 2) &&(std::abs(value_diffs[t - spk_len + 1]) < critical_values(time_diffs[t - spk_len + 1] - 1, month_locs[t]) / 2.)) ||
											((spk_len == 3) && (std::abs(value_diffs[t - spk_len + 1]) < critical_values(time_diffs[t - spk_len + 1] - 1,month_locs[t]) / 2.) &&
											(std::abs(value_diffs[t - spk_len + 2]) < critical_values(time_diffs[t - spk_len + 2] - 1, month_locs[t] / 2.))))
										{
											//check if following value is valid
											if (value_diffs[t] != Cast<float>(st_var.getMdi()) &&
												critical_values(time_diffs[t] - 1, month_locs[t]) != Cast<float>(st_var.getMdi()) &&
												value_diffs[t] != st_var.getFdi())

											{
												//and if at least critical value 
												if (std::abs(value_diffs[t]) >= critical_values(time_diffs[t] - 1, month_locs[t]))
												{
													//test if surrounding differences below 1/2 critical value

													if (std::abs(value_diffs[t - spk_len - 1]) <= critical_values(time_diffs[t - spk_len - 1] - 1, month_locs[t]) / 2.)
													{
														if (std::abs(value_diffs[t + 1]) <= critical_values(time_diffs[t + 1] - 1, month_locs[t]) / 2.)
														{
															//set the flags
															flags[std::slice(t - spk_len + 1, spk_len , 1)] = 1;

														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			station.setQc_flags(flags, good_to_uncompress, flag_col[v]);

			varraysize flag_locs = npwhere(station.getQc_flags()[flag_col[v]], "!", float(0));

			print_flagged_obs_number(logfile, "Spike", variable, flag_locs.size() - prev_flag_number); // additional flags

			//copy flags into attribute
			st_var.setFlags(flag_locs, float(1));

			v++;
		}
		append_history(station, "Spike check");
	}
}