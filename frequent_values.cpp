#include "frequent_values.h"

using namespace std;
using namespace UTILS;
using namespace PYTHON_FUNCTION;

namespace INTERNAL_CHECKS
{
	void  fvc(CStation& station, std::vector<std::string> variable_list, std::vector<int>  flag_col, boost::gregorian::date start, boost::gregorian::date end, std::ofstream &logfile)
	{

		const int MIN_DATA_REQUIRED = 500;// to create histogram for complete record
		const int MIN_DATA_REQUIRED_YEAR = 100; // to create histogram
		std::vector<pair<int,int>> month_ranges = month_starts_in_pairs(start, end);
		/*  creer un array de 3 dimensions à partir de month_ranges*/
		vector<vector<pair<int, int>>> month_ranges_years;
		int iteration = 1;
		vector<pair<int, int>> month;
		int index = 0;
		for (const pair<int, int> month_it : month_ranges)
		{
			if (iteration <= 12)
			{
				month.push_back(make_pair(month_it.first, month_it.second));
				iteration++;
			}
			else
			{
				month_ranges_years.push_back(month);
				index++;
				month.clear();
				month.push_back(make_pair(month_it.first, month_it.second));
				iteration = 2;
			}

		}
		month_ranges_years.push_back(month);
	
		int v = 0;
		for (string variable : variable_list)
		{
			CMetVar & st_var = station.getMetvar(variable);
			CMaskedArray<float> filtered_data = apply_filter_flags(st_var);
			float  reporting_accuracy = UTILS::reporting_accuracy(filtered_data);

			//apply flags - for detection only   (filtered_data)
			CMaskedArray<float> season_data;
			int thresholds[3];
			for (int season : {0, 1, 2, 3, 4})
			{
				if (season == 0)
				{
					//all year
					season_data = PYTHON_FUNCTION::ma_masked_values<float>(filtered_data.compressed(), st_var.getFdi());
					thresholds[0] = 30;
					thresholds[1] = 20;
					thresholds[2] = 10;
				}
				else
				{
					thresholds[0] = 20;
					thresholds[1] = 15;
					thresholds[2] = 10;

					season_data.resize(0);
				
					
					for (const vector<pair<int, int>> year : month_ranges_years)
					{
								//churn through months extracting data, accounting for fdi and concatenating together

						/****** Indices où extraire les données selon les saisons *****/
						
						// Saison 0
						varraysize indices0(year[11].second - year[0].first);
						for (size_t i = 0; i < indices0.size(); i++)
						{
							indices0[i] = year[0].first + i;
						}
						// Saison 1
						
						varraysize indices1(year[4].second - year[2].first);

						for (size_t i = 0; i < indices1.size(); i++)
						{
							indices1[i] = year[2].first + i;
						}
						// Saison 2
						varraysize indices2(year[7].second - year[5].first);

						for (size_t i = 0; i < indices2.size(); i++)
						{
							indices2[i] = year[5].first + i;
						}
						// Saison 3
						varraysize indices3(year[10].second - year[8].first);

						for (size_t i = 0; i < indices3.size(); i++)
						{
							indices3[i] = year[8].first + i;
						}
						// Saison 4
						varraysize indices40(year[1].second - year[0].first);

						for (size_t i = 0; i < indices40.size(); i++)
						{
							indices40[i] = year[0].first + i;
						}
						varraysize indices41(year[11].second - year[11].first);

						for (size_t i = 0; i < indices41.size(); i++)
						{
							indices41[i] = year[11].first + i;
						}
			/***************************************************************************/
						if (season == 1) // mars,avril,mai
						{
							
							CMaskedArray<float> new_filtered_data = filtered_data[indices1];

							new_filtered_data = ma_masked_values(new_filtered_data, st_var.getFdi());

							Concatenate(season_data , new_filtered_data); 
						}
						else if (season == 2) //june, july, august
						{
							
							CMaskedArray<float> new_filtered_data = filtered_data[indices2];

							new_filtered_data = ma_masked_values(new_filtered_data, st_var.getFdi());

							Concatenate(season_data, new_filtered_data);

						}
						else if (season == 3) // september,october,november
						{
							
							CMaskedArray<float> new_filtered_data = filtered_data[indices3];
							
							new_filtered_data = ma_masked_values(new_filtered_data, st_var.getFdi());

							Concatenate(season_data, new_filtered_data);
						}
						else //december + january,februay
						{
							
							CMaskedArray<float> new_filtered_data = filtered_data[indices40];
						
							new_filtered_data = ma_masked_values(new_filtered_data, st_var.getFdi());
							Concatenate(season_data, new_filtered_data);
														
							new_filtered_data = filtered_data[indices41];
							
							new_filtered_data = ma_masked_values(new_filtered_data, st_var.getFdi());
							Concatenate(season_data, new_filtered_data);
						}
					}
				}
				varrayfloat season_datas = season_data.compressed();
				
				varrayfloat seven_bins;
				if (season_datas.size() > MIN_DATA_REQUIRED)
				{
					varrayfloat bins, bincenters;

					if (0 < reporting_accuracy  && reporting_accuracy <= 0.5)  //-1 used as missing value
						UTILS::create_bins(season_datas, 0.5, bins, bincenters);
					else
						UTILS::create_bins(season_datas, 1, bins, bincenters);

					varrayfloat binEdges(bins);
					varrayfloat hist = PYTHON_FUNCTION::histogram(season_datas, binEdges);

					varrayfloat bad_bin(0., hist.size());
					//scan through bin values and identify bad ones
					int e = 0;
					for (float element : hist)
					{
						if (e > 3 && e <= hist.size() - 3) //don't bother with first three or last three bins
						{
							seven_bins = hist[std::slice(e - 3, 7, 1)];
							if (seven_bins[3] == seven_bins.max() && seven_bins[3] != 0) //is local maximum and != zero
							{
								if (seven_bins[3] / seven_bins.sum() >= 0.5 && seven_bins[3] >= thresholds[0])
								{
									//contains >50 % of data and is greater than threshold
									bad_bin[e] = 1;
								}
							}
						}
						e++;
					}
					//having identified possible bad bins, check each year in turn
					
					CMaskedArray<float> year_data;
					varrayfloat year_flags;
					for(const vector<pair<int, int>> year : month_ranges_years)
					{
						/***************** Indices où extraire les données selon les saisons **************/
						// Saison 0
						varraysize indices0(year[11].second - year[0].first);
						for (size_t i = 0; i < indices0.size(); i++)
						{
							indices0[i] = year[0].first + i;
						}
						// Saison 1

						varraysize indices1(year[4].second - year[2].first);

						for (size_t i = 0; i < indices1.size(); i++)
						{
							indices1[i] = year[2].first + i;
						}
						// Saison 2
						varraysize indices2(year[7].second - year[5].first);

						for (size_t i = 0; i < indices2.size(); i++)
						{
							indices2[i] = year[5].first + i;
						}
						// Saison 3
						varraysize indices3(year[10].second - year[8].first);

						for (size_t i = 0; i < indices3.size(); i++)
						{
							indices3[i] = year[8].first + i;
						}
						// Saison 4
						varraysize indices40(year[1].second - year[0].first);

						for (size_t i = 0; i < indices40.size(); i++)
						{
							indices40[i] = year[0].first + i;
						}
						varraysize indices41(year[11].second - year[11].first);

						for (size_t i = 0; i < indices41.size(); i++)
						{
							indices41[i] = year[11].first + i;
						}

						/***************************************************************/
						if (season == 0)
						{
							
							year_data.resize(indices0.size());
							year_data = st_var.getAllData()[indices0];
							year_data = ma_masked_values(year_data, st_var.getFdi());

							year_flags.resize(indices0.size());
							year_flags = station.getQc_flags()[flag_col[v]][indices0];
						}
						if (season == 1) // mars,avril,may
						{
							year_data.resize(indices1.size());
							year_data = st_var.getAllData()[indices1];
							year_data = ma_masked_values(year_data, st_var.getFdi());

							year_flags.resize(indices1.size());
							year_flags = station.getQc_flags()[flag_col[v]][indices1];
						}
						else if (season == 2) //june, july, august
						{

							year_data.resize(indices2.size());
							year_data = st_var.getAllData()[indices2];
							year_data = ma_masked_values(year_data, st_var.getFdi());

							year_flags.resize(indices2.size());
							year_flags = station.getQc_flags()[flag_col[v]][indices2];
						}
						else if (season == 3) // september,october,november
						{

							year_data.resize(indices3.size());
							year_data = st_var.getAllData()[indices3];
							year_data = ma_masked_values(year_data, st_var.getFdi());

							year_flags.resize(indices3.size());
							year_flags = station.getQc_flags()[flag_col[v]][indices3];
						}
						else if (season == 4)//december + january,februay
						{
							
							year_data.resize(indices40.size());
							year_data = st_var.getAllData()[indices40];
							year_data = ma_masked_values(year_data, st_var.getFdi());

							CMaskedArray<float> dum = st_var.getAllData()[indices41];
							dum = ma_masked_values(dum, st_var.getFdi());
							
							Concatenate(year_data, dum);
							
							//year_flags.resize(indices.size());
							Concatenate(year_flags,station.getQc_flags()[flag_col[v]][indices40]);
							Concatenate(year_flags,station.getQc_flags()[flag_col[v]][indices41]);

						}
						if (year_data.compressed().size() > MIN_DATA_REQUIRED_YEAR)
						{
							hist = PYTHON_FUNCTION::histogram(year_data.compressed(), binEdges);
							int e = 0;
							for (float element : hist)
							{
								if (bad_bin[e] == 1) //only look at pre-identified bins
								{
									if (e >= 3 && e <= hist.size() - 3) //don't bother with first three or last three bins
									{
										seven_bins = hist[std::slice(e - 3, 7, 1)];
										if (seven_bins[3] == seven_bins.max() && seven_bins[3] != 0) //is local maximum and != zero
										{
											if ((seven_bins[3] / seven_bins.sum() >= 0.5 && seven_bins[3] >= thresholds[1])
												|| (seven_bins[3] / seven_bins.sum() >= 0.9 && seven_bins[3] >= thresholds[2]))
											{
												//contains  >50% or >90% of data and is greater than appropriate threshold
												// Falg these data
												bad_bin[e] = 1;
												//np_where avec deux conditions
												std::valarray<bool> bad_points(false, year_flags.size());
												for (int i = 0; i < year_data.size(); i++)
												{
													if (year_data[i] >= binEdges[e] && year_data[i] < binEdges[e + 1])
														bad_points[i] = true;
												}
												year_flags[bad_points] = 1;
											}
										}
									}
								}
								e++;
							}

						}
						//copy flags back
			
						if (season == 0)
						{
							
							station.setQc_flags(year_flags, indices0, flag_col[v]);
						}
						else if (season == 1)
							station.setQc_flags(year_flags, indices1, flag_col[v]);
						else if (season == 2)
							station.setQc_flags(year_flags, indices2, flag_col[v]);
						else if (season == 3)
							station.setQc_flags(year_flags, indices3, flag_col[v]);
						else if (season == 4)
						{
							
							station.setQc_flags(year_flags[indices40], indices40, flag_col[v]);
							station.setQc_flags(year_flags[std::slice(indices40.size(),indices41.size(),1)], indices41, flag_col[v]);
						}
					}
					
				}
			}
			
			varraysize flag_locs = npwhere<float>(station.getQc_flags()[flag_col[v]], "!",0);
			UTILS::print_flagged_obs_number(logfile, "Frequent value", variable, flag_locs.size());
			//copy flags into attribute
			st_var.setFlags(flag_locs, 1);
			v++;
		}
		UTILS::append_history(station, "Frequent values Check");
	}
}