#include "variance.h"

using namespace std;
using namespace UTILS;
using namespace boost;
using namespace PYTHON_FUNCTION;
using namespace boost::numeric::ublas;


namespace INTERNAL_CHECKS
{
	void evc(CStation& station, std::vector<std::string> variable_list, std::vector<int> flag_col, boost::gregorian::date start, boost::gregorian::date  end, std::ofstream& logfile, bool idl )
	{
		//very similar to climatological check - ensure that not duplicating
		int v = 0;
		for (string variable : variable_list)
		{
			CMetVar& st_var = station.getMetvar(variable);
			float mdi = Cast<float>(st_var.getMdi());
			float reporting_resolution = reporting_accuracy(apply_filter_flags(st_var));
			int  reporting_freq = reporting_frequency(apply_filter_flags(st_var));
			std::vector<pair<int, int>> month_ranges = month_starts_in_pairs(start, end);
			
			std::vector<std::valarray<std::pair<int, int>>> month_ranges_years = L_reshape3(month_ranges, 12);
			
			for (size_t month = 0; month < 12; ++month)
			{
				//set up hourly climatologies
				varrayfloat hourly_clims (st_var.getAllData().m_fill_value, 24);
				std::vector<int> year_ids;
				std::vector<CMaskedArray<float>> this_month;
				valarray<valarray<int>> month_data_count(month_ranges_years.size());

				for (size_t i = 0; i < month_data_count.size(); i++)
				{
					valarray<int> dummy(0, month_ranges_years[0].size());
					month_data_count[i] = dummy;
				}

				month_data_count[month]=concatenate_months(month_ranges_years[v], st_var.getAllData(), this_month, year_ids,mdi, true);

				//winsorize and get hourly climatology 
				for (size_t h = 0; h < 24; h++)
				{
					CMaskedArray<float> this_hour = ExtractColumn(this_month, h);
					//need to have data if this is going to work!
					if (this_hour.compressed().size() > 0)
					{
						//winsorize & climatologies - done to match IDL
						if (idl)
						{
							varrayfloat dummy(this_hour.compressed().size() + 1);

							for (size_t i = 0; i < dummy.size() - 1; ++i)
								dummy[i] = this_hour.compressed()[i];
							dummy[dummy.size() - 1] = float(-999999);
							winsorize(dummy, 0.05, idl);
							CMaskedArray<float> this_hour_winsorized = CMaskedArray<float>(dummy);
							hourly_clims[h] = this_hour_winsorized.ma_sum() / (this_hour_winsorized.size() - 1);

						}
						else
						{
							varrayfloat dummy = this_hour.compressed();
							winsorize(dummy, 0.05, idl);
							CMaskedArray<float> this_hour_winsorized = CMaskedArray<float>(dummy);
							this_hour_winsorized.masked(mdi);
							hourly_clims[h] = this_hour_winsorized.ma_mean();
						}
					}
				}

				valarray<CMaskedArray<float>> anomalies(this_month.size());

				for (size_t i = 0; i < this_month.size(); i++)
				{
					anomalies[i] = this_month[i];
					anomalies[i].m_data -= hourly_clims;
		
				}
				//extract IQR of anomalies(using 1 / 2 value to match IDL)
				float iqr;
				varrayfloat comp_anomalies = CompressedMatrice(anomalies);
				if (comp_anomalies.size() >= 10)
				{
					iqr = IQR(comp_anomalies) / 2.; // to match IDL
					if (iqr < 1.5)  iqr = 1.5;
				}
				else iqr = mdi;

				valarray<CMaskedArray<float>> normed_anomalies = anomalies;

				for (size_t i = 0; i < normed_anomalies.size(); i++)
				{
					normed_anomalies[i] /= iqr;
				}

				CMaskedArray<float> variances(mdi, month_ranges_years[0].size());
				variances.m_mask = false;
				varrayfloat rep_accuracies(mdi, month_ranges_years[0].size());
				varrayfloat rep_freqs(mdi, month_ranges_years[0].size());

				//extract variance of normalised anomalies for each year
				
				for (int y = 0; y < month_ranges_years[0].size(); y++)
				{
					varraysize year_locs = npwhere(year_ids, "=", y);

					valarray<CMaskedArray<float>> this_year = normed_anomalies[year_locs];
					CMaskedArray<float> this_year_1D = Shape(this_year);
					this_year_1D.masked(mdi);
					varrayfloat Compress = this_year_1D.compressed();
					if (Compress.size() >= 30)
					{
						variances.m_data[y] = mean_absolute_deviation(Compress, true);
						rep_accuracies[y] = reporting_accuracy(this_year_1D);
						rep_freqs[y] = reporting_frequency(this_year_1D);
					}
					else variances.m_mask[y] = true;
				}
				varraysize good = npwhere(month_data_count[month], ">=", 100);
				float median_variance,iqr_variance;
				//get median and IQR of variance for all years for this month

				if (good.size() > 10)
				{
					median_variance = idl_median(variances[good].compressed());
					iqr_variance = IQR(variances[good].compressed()) / 2.; // to match IDL

					if (iqr_variance < 0.01) iqr_variance = 0.01;
				}
				else
				{
					median_variance = mdi;
					iqr_variance = mdi;
				}

				//if SLP, then get median and MAD of SLP and windspeed for month

				CMetVar& winds = station.getMetvar("windspeeds");
				CMetVar& slp = station.getMetvar("slp");

				float median_wind ;
				float median_slp;

				float wind_MAD ;
				float slp_MAD;

				if (variable == "slp" || variable == "windspeeds")
				{
										//refactor this as similar in style to how target data extracted 
					CMaskedArray<float> winds_year;
					std::vector<CMaskedArray<float >> winds_month;
					CMaskedArray<float> slp_year;
					std::vector<CMaskedArray<float >> slp_month;
					for (int y = 0; y < month_ranges_years[0].size(); y++)
					{
						
						
						if (y == 0)
						{
							winds_year = winds.getAllData()(month_ranges_years[month][y].first, month_ranges_years[month][y].second);
							winds_month=C_reshape(winds_year,24);

							slp_year = slp.getAllData()(month_ranges_years[month][y].first, month_ranges_years[month][y].second);
							slp_month=C_reshape(slp_year, 24);

						}
						else
						{

							winds_year = winds.getAllData()(month_ranges_years[month][y].first, month_ranges_years[month][y].second);
							Concatenate(winds_month,C_reshape(winds_year, 24));

							slp_year = slp.getAllData()(month_ranges_years[month][y].first, month_ranges_years[month][y].second);
							Concatenate(slp_month,C_reshape(slp_year, 24));
						}
					}
					varrayfloat winds_month_compress = CompressedMatrice(winds_month);
					varrayfloat slp_month_compress = CompressedMatrice(slp_month);
					
					 median_wind = idl_median(winds_month_compress);
					 median_slp = idl_median(slp_month_compress);

					 wind_MAD = mean_absolute_deviation(winds_month_compress);
					 slp_MAD = mean_absolute_deviation(slp_month_compress);


				}
				//now test to see if variance exceeds expected range
				for (int y = 0; y < month_ranges_years.size(); y++)
				{
					if (variances[y] != mdi && iqr_variance != mdi && median_variance != mdi && month_data_count[y][month] >= DATA_COUNT_THRESHOLD)
					{
						//if SLP, then need to test if deep low pressure ("hurricane/storm") present
						//   as this will increase the variance for this month + year
						float iqr_threshold;
						if (variable == "slp" || variable == "windspeeds")
						{
							iqr_threshold = 6;
							//increase threshold if reporting frequency and resolution of this
							//  year doesn't match average
							if (rep_accuracies[y] != reporting_resolution && rep_freqs[y] != reporting_freq)
							{
								iqr_threshold = 8;
							}
							if (std::abs((variances[y] - median_variance) / iqr_variance) > iqr_threshold)
							{
								//check for storms     
								
								CMaskedArray<float> winds_month_storm = winds.getAllData()(month_ranges_years[y][month].first, month_ranges_years[y][month].second);
								CMaskedArray<float>  slp_month_storm = slp.getAllData()(month_ranges_years[y][month].first, month_ranges_years[y][month].second);

								bool storm = false;
								if (winds_month_storm.compressed().size() >= 1 && slp_month_storm.compressed().size() >=1)
								{
									// locations where winds greater than threshold
									varraysize high_winds = npwhere((winds_month_storm.m_data - median_wind) / wind_MAD,">", MAD_THRESHOLD);
									varraysize low_slps =  npwhere((median_slp - slp_month_storm.m_data) / slp_MAD ,">", MAD_THRESHOLD);

									varraysize match_loc= high_winds[in1D(high_winds,low_slps)];

									if (match_loc.size() > 0) storm = true;

								}
								else cout << " Write spurious" << endl;

								//check the SLP first difference series to ensure a drop down
								//  and climb out of minimum SLP / or climb up and down from maximum wind speed
								varrayfloat diffs;
								if (variable == "slp")
									diffs = npDiff(slp_month_storm.compressed());
								else if (variable == "windspeeds")
									diffs = npDiff(winds_month_storm.compressed());

								int negs = 0, poss = 0;
								int biggest_neg = 0, biggest_pos = 0;

								for (float diff :diffs )
								{
									if (diff > 0)
									{
										if (negs > biggest_neg)  biggest_neg = negs;
										negs = 0;
										poss += 1;
									}
									else
									{
										if (poss > biggest_pos) biggest_pos = poss;
										poss = 0;
										negs += 1;
									}
								}
								if ((biggest_neg < 10) && (biggest_pos < 10) && !storm)
								{
									// not a hurricane, so mask
									varraysize idx(month_ranges_years[y][month].second - month_ranges_years[y][month].first);
									std::iota(std::begin(idx), std::end(idx), month_ranges_years[y][month].first);

									station.setQc_flags(1, idx, flag_col[v]);
										
									logfile << "No Storm or Hurricane in  " << month + 1 << "  " << y + start.year() << "  . flagging " << endl;
								}
								else
								{
									//hurricane
									logfile << "Storm or Hurricane in  " << month + 1 << "  " << y + start.year() << "   . Not flagging " << endl;
									
								}
							}
						}
						else
						{
							iqr_threshold = 8;

							if ((rep_accuracies[y] != reporting_resolution) && (rep_freqs[y] != reporting_freq))
								iqr_threshold = 10;

							if (std::abs(variances[y] - median_variance) / iqr_variance > iqr_threshold)
							{
								//remove the data
								varraysize idx(month_ranges_years[y][month].second - month_ranges_years[y][month].first);
								std::iota(std::begin(idx), std::end(idx), month_ranges_years[y][month].first);

								station.setQc_flags(1, idx, flag_col[v]);
							}
						}
					}
				}

			}
			varraysize flag_locs = npwhere(station.getQc_flags(flag_col[v]), "!", float(0));

			print_flagged_obs_number(logfile, " Variance", variable, flag_locs.size());

			//copy flags into attribute
			st_var.setFlags(flag_locs, 1);

			v++;
		}
		append_history(station, " Excess Variance Check ");
	}
}