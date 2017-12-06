#include "climatological.h"


using namespace std;
using namespace UTILS;
using namespace boost;
using namespace PYTHON_FUNCTION;

namespace INTERNAL_CHECKS
{

	float coc_get_weights(const CMaskedArray<float>& monthly_vqvs, varraysize monthly_subset, varraysize filter_subset)
	{
		float weights;
		varrayfloat filterweights{ 1., 2., 3., 2., 1. };
		varrayfloat dummy = monthly_vqvs.m_data[monthly_subset];
		dummy -= monthly_vqvs.m_data[monthly_subset].apply(MyApplyRoundFunc);
		dummy = dummy.apply(MyApplyCeilFunc);
		dummy *= filterweights[filter_subset];

		if (dummy.sum() == 0)
		{
			weights = 0;
		}
		else
		{
			varrayfloat dum = filterweights[filter_subset];
			dum *= monthly_vqvs.m_data[monthly_subset];
			weights = dum.sum() / dummy.sum();	
		}

		return weights;
	}
	void coc_low_pass_filter(valarray<CMaskedArray<float>>& normed_anomalies,const vector<int>& year_ids,const CMaskedArray<float>& monthly_vqvs, int years)
	{
		varraysize monthly_range;
		varraysize filter_range;

		for (int year = 0; year < years;year++)
		{
			if (year == 0)
			{
				monthly_range = Arange(3, 0);
				filter_range = Arange(5, 2);
			}
			else if (year == 1)
			{
				monthly_range = Arange(4, 0);
				filter_range = Arange(5, 1);
			}
			else if (year == years - 2)
			{
				monthly_range = Arange(monthly_vqvs.size(), -4 + monthly_vqvs.size());
				filter_range = Arange(4, 0);
			}
			else if (year == years - 1)
			{
				monthly_range = Arange(monthly_vqvs.size(), -3 + monthly_vqvs.size());
				filter_range = Arange(3, 0);
			}
			else
			{
				monthly_range = Arange(year + 3, year - 2);
				filter_range = Arange(5);
			}
			
			varrayfloat dummy = monthly_vqvs.m_data[monthly_range].apply(MyApplyAbsFunc);

			if (dummy.sum() != 0)

			{
				float weights = coc_get_weights(monthly_vqvs, monthly_range, filter_range);
	
				varraysize year_locs = npwhere(year_ids,"=", year);

				for (size_t i = 0; i < year_locs.size(); i++)
				{
					varrayfloat indata = normed_anomalies[year_locs[i]].m_data;
					indata -= weights;
					normed_anomalies[year_locs[i]].m_data = indata;
				}
			}
		}
	}

	
	void coc_find_and_apply_flags(std::valarray<std::pair<int, int>>& month_ranges, valarray<CMaskedArray<float>>& normed_anomalies, varrayfloat& flags, const vector<int>& year_ids, float threshold, int gap_start, bool upper)
	{
		int y = 0;
		for (pair<int, int> year : month_ranges)
		{
			varraysize year_locs = npwhere(year_ids, "=", y);
			valarray<CMaskedArray<float>> this_year_data(year_locs.size());
			for (size_t i = 0; i < year_locs.size(); i++)
			{
				this_year_data[i] = normed_anomalies[year_locs[i]]; 
			}
			varrayfloat dummy = flags[std::slice(year.first, year.second - year.first, 1)];
			vector<varrayfloat> this_year_flags = C_reshape(dummy, 24);
			valarray<pair<size_t,size_t>> tentative_locations;

			if (upper)
			{
				tentative_locations = np_ma_where(this_year_data, ">", threshold);
			}
			else
			{
				
				tentative_locations = np_ma_where(this_year_data, "<", -threshold);
			}
			
			for (size_t i = 0; i < tentative_locations.size(); i++)
			{
				this_year_flags[tentative_locations[i].first][tentative_locations[i].second] = 2;
			}

			if (gap_start != 0)
			{
				valarray<pair<size_t, size_t>> gap_cleaned_locations;
				if (upper)
					gap_cleaned_locations = np_ma_where(this_year_data, ">", float(gap_start));
				else 
					gap_cleaned_locations = np_ma_where(this_year_data, "<", float(gap_start));

				

				for (size_t i = 0; i < gap_cleaned_locations.size(); i++)
				{
					this_year_flags[gap_cleaned_locations[i].first][gap_cleaned_locations[i].second] = 1;
				}		
			}
			varrayfloat year_flag = Shape(this_year_flags);
			flags[std::slice(year.first, year.second - year.first, 1)] = year_flag;
			y++;

		}
	}
		
	void coc(CStation& station, std::vector<std::string> variable_list, std::vector<int> flag_col, boost::gregorian::date  start, boost::gregorian::date end, std::ofstream&  logfile, bool idl )
	{
		int v = 0;
		for (string variable : variable_list)
		{
			CMetVar& st_var = station.getMetvar(variable);
			CMaskedArray<float> all_filtered = apply_filter_flags(st_var);

			vector<pair<int, int>> month_ranges = month_starts_in_pairs(start, end);
			// array de taille 12*43*2 où chaque ligne correspond à un mois
			std::vector<std::valarray<pair<int, int>>> month_ranges_years = L_reshape3(month_ranges, 12);		

			for (int month = 0; month < 12; ++month)
			{
				varrayfloat  hourly_climatologies(Cast<float>(st_var.getMdi()), 24);
				//append all e.g.Januaries together

				vector<int> year_ids;
				varrayInt datacount(month_ranges.size());
				vector<CMaskedArray<float>> this_month;
				vector<CMaskedArray<float>> this_month_filtered;
				datacount = concatenate_months(month_ranges_years[month], st_var.getAllData(), this_month_filtered, year_ids, Cast<float>(st_var.getMdi()), true);
				year_ids.clear();
				datacount = concatenate_months(month_ranges_years[month], all_filtered, this_month, year_ids, Cast<float>(st_var.getMdi()), true);


				//if fixed climatology period, sort this here


				//get hourly climatology for each month
				for (int hour = 0; hour < 24; ++hour)
				{
					CMaskedArray<float> this_hour = ExtractColumn(this_month, hour);
					//need to have data if this is going to work!
					if (this_hour.compressed().size()>0)
					{
						//winsorize & climatologies - done to match IDL
						if (idl)
						{
							varrayfloat dummy(this_hour.compressed().size() + 1);

							for (size_t i = 0; i < dummy.size() - 1; ++i)
								dummy[i] = this_hour.compressed()[i];
							dummy[dummy.size() - 1] = float(-999999);

							winsorize(dummy, 0.05, idl);
							hourly_climatologies[hour] = dummy.sum() / (dummy.size() - 1);

						}
						else
						{
							varrayfloat dummy = this_hour.compressed();
							winsorize(dummy, 0.05, idl);

							hourly_climatologies[hour] = dummy.sum() / dummy.size();
						}
					}
				}
				//Verifier si this_month contient des valeurs non masquées
				valarray<bool> v_mask(this_month.size());
				for (size_t i = 0; i < this_month.size(); ++i)
					v_mask[i] = this_month[i].m_mask.max();
				if (v_mask.max()) //equivalent de len(this_month.compressed()) > 0
				{
					//can get stations with few obs in a particular variable.

					//anomalise each hour over month appropriately
					valarray<CMaskedArray<float>> anomalies(this_month.size()), anomalies_filtered(this_month_filtered.size());

					for (size_t i = 0; i < this_month.size(); ++i)
					{
						anomalies[i] = this_month[i];
						anomalies[i].m_data -= hourly_climatologies;
						//anomalies[i].masked(this_month[i].m_fill_value);

						anomalies_filtered[i] = this_month_filtered[i];
						anomalies_filtered[i].m_data -= hourly_climatologies;

					}

					varrayfloat Compress = CompressedMatrice(anomalies);
					float iqr;
					if (Compress.size() >= 10)
					{
						iqr = IQR(Compress) / 2;
						if (iqr < 1.5)  iqr = 1.5;
					}
					else
						iqr = Cast<float>(st_var.getMdi());

					valarray<CMaskedArray<float>> normed_anomalies = anomalies;
					valarray<CMaskedArray<float>> normed_anomalies_filtered = anomalies_filtered;

					for (size_t i = 0; i < normed_anomalies.size(); i++)
					{
						normed_anomalies[i] /= iqr;
						normed_anomalies_filtered[i] /= iqr;
					}
					// get average anomaly for year
					CMaskedArray<float> monthly_vqvs(month_ranges_years[0].size());
					monthly_vqvs.m_mask = false;
					for (int year = 0; year < month_ranges_years[0].size(); year++)
					{
						varraysize year_locs = npwhere(year_ids, "=", year);
						valarray<CMaskedArray<float>> this_year = normed_anomalies_filtered[year_locs];
						varrayfloat Compress = CompressedMatrice(this_year);

						if (Compress.size() > 0)
						{
							//need to have data for this to work!

							monthly_vqvs.m_data[year] = idl_median(Compress);
						}
						else
							monthly_vqvs.m_mask[year] = true;

					}
					//low pass filter
					coc_low_pass_filter(normed_anomalies, year_ids, monthly_vqvs, int(month_ranges_years[0].size()));
					varrayfloat bins, bincenters;
					UTILS::create_bins(normed_anomalies, 1, bins, bincenters);

					varrayfloat binEdges(bins);
					Compress = CompressedMatrice(normed_anomalies);
					double mu = Compress.sum() / Compress.size();
					WBSF::CStatistic stat;
					for (size_t i = 0; i < Compress.size(); i++)
					{
						stat += Compress[i];
					}


					varrayfloat hist = PYTHON_FUNCTION::histogram(Compress, binEdges);
					varrayfloat gaussian = fit_gaussian(bincenters, hist, hist.max(), mu, stat[WBSF::STD_DEV]);

					float minimum_threshold = std::round(1. + invert_gaussian(FREQUENCY_THRESHOLD, gaussian));

					size_t uppercount = np_ma_where(normed_anomalies, ">", minimum_threshold).size();
					size_t lowercount = np_ma_where(normed_anomalies, "<", -minimum_threshold).size();

					varrayfloat these_flags = station.getQc_flags(flag_col[v]);

					int gap_start = dgc_find_gap(hist, binEdges, minimum_threshold, 1);

					coc_find_and_apply_flags(month_ranges_years[month], normed_anomalies, these_flags, year_ids, minimum_threshold, gap_start, true);
					
					gap_start = dgc_find_gap(hist, binEdges, -minimum_threshold, 1);
					coc_find_and_apply_flags(month_ranges_years[month], normed_anomalies, these_flags, year_ids, minimum_threshold, gap_start, false);

					station.setQc_flags(these_flags, flag_col[v]);

				}
			}
			varraysize flag_locs = npwhere(station.getQc_flags(flag_col[v]), "!", float(0));
			//copy flags into attribute
			st_var.setFlags(flag_locs, float(1));

			print_flagged_obs_number(logfile, "Climatological", variable, flag_locs.size());
			logfile << "where " << endl;
			size_t nflags = npwhere(station.getQc_flags(flag_col[v]), "=", float(1)).size();

			print_flagged_obs_number(logfile, "  Firm Clim", variable, nflags);

			nflags = npwhere(station.getQc_flags(flag_col[v]), "=", float(2)).size();

			print_flagged_obs_number(logfile, "  Tentative Clim", variable, nflags);

			v++;
		}

		append_history(station, "Climatological Check");
	}
}