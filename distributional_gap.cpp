#include "distributional_gap.h"

using namespace std;
using namespace UTILS;
using namespace boost;
using namespace PYTHON_FUNCTION;
using namespace boost::numeric::ublas;


namespace
{
	//# mean vs median
	//# SD vs IQR
	//IMAGELOCATION = "blank"
	//# need to explain all of these
	//static const double OBS_LIMIT = 50;
	//static const double VALID_MONTHS = 5;
	//static const double SPREAD_LIMIT = 2;
	//static const double MONTH_LIMIT = 60;
	//static const double LARGE_LIMIT = 5; // IQR
	//static const double MEAN = false;
	//static const double NON_GAP_FLAG = false;
	//static const double FREQUENCY_THRESHOLD = 0.1;
	//static const double GAP_SIZE = 2;
	//static const double BIN_SIZE = 0.25;
}
//#************************************************************************
namespace INTERNAL_CHECKS
{
	float dgc_get_monthly_averages( CMaskedArray<float>& data, int limit, double mdi, bool MEAN )
	{
		float average = mdi;
		varrayfloat Compress = data.compressed();
		if (Compress.size() >= limit)
		{
			if (MEAN)
				average = Compress.sum() / Compress.size();
			else
				average = idl_median(Compress);
		}
		return average;
	}
	float dgc_get_monthly_averages(varrayfloat& Compress, int limit, double mdi, bool MEAN)
	{
		float average = mdi;
		if (Compress.size() >= limit)
		{
			if (MEAN)
				average = Compress.sum() / Compress.size();
			else
				average = idl_median(Compress);
		}
		return average;
	}
	//#************************************************************************
	int dgc_find_gap(varrayfloat& hist, varrayfloat& bins, float threshold, int gap_size)
	{
		size_t start = size_t(np_argmax(hist));
		bool positive;
		if (bins[start] < threshold) positive = true;
		else positive = false;

		int n = 0;
		int	gap_length = 0;
		int	gap_start = 0;

		while (true)
		{
			if (hist[start + n] == 0)
			{
				gap_length += 1;
				if (gap_start == 0)
				{
					// plus 1 to get upper bin boundary
					if (positive && bins[start + n + 1] >= threshold)
						gap_start = bins[start + n + 1];
					else if (!positive && bins[start + n] <= threshold)
						gap_start = bins[start + n];
				}
			}

			else
			{
				if (gap_length < gap_size)
					gap_length = 0;

				else if (gap_length >= gap_size && gap_start != 0)
					break;
			}

			if (start + n == hist.size() - 1 || (start + n == 0))
				break;

			if (positive)  n += 1;
			else n -= 1;
		}
		return gap_start;
	}


	void dgc_monthly(CStation& station, string variable, varrayfloat& flags, gregorian::date start, gregorian::date end)
	{


		CMetVar& st_var = station.getMetvar(variable);

		std::vector<pair<int, int>> month_ranges = month_starts_in_pairs(start, end);

		// get monthly averages
		CMaskedArray<float> month_average(Cast<float>(st_var.getMdi()), month_ranges.size());

		CMaskedArray<float> month_average_filtered(Cast<float>(st_var.getMdi()), month_ranges.size());

		CMaskedArray<float>  all_filtered = apply_filter_flags(st_var);
		int m = 0;
		for (pair<int, int>  month : month_ranges)
		{
			
			CMaskedArray<float> data = st_var.getAllData()(month.first,month.second);

			CMaskedArray<float> filtered = all_filtered(month.first, month.second );

			month_average.m_data[m] = dgc_get_monthly_averages(data, OBS_LIMIT, Cast<float>(st_var.getMdi()), MEAN);

			month_average_filtered.m_data[m] = dgc_get_monthly_averages(filtered, OBS_LIMIT, Cast<float>(st_var.getMdi()), MEAN);
			m += 1;
		}

		// get overall monthly climatologies - use filtered data

		std::vector<std::valarray<float>> list_month_average = C_reshape(month_average.m_data, 12);
		std::vector<std::valarray<float>>  list_month_average_filtered = L_reshape(month_average_filtered.m_data, 12);

		/***Initialisation de standardised_months*/

		std::valarray<varrayfloat> standardised_months(list_month_average.size());
		// Initialisation
		for (size_t i = 0; i < standardised_months.size(); i++)
		{
			varrayfloat dummy(Cast<float>(st_var.getMdi()), 12);
			standardised_months[i]=dummy;
		}
		for (size_t m = 0; m < 12; m++)
		{
			varraysize valid_filtered = npwhere(list_month_average_filtered[m], "!", Cast<float>(st_var.getMdi()));

			if (valid_filtered.size() >= VALID_MONTHS)
			{
				WBSF::CStatistic valid_data;
				
				
				for (size_t i = 0; i < valid_filtered.size(); i++)
				{
					valid_data += list_month_average_filtered[valid_filtered[i]][m];
				}
				
				float clim;
				float spread;
				if (MEAN)
				{
					clim = valid_data[WBSF::MEAN];
					spread = valid_data[WBSF::STD_DEV];
				}
				else
				{
					clim = valid_data[WBSF::MEDIAN];

					if (spread <= float(SPREAD_LIMIT))
						spread = float(SPREAD_LIMIT);

					spread = valid_data[WBSF::INTER_Q];

					if (spread <= SPREAD_LIMIT)  spread = SPREAD_LIMIT;
						
				}

				for (size_t i = 0; i < valid_filtered.size(); i++)
				{
					standardised_months[valid_filtered[i]][m] = (list_month_average[valid_filtered[i]][m] - clim)/spread;
				}
				
			}
		}
		CMaskedArray<float> new_standardised_months = Shape(standardised_months);
		
		varraysize good_months = npwhere(new_standardised_months, "!", Cast<float>(st_var.getMdi()));

		// must be able to do this with masked arrays
		

		// remove all months with a large standardised offset

		if (good_months.size() >= MONTH_LIMIT)
		{
			new_standardised_months.m_fill_value = Cast<float>(st_var.getMdi());
			new_standardised_months.masked(Cast<float>(st_var.getMdi()));

			varraysize large_offsets = npwhere(new_standardised_months ,">=", float(LARGE_LIMIT));

			if (large_offsets.size() > 0)
			{
				for (size_t lo : large_offsets)
				{
					flags[std::slice(month_ranges[lo].first, month_ranges[lo].second - month_ranges[lo].first,1)] = 1;
				}

				
			}
			// walk distribution from centre and see if any assymetry

			varrayfloat standardised_months_sort = new_standardised_months.m_data[good_months];

			// initialize original index locations
			varraysize sort_order(standardised_months_sort.size());
			std::iota(std::begin(sort_order), std::end(sort_order), 0);

			// sort indexes based on comparing values in v
			sort(std::begin(sort_order), std::end(sort_order),
				[&standardised_months_sort](size_t i1, size_t i2) {return standardised_months_sort[i1] < standardised_months_sort[i2]; });

			std::sort(std::begin(standardised_months_sort), std::end(standardised_months_sort));
			
			size_t mid_point = (good_months.size()) / 2;

			bool good = true;
			int iter = 1;
			while (good)
			{
				if (standardised_months_sort[mid_point - iter] != standardised_months_sort[mid_point + iter])
				{
					//using IDL notation
					varrayfloat tempvals{ std::abs(standardised_months_sort[mid_point - iter]), std::abs(standardised_months_sort[mid_point + iter]) };
					
					if (tempvals.min() != 0)
					{
						if ((tempvals.max() / tempvals.min()) >= 2. && tempvals.min() >= 1.5)
						{
							// substantial asymmetry in distribution - at least 1.5 from centre and difference of 2.
							varraysize bad = good_months[sort_order];
							if (tempvals[0] == tempvals.max())
							{
								// LHS
								varraysize dum = bad[std::slice(0, mid_point - iter, 1)];
								bad = dum;
							}
							else if (tempvals[1] == tempvals.max())
							{
								//RHS
								varraysize dum = bad[std::slice(mid_point + iter, bad.size() - mid_point - iter, 1)];
								bad = dum;
							}

							for (size_t b : bad)
							{
								flags[std::slice(month_ranges[b].first, month_ranges[b].second - month_ranges[b].first, 1)] = 1;
							}

							good = false;
						}
					}
				}

				iter += 1;
				if (iter == mid_point) break;
					
			}

		}
		cout << endl;
	}

	//#************************************************************************
	varrayfloat dgc_all_obs(CStation& station, string variable, varrayfloat& flags, gregorian::date start, gregorian::date end, bool windspeeds, bool GH)
	{

		// '''RJHD addition working on all observations'''


		CMetVar& st_var = station.getMetvar(variable);

		std::vector<pair<int, int>> month_ranges = month_starts_in_pairs(start, end);

		std::vector<std::valarray<std::pair<int, int>>> month_ranges_years = L_reshape3(month_ranges, 12);

		CMaskedArray<float> all_filtered = apply_filter_flags(st_var);
		
		for (int month = 0; month < 12; month++)
		{
			CMaskedArray<float> windspeeds_month;
			float windspeeds_month_average;
			float windspeeds_month_mad;
			if (windspeeds == true)
			{
				CMetVar& st_var_wind = station.getMetvar("windspeeds");
				
				// get monthly averages
				CMaskedArray<float> windspeeds_month(st_var_wind.getData().size());
				int y = 0;
				for (pair<int,int> year : month_ranges_years[month])
				{
					if (y == 0)
						windspeeds_month = st_var_wind.getAllData()(year.first,year.second);
					else
						Concatenate(windspeeds_month, st_var_wind.getAllData()(year.first, year.second));
					y++;
				}
				windspeeds_month.masked(st_var_wind.getAllData().m_fill_value);
				windspeeds_month_average = dgc_get_monthly_averages(windspeeds_month, OBS_LIMIT, Cast<float>(st_var_wind.getMdi()), MEAN);

				windspeeds_month_mad = mean_absolute_deviation(windspeeds_month.compressed(),true);
			}

			std::vector<CMaskedArray<float>> this_month_data;
			std::vector<CMaskedArray<float>> this_month_filtered;
			std::vector<int> year_ids_dummy;
			 
			valarray<int> datacount_dummy = concatenate_months(month_ranges_years[month], st_var.getAllData(), this_month_data, year_ids_dummy, false);

			datacount_dummy = concatenate_months(month_ranges_years[month], all_filtered, this_month_filtered, year_ids_dummy, false);

			varrayfloat Compress = CompressedMatrice(this_month_filtered);
			WBSF::CStatisticEx this_month_filtered_data;

			for (size_t i = 0; i < Compress.size(); i++)
			{
				this_month_filtered_data += Compress[i];
			}
			float monthly_median;
			float iqr;

			if (Compress.size() > OBS_LIMIT)
			{
				monthly_median = this_month_filtered_data[WBSF::MEDIAN];

				iqr = this_month_filtered_data[WBSF::INTER_Q];

				if (iqr == 0.0)
				{
					// to get some spread if IQR too small                   
					iqr = IQR(Compress, 0.05);
					cout << "Spurious_stations file not yet sorted" << endl;;
				}
				
				if (iqr != 0.0)
				{
					
					CMaskedArray<float> monthly_values(CompressedMatrice(this_month_data));
					monthly_values.m_fill_value = this_month_data[0].m_fill_value;
					monthly_values.masked(this_month_data[0].m_fill_value);

					monthly_values.m_data -= monthly_median;
					monthly_values.m_data /= iqr;
					
					
					
					varrayfloat Compress_month_values = monthly_values.compressed();

					varrayfloat bins, bincenters, plot_bincenters;

					create_bins(Compress_month_values, float(BIN_SIZE/10), bins, plot_bincenters);

					create_bins(Compress_month_values, BIN_SIZE, bins, bincenters);

					varrayfloat binEdges(bins);
					
					varrayfloat hist = PYTHON_FUNCTION::histogram(Compress_month_values, binEdges);

					WBSF::CStatistic stat;
					for (size_t i = 0; i < Compress_month_values.size(); i++)
					{
						stat += Compress_month_values[i];
					}

					float u_minimum_threshold, l_minimum_threshold;

					varrayfloat plot_gaussian;
					if (GH)
					{
						// Use Gauss-Hermite polynomials to add skew and kurtosis to Gaussian fit - January 2015 ^RJHD

					}
					else
					{
						if (month == 8)
							cout << endl;
						varrayfloat gaussian = fit_gaussian(bincenters, hist, hist.max(),stat[WBSF::MEAN],stat[WBSF::STD_DEV]);

						// assume the same threshold value
						u_minimum_threshold = 1 + round(invert_gaussian(FREQUENCY_THRESHOLD, gaussian));
						l_minimum_threshold = -u_minimum_threshold;

						plot_gaussian = Gaussian(plot_bincenters, gaussian);
					}
					

					size_t uppercount = npwhere(Compress_month_values,">", u_minimum_threshold).size();
					size_t lowercount = npwhere(Compress_month_values, "<", l_minimum_threshold).size();
					
					// this needs refactoring - but lots of variables to pass in
					int gap_start;
					CMaskedArray<float> this_year_data;
					varrayfloat this_year_flags;
					varraysize gap_cleaned_locations;
					if (uppercount > 0)
					{
						gap_start = dgc_find_gap(hist, binEdges, u_minimum_threshold);
						

						if (gap_start != 0)
						{

							for (pair<int, int> year : month_ranges_years[month])
							{
								this_year_data = all_filtered(year.first, year.second);

								this_year_flags = flags[std::slice(year.first, year.second - year.first, 1)];

								gap_cleaned_locations = npwhere(((this_year_data.m_data - monthly_median) / iqr), ">", float(gap_start));
								
								this_year_flags[gap_cleaned_locations] = 1;
								varraysize idx(year.second - year.first);
								std::iota(std::begin(idx), std::end(idx), year.first);
								flags[idx] = this_year_flags;
								
							}	
						}
					}
					
					
					if (lowercount > 0)
					{
						gap_start = dgc_find_gap(hist, binEdges, l_minimum_threshold);

						if (gap_start != 0)
						{

							for (pair<int, int> year : month_ranges_years[month])
							{
								

								this_year_data = all_filtered(year.first, year.second);

								this_year_flags = flags[std::slice(year.first, year.second - year.first, 1)];

								gap_cleaned_locations = npwhereAnd(((this_year_data.m_data - monthly_median) / iqr),"<", float(gap_start), this_year_data.m_mask,"!",true);
								
								this_year_flags[gap_cleaned_locations] = 1;
								varraysize idx(year.second - year.first);
								std::iota(std::begin(idx), std::end(idx), year.first);
								flags[idx] = this_year_flags;

								if (windspeeds)
								{
									this_year_flags[gap_cleaned_locations] = 2; // tentative flags

									//float slp_average = dgc_get_monthly_averages(CompressedMatrice(this_month_data), OBS_LIMIT, Cast<float>(st_var.getMdi()), MEAN);

									//float slp_mad = mean_absolute_deviation(CompressedMatrice(this_month_data), true);

									//varraysize storms = npwhere(((windspeeds_month.m_data - windspeeds_month_average) / windspeeds_month_mad),">", float(4.5));
									//	/*(this_month_data.c - slp_average) / slp_mad) > 4.5));*/

									//if (storms.size() >= 2)
									//{
									//	varraysize storm_1diffs = npDiff(storms);
									//	varraysize separations = npwhere(storm_1diffs,"!",size_t(1));
									//}
								}
							}
						}
						//for sep in separations:
					}
				}
			}
		}

		return flags; // dgc_all_obs
	}

	//#************************************************************************
	void dgc(CStation& station, std::vector<std::string> variable_list, std::vector<int> flag_col, boost::gregorian::date start,
		boost::gregorian::date  end, std::ofstream& logfile, bool GH)
	{

		int v = 0;
		for (string variable : variable_list)
		{
			
			varrayfloat flags= station.getQc_flags(flag_col[v]);
			dgc_monthly(station, variable,flags, start, end);
			station.setQc_flags(flags, flag_col[v]);

			if (variable == "slp")
			{// need to send in windspeeds too     
				flags = station.getQc_flags(flag_col[v]);
				dgc_all_obs(station, variable, flags, start, end, true, false);
				station.setQc_flags(flags, flag_col[v]);
			}
			else
			{
				flags = station.getQc_flags(flag_col[v]);
				dgc_all_obs(station, variable, flags, start, end, false, false);
				station.setQc_flags(flags, flag_col[v]);;
			}

			varraysize flag_locs = npwhere(station.getQc_flags()[flag_col[v]],"!", float(0));

			
			print_flagged_obs_number(logfile, "Distributional Gap", variable, flag_locs.size());

			// copy flags into attribute

			CMetVar& st_var = station.getMetvar(variable);
			st_var.setFlags(flag_locs,1);

			//# MATCHES IDL for 030660-99999, 2 flags in T, 30-06-2014

			v++;
		}
		append_history(station, "Distributional Gap Check");
	}
}
