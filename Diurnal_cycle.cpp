#include "Diurnal_cycle.h"


using namespace std;
using namespace UTILS;
using namespace PYTHON_FUNCTION;

namespace INTERNAL_CHECKS
{


	bool dcc_quartile_check( CMaskedArray<float>& day)
	{
		valarray<int> quartile_has_data(0, 4);
		varraysize quart_start; //quart_end;
		quart_start = PYTHON_FUNCTION::Arange<int>(HOURS, 0, int(HOURS / 4.));
		/*varraysize quart_end = quart_start;
		quart_end += 6;*/
		for (int q : {0, 1, 2, 3})
		{
			varraysize indice{ quart_start[q], quart_start[q] + 1, quart_start[q] + 2, quart_start[q] + 3, quart_start[q] + 4, quart_start[q] + 5 };
			
			if (day[indice].compressed().size() > 0)
				quartile_has_data[q] = 1;
		}
		if (quartile_has_data.sum() >= 3) return true;
		else return false;
	}

	


	void dcc(CStation& station, vector<string> variable_list, vector<string>& full_variable_list, vector<int> flag_col,
		ofstream &logfile)
	{
		//list of flags for each variable
		vector<varrayfloat> diurnal_flags;
		int v = 0;
		for (string variable : variable_list)
		{
			CMetVar& st_var = station.getMetvar(variable);
			CMaskedArray<float> filtered_data = UTILS::apply_filter_flags(st_var);
			vector<CMaskedArray<float>> v_filtered_data = C_reshape(filtered_data, 24);

			int number_of_days = v_filtered_data.size();

			valarray<int> diurnal_best_fits(INTMDI, number_of_days);
			valarray<int> diurnal_uncertainties(INTMDI, number_of_days);

			int d = 0;
			for (CMaskedArray<float> day : v_filtered_data)
			{

				//enough observations and have large enough diurnal range
				if (day.compressed().size() > OBS_PER_DAY)
				{
					float obs_daily_range = day.compressed().max() - day.compressed().min();
					if (obs_daily_range >= DAILY_RANGE)
					{
						varrayfloat scaled_sine, diffs;
						if (dcc_quartile_check(day))
						{
							scaled_sine = dcc_make_sine().apply([](float val){ return ((val + 1) / 2); });
							scaled_sine *= obs_daily_range;
							scaled_sine += day.compressed().min();
							diffs.resize(HOURS, 0.);
							//Find differences for each shifted sine --> cost function
							for (int h = 0; h < HOURS; h++)
							{

								varrayfloat data = (day - scaled_sine).compressed();
								data = data.apply([](float val){ return std::abs(val); });
								diffs[h] = data.sum(); // check this array
								PYTHON_FUNCTION::np_roll<float>(scaled_sine, 1);
							}
							diurnal_best_fits[d] = PYTHON_FUNCTION::np_argmin(diffs);
							//default uncertainty is the average time resolution of the data
							diurnal_uncertainties[d] = int(round(float(HOURS) / day.size()));
							if (DYNAMIC_DIURNAL)
							{
								float critical_value = diffs.min() + ((diffs.max() - diffs.min())*0.33);
								//centre so minimum in middle

								PYTHON_FUNCTION::np_roll<float>(diffs, 11 - diurnal_best_fits[d]);
								int  uncertainty = 1;
								while (uncertainty < 11)
								{
									if (diffs[11 - uncertainty] > critical_value && diffs[11 + uncertainty] > critical_value)
									{
										//break if both sides greater than critical difference when counting outwards
										break;
									}
									uncertainty++;
								}
								//check if uncertainty greater than time resolution for day
								if (uncertainty > diurnal_uncertainties[d])
									diurnal_uncertainties[d] = uncertainty;
							}
						}

					}
				}
				d++;
			}
			//For each uncertainty range (1-6h) find median of cycle offset
			valarray<int> best_fits(-9, 6);
			for (int h = 0; h < 6; h++)
			{
				varraysize locs = npwhere<int>(diurnal_uncertainties, "=", h + 1);
				if (locs.size() > 300)
				{
					valarray<int> fits = diurnal_best_fits[locs];
					best_fits[h] = UTILS::idl_median(fits);  //A tester avec de nouvelles données
				}
			}
			//Build up range of cycles incl, uncertainty to find where best of best located
			valarray<int> hours = PYTHON_FUNCTION::arange<int>(24);
			valarray<int> hour_matches(0, 24);
			int diurnal_peak = -9;
			int number_estimates = 0;
			for (int h = 0; h < 6; h++)
			{
				if (best_fits[h] = !- 9)
				{
					//Store lowest uncertainty best fit as first guess
					if (diurnal_peak == -9)
					{
						diurnal_peak = best_fits[h];
						PYTHON_FUNCTION::np_roll<int>(hours, 11 - diurnal_peak);
						hour_matches[slice(11 - (h + 1), 2 * h + 4, 1)] = 1;
						number_estimates += 1;
					}
					varraysize centre = PYTHON_FUNCTION::npwhere<int>(hours, "=", best_fits[h]);
					if ((centre[0] - h + 1) >= 0)
					{
						if ((centre[0] + h + 1) <= 23)
						{
							//A revoir
							valarray<int> dummy = hour_matches[slice((centre[0] - (h + 1)), 2 * h + 4, 1)];
							dummy += 1;
							hour_matches[slice((centre[0] - (h + 1)), 2 * h + 4, 1)] = dummy;
							dummy.free();
						}
						else
						{
							valarray<int> dummy = hour_matches[slice((centre[0] - (h + 1)), hour_matches.size() - (centre[0] - (h + 1)), 1)];
							dummy += 1;
							hour_matches[slice((centre[0] - (h + 1)), hour_matches.size() - (centre[0] - (h + 1)), 1)] = dummy;
							dummy = hour_matches[slice(0, (centre[0] + h + 2 - 24), 1)];
							dummy += 1;
							hour_matches[slice(0, (centre[0] + h + 2 - 24), 1)];
							dummy.free();
						}
					}
					else
					{
						valarray<int> dummy = hour_matches[slice(0, (centre[0] + h + 2), 1)];
						dummy += 1;
						hour_matches[slice(0, (centre[0] + h + 2), 1)] = dummy;
						dummy = hour_matches[slice((centre[0] - (h + 1)), hour_matches.size() - (centre[0] - (h + 1)), 1)];
						dummy += 1;
						hour_matches[slice((centre[0] - (h + 1)), hour_matches.size() - (centre[0] - (h + 1)), 1)] = dummy;
						dummy.free();
					}
					number_estimates += 1;
				}
			}
			//If value at lowest uncertainty not found in all others, then see what value is found by all others
			if (hour_matches[11] != number_estimates)
			{
				varraysize all_match = npwhere<int>(hour_matches, "=", number_estimates);
				//if one is, then use it
				if (all_match.size() > 0)
					diurnal_peak = all_match[0];
				else
					diurnal_peak = -9;
			}
			// Now have value for best fit diurnal offset
			valarray<int> potentially_spurious(INTMDI, number_of_days);
			varrayfloat diurnal_flags(filtered_data.size());
			if (diurnal_peak != -9)
			{
				valarray<int> hours = PYTHON_FUNCTION::arange<int>(24);
				PYTHON_FUNCTION::np_roll<int>(hours, 11 - diurnal_peak);
				for (int d = 0; d < number_of_days; d++)
				{
					if (diurnal_best_fits[d] != INTMDI)
					{
						/*Checks if global falls inside daily value+/-range
							rather than seeing if each day falls in global value + / -range'''*/
						int min_range = 11 - diurnal_uncertainties[d];
						int max_range = 11 + diurnal_uncertainties[d];
						int maxloc = (PYTHON_FUNCTION::npwhere<int>(hours, "=", diurnal_best_fits[d]))[0];

						if (maxloc < min_range || maxloc > max_range)
							potentially_spurious[d] = 1;
						else
							potentially_spurious[d] = 0;
					}
				}
				//count number of good, missing and not-bad days
				int n_good = 0;
				int	n_miss = 0;
				int	n_not_bad = 0;
				int	total_points = 0;
				int	total_not_miss = 0;
				valarray<int> to_flag(0, number_of_days);
				for (int d = 0; d < number_of_days; d++)
				{
					if (potentially_spurious[d] == 1)
					{
						n_good = 0;
						n_miss = 0;
						n_not_bad = 0;
						total_points += 1;
						total_not_miss += 1;
					}
					else
					{
						if (potentially_spurious[d] == 0)
						{
							n_good += 1;
							n_not_bad += 1;
							if (n_miss != 0)
								n_miss = 0;
							total_not_miss += 1;
						}
						if (potentially_spurious[d] == -999)
						{
							n_miss += 1;
							n_not_bad += 1;
							if (n_good != 0)
								n_good = 0;
							total_not_miss += 1;
						}
						if ((n_good == 3) || (n_miss == 3) || (n_not_bad >= 6))
						{
							if (total_points >= 30)
								if (float(total_not_miss) / total_points >= 0.5)
								{
									valarray<int> dummy = to_flag[slice(d - total_points, total_points, 1)];
									dummy += 1;
									to_flag[slice(d - total_points, total_points, 1)] = dummy;
									dummy.free();
								}
							n_good = 0;
							n_miss = 0;
							n_not_bad = 0;
							total_points = 0;
							total_not_miss = 0;
						}
					}
				}


				for (int d = 0; d < number_of_days; ++d)
				{
					if (to_flag[d] == 1)
					{
						varraysize good = npwhere(v_filtered_data[d].m_mask, "=", false);

						if (good.size() > 0)
						{
							good += 24 * d;

							diurnal_flags[good] = 1;
						}
					}
				}

				cout << "currently matches IDL, but should all hours in days have flags set, not just the missing/flagged ones?" << endl;
			}
			station.setQc_flags(diurnal_flags, flag_col[v]);

			varraysize flag_locs = npwhere(station.getQc_flags(flag_col[v]), "!", float(0));

			print_flagged_obs_number(logfile, "Diurnal Cycle", variable, flag_locs.size());

			//copy flags into attribute
			st_var.setFlags(flag_locs, float(1));

			// A Tester avec la station   030660-99999, 30-06-2014, 855 flagged RJHD
		}
		UTILS::apply_flags_all_variables(station, full_variable_list, 0, logfile, "Diurnal Cycle");
		UTILS::append_history(station, "Diurnal Cycle Check");
	}
}