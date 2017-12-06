#include "humidity.h"

using namespace std;
using namespace UTILS;
using namespace boost;
using namespace PYTHON_FUNCTION;
using namespace boost::numeric::ublas;

namespace INTERNAL_CHECKS
{

	void hcc_sss(varrayfloat & T, varrayfloat& D, std::vector<pair<int, int>> month_ranges, std::ofstream& logfile, varrayfloat& flags)
	{
		
		//flag each location where D > T

		
		
		for (pair<int,int> month : month_ranges)
		{
			float data_count = 0, sss_count = 0;
			for (size_t t = month.first; t < month.second; t++)
			{
				data_count += 1;
				if (D[t] > T[t])
				{
					sss_count += 1;
					flags[t] = 1;
				}		
			}
			//test whole month
			// if more than 20% flagged, flag whole month
			if (sss_count/data_count >= SSS_MONTH_FRACTION)
				flags[std::slice(month.first, month.second - month.first, 1)] = 1;
		}
		size_t nflags = npwhere(flags, "!", float(0)).size();
		print_flagged_obs_number(logfile, "Supersaturation", "temperature", nflags);
		
	}

	varrayfloat hcc_dpd(varrayfloat& times, CMaskedArray<float>& T, CMaskedArray<float>& D, CMaskedArray<float>& P, CMaskedArray<float>& C, CMaskedArray<float>& SX, std::ofstream& logfile)
	{

		varrayfloat flags(float(0),T.size());
		CMaskedArray<float> dpds = T - D;
		float last_dpds = -9999;
		float string_start_time = times[0];
		int start_loc = 0;
		int t = 0;
		for (float tt : times)
		{
			if (tt > 0 && tt < times[times.size() - 1])
			{
				if (dpds.m_mask[t] == false)
				{
					//if change in DPD, examine previous string
					if (dpds[t] != last_dpds)
					{
						// if long enough
						int j;
						(t < 1) ? j = times.size() - 1: j = t-1;
						if (times[j] - string_start_time >= 24)
						{
							float abs_diff;
							CMaskedArray<float> these_dpds = dpds(start_loc, t);
							varraysize good = npwhere(these_dpds.m_mask, "=", false);
							if (T[t] >= 0)
								abs_diff = 0.25;
							else
								abs_diff = 1;
							//has enough data and is small enough
							if ((good.size() >= 4) && std::abs(last_dpds) <= abs_diff)
							{
								//check if weather event could explain it
								CMaskedArray<float> these_sigwx = SX(start_loc, t);
								CMaskedArray<float> these_P = P(start_loc, t);
								CMaskedArray<float> these_CB = C(start_loc, t);

								//use past significant weather, precipitation or low cloud base
								
								std::vector<bool> fog;
								for (size_t i = 0; i < good.size(); ++i)
								{
									if (these_sigwx[good][i] >= 4 || these_P[good][i]>0 || (these_CB[good][i] > 0 && these_CB[good][i] < 1000))
									{
										fog.push_back(true);
									}
								}

								if (fog.size() >= 1)
								{
									if (1.*fog.size() / good.size() < 0.333)
									{
										flags[std::slice(start_loc, t, 1)] = 1;  
									}		
								}
							}
						}
						string_start_time = tt;
						start_loc = t;
						last_dpds = dpds[t];
					}
				}
			}
			t++;
		}	
		size_t nflags = npwhere(flags,"!", float(0) ).size();
		print_flagged_obs_number(logfile, "Dewpoint Depression", "temperature", nflags);
		return flags;
	}

	varrayfloat hcc_cutoffs(CMaskedArray<float>& T, CMaskedArray<float>& D, std::vector<pair<int, int>> month_ranges, std::ofstream& logfile)
	{
		varrayfloat flags(float(0), T.size());

		int binwidth = 10;
		varraysize bins = Arange(70, 7 - 90, binwidth);
		int m = 0;
		
		for (pair<int,int> month : month_ranges)
		{
			CMaskedArray<float> this_month_T = T(month.first, month.second);
			CMaskedArray<float> this_month_D = D(month.first, month.second);

			varraysize goodT = npwhere(this_month_T.m_mask, "=", false);
			varraysize goodD = npwhere(this_month_D.m_mask, "=", false);

			//check if more than 112 obs(4 / d * 28 d)
			if (goodT.size() > 112)
			{
				//run through each bin
				for (size_t bin : bins)
				{
					varraysize bin_locs = npwhereAnd(this_month_T.m_data, float(bin), float(bin + binwidth), ">=", "<");
					// if data in this bin
					if (bin_locs.size() > 20)
					{
						//get the data for this bin 
						CMaskedArray<float> binT = this_month_T[bin_locs];
						CMaskedArray<float> binD = this_month_D[bin_locs];

						// find good temperatures
						varraysize good_binT = npwhere(binT.m_mask, "=", false);
						varraysize good_binD = npwhere(binD.m_mask, "=", false);
							
						//and the bad dewpoints coincident with the good temperatures
						valarray<bool> dummy = binD.m_mask[good_binT];
						varraysize bad_binD_T = npwhere(dummy, "=", false);

						if (bad_binD_T.size() != 0)
						{
							double bad_D_fraction = 1. - bad_binD_T.size() / good_binT.size();
							//if more than 50% missing
							if (bad_D_fraction >= 0.5)
							{
								//check the temporal resolution - if OK, then flag
								//This is a better temporal resolution calculation than IDL - will pick up more months
								size_t T_resoln = idl_median(npDiff(goodT));
								size_t D_resoln = idl_median(npDiff(goodD));

								//only flag if resolutions the same and number of observations in total are similar
								if ((T_resoln != D_resoln) && goodD.size()/goodT.size() < 0.666)
									continue;
								else
								{
									//flags[int(month[0]):int(month[1])][goodD] = 1
									/*
										break this loop testing bins as whole month flagged
										no point doing any more testing
										return to month loop and move to next month'''
									*/
									break;
								}
							}
						}
					}
				}
			}
			m++;
		}
		size_t nflags = npwhere(flags, "!", float(0)).size();
		print_flagged_obs_number(logfile, "Dewpoint Cut-off", "temperature", nflags);

		return flags;
	}

	void hcc(CStation &station, std::vector<int> flag_col, boost::gregorian::date start, boost::gregorian::date  end, std::ofstream& logfile)
	{
		CMetVar& temperatures = station.getMetvar("temperatures");
		CMetVar& dewpoints = station.getMetvar("dewpoints");
		std::vector<pair<int, int>> month_ranges = month_starts_in_pairs(start, end);

		//Supersaturation

		varrayfloat flags(float(0), temperatures.getData().size());
		hcc_sss(temperatures.getData(), dewpoints.getData(), month_ranges,logfile,flags);
		station.setQc_flags(flags, flag_col[0]);

		//Dew point depression
		
		CMetVar& precip = station.getMetvar("precip1_depth");
		CMetVar& cloudbase = station.getMetvar("cloud_base");
		CMetVar& past_sigwx = station.getMetvar("past_sigwx1");
		
		varrayfloat times=station.getMetvar("time").getData();
		
		flags.resize(temperatures.getAllData().size());
		flags = hcc_dpd(times, temperatures.getAllData(), dewpoints.getAllData(), precip.getAllData(), cloudbase.getAllData(), past_sigwx.getAllData(), logfile);
		station.setQc_flags(flags, flag_col[1]);
	
		//Dew point cutoffs
		flags = hcc_cutoffs(temperatures.getAllData(), dewpoints.getAllData(), month_ranges, logfile);
		station.setQc_flags(flags, flag_col[2]);
		
		for (size_t col : {0, 1, 2})
		{
			varraysize flag_locs = npwhere(station.getQc_flags(flag_col[col]), "!", float(0));
			dewpoints.setFlags(flag_locs, float(1));
		}
		append_history(station, "Temperature-Humidity Cross Check");
	}

}