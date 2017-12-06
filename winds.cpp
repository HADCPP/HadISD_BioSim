#include "winds.h"


using namespace std;
using namespace UTILS;
using namespace boost;
using namespace PYTHON_FUNCTION;
using namespace boost::numeric::ublas;

namespace INTERNAL_CHECKS
{
	void logical_checks(CStation &station, std::vector<int> flag_col, std::ofstream& logfile)
	{
		CMetVar& speed = station.getMetvar("windspeeds");
		CMetVar& direction = station.getMetvar("winddirs");

		//recover direction information where the speed is Zero

		varraysize fix_zero_direction = npwhereAnd<float,bool>(speed.getData(), "=", float(0), direction.getAllData().m_mask, "=", true);
		direction.getAllData().m_data[fix_zero_direction] = 0;
		direction.getAllData().m_mask[fix_zero_direction] = false;
		station.setQc_flags(float(-1),fix_zero_direction, flag_col[1]);

		direction.getAllData().masked();
		//negative speeds

		varraysize  negative_speed = ma_masked_where(speed.getAllData(), "<", float(0));
		station.setQc_flags(float(1), negative_speed, flag_col[0]);

		//negative directions (don't try to adjust)

		varraysize  negative_direction = ma_masked_where(direction.getAllData(), "<", float(0));
		station.setQc_flags(float(1), negative_direction, flag_col[1]);

		//  wrapped directions (don't try to adjust)
		varraysize  wrapped_direction = ma_masked_where(direction.getAllData(), ">", float(360));
		station.setQc_flags(float(1), wrapped_direction, flag_col[1]);

		//no direction possible if speed == 0
		varraysize  bad_direction = npwhereAnd(speed.getData(), "=", float(0), direction.getData(), "!", float(0));
		station.setQc_flags(float(1), bad_direction, flag_col[1]);

		//northerlies given as 360, not 0 --> calm

		varraysize  bad_speed = npwhereAnd(direction.getData(), "=", float(0), speed.getData(), "!", float(0));
		station.setQc_flags(float(1), bad_speed, flag_col[0]);

		//and output to file/screen
		varraysize flag_locs0 = npwhere(station.getQc_flags(flag_col[0]), ">", float(0));   // in case of direction fixes
		varraysize flag_locs1 = npwhere(station.getQc_flags(flag_col[1]), ">", float(0));  // in case of direction fixes

		print_flagged_obs_number(logfile, "Wind Logical Checks", "windspeeds", flag_locs0.size());
		print_flagged_obs_number(logfile, "Wind Logical Checks", "winddirs", flag_locs1.size());

		//copy flags into attribute
		speed.setFlags(flag_locs0, float(1));
		direction.setFlags(flag_locs1, float(1));
	}
	void wind_create_bins(CMaskedArray<float>&  indata, varrayfloat& binEdges, varrayfloat& bincenters)
	{
		//set up the bins
		float bmins = ma_min(indata);
		float bmax = ma_max(indata);
		float binwidth = (bmax - bmins)/20;
		//get power of ten
		int decimal_places = to_string(int(1 / binwidth)).size();
		binwidth = Round(binwidth, decimal_places);
		binEdges = arange(bmax+3*binwidth, float(0), binwidth);
		bincenters = binEdges[std::slice(1, binEdges.size() - 1, 1)];
		bincenters += binEdges[std::slice(0, binEdges.size() - 1, 1)];
		bincenters *= float(0.5);

	}

	float get_histogram_norm(varrayfloat&  indata, varrayfloat& binEdges)
	{
		varrayfloat natural_hist = histogram(indata, binEdges);
		varrayfloat normed_hist = histogram(indata, binEdges,true);

		size_t maxloc = np_argmax(natural_hist);
		return float(natural_hist[maxloc] / normed_hist[maxloc]);
	}

	void wind_rose_check(CStation &station, int flag_col, boost::gregorian::date start, boost::gregorian::date  end,std::ofstream& logfile)
	{
		CMaskedArray<float> speed = station.getMetvar("windspeeds").getAllData();
		CMaskedArray<float> direction = station.getMetvar("winddirs").getAllData();
		varrayfloat flags = station.getQc_flags()[flag_col];
		std::vector<pair<int, int>> month_ranges = month_starts_in_pairs(start, end);
		std::vector<std::valarray<pair<int, int>>> month_ranges_years= L_reshape3(month_ranges,12);
		
		//histogram of wind directions ( ~ unravelled wind-rose)
		int bw = 20;
		varrayfloat binEdges = arange(float(360 + bw), float(0), float(bw));
		varrayfloat full_hist = histogram(direction.compressed(), binEdges,true);

		//use rmse as this is known(Chi - sq remains just in case)
		CMaskedArray<float> rmse = CMaskedArray<float>(float(- 1), month_ranges_years[0].size());
		CMaskedArray<float> chisq = CMaskedArray<float>(float(-1), month_ranges_years[0].size());

		//run through each year to extract RMSE's
		varrayfloat hist(binEdges.size());
		
		for (int y = 0; y < rmse.size();y++)
		{
			size_t taille = month_ranges_years[11][y].first - month_ranges_years[0][y].first;

			varraysize indices(taille);
			for (size_t i = 0; i < indices.size();i++)
			{
				indices[i] = month_ranges_years[0][y].first+i;
			}
			CMaskedArray<float> _direction = direction[indices];

			if (_direction.compressed().size() > 0)
			{
				hist = histogram(_direction.compressed(), binEdges, true);
				chisq.m_data[y] = ((full_hist - hist)*(full_hist - hist) / (full_hist + hist)).sum() / 2.;
				rmse.m_data[y] = std::sqrt(((full_hist - hist)*(full_hist - hist)).sum() / hist.size());
			}
			else rmse.m_mask[y] = true;
			
		}
		//Now to bin up the differences and see what the fit is.
		//need to have values spread so can bin!
		valarray<float> rmse_unique = rmse.compressed()[rmse.compressed() != float(-1)];  // array contenant les éléments != de -1
		if (rmse_unique.size() > 0)
		{
			varrayfloat bincenters(binEdges.size());
			wind_create_bins(rmse, binEdges, bincenters);
			hist = histogram(rmse.m_data, binEdges);

			float norm = get_histogram_norm(rmse.m_data, binEdges);

			//inputs for fit
			double mu = rmse.m_data.sum() / rmse.size();
			double std = stdev(rmse.m_data, mu);

			//try to get decent fit to bulk of obs.

			varrayfloat gaussian = fit_gaussian(bincenters, hist, hist.max(), mu, std);
			//invert Gaussian to find initial threshold, then hunt for first gap beyond
			float threshold = invert_gaussian(PROB_THRESHOLD, gaussian);
			/*
			 if dist_pdf[-1] < PROB_THRESHOLD:
            # then curve has dropped below the threshold, so can find some updated ones.
            threshold = -np.where(dist_pdf[::-1] > PROB_THRESHOLD)[0][0]
			 else:
            threshold = bincenters[-1]*/

			int n = 0;
			int center = np_argmax(hist);
			float gap = bincenters[bincenters.size() - 1];  //nothing should be beyond this
			while (true)
			{
				if (center + n + 1 == bincenters.size())  //gone beyond edge - nothing to flag, so just break
					break;
				if (bincenters[center + n] < threshold)
				{
					n += 1;
					continue; // continue moving outwards
				}
				if (hist[center + n] == float(0))
				{
					//found one
					if (center + n + 1 == bincenters.size()) //gone beyond edge - nothing to flag, so just break
						break;
					else if (hist[center + n + 1] == float(0))
					{
						// has to be two bins wide ?
						gap = bincenters[center + n];
						break;
					}
				}
				n += 1;
			}
			//run through each year to extract RMSE's
			int y = 0;
			for (int y = 0; y < rmse.size(); y++)
			{
				if (rmse[y] > gap)
				{
					//only flag where there are observations
					size_t taille = month_ranges_years[y][0].second - month_ranges_years[y][0].first;
					varraysize indices(taille);
					for (size_t i = 0; i < indices.size(); i++)
					{
						indices[i] = month_ranges_years[y][0].first + i;
					}
					valarray<bool> dummy1 = direction.m_mask[indices];
					valarray<bool> dummy2 = speed.m_mask[indices];
					varraysize good = npwhereOr(dummy1, "=", false, dummy2, "=", false);

					varrayfloat flags_dum(float(0),indices.size());
					flags_dum[good] = 1;
					
					flags[indices] = flags_dum;
				}
				else rmse.m_mask[y] = false;
				y++;
			}
		}
		//and apply the flags and output text
		
		varraysize flag_locs = npwhere(flags, "!", float(0));
		print_flagged_obs_number(logfile, "Wind Rose Check"," windspeeds / dirs", flag_locs.size());

		station.setQc_flags(flags,flag_col);
		//and flag the variables
		
		station.getMetvar("windspeeds").setFlags(flag_locs,float(1));
		station.getMetvar("winddirs").setFlags(flag_locs, float(1));
	}
		
	void wdc(CStation &station, std::vector<int> flag_col, boost::gregorian::date start, boost::gregorian::date  end,
		std::ofstream& logfile)
	{
		//what to do about synergistic flagging - can there be more speed obs than dir obs or vv?

		cout << "running wind checks" << endl;
		logical_checks(station, { flag_col[0], flag_col[1] }, logfile);

		wind_rose_check(station, flag_col[2], start, end, logfile);
	}
}