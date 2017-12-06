#include "clouds.h"


using namespace std;
using namespace UTILS;
using namespace boost;
using namespace PYTHON_FUNCTION;
using namespace boost::numeric::ublas;


namespace INTERNAL_CHECKS
{

	void unobservable(CStation& station, std::vector<int> flag_col, std::ofstream& logfile)
	{
		//for each cloud variable, find bad locations and flag
		int c = 0;
		for (string cloud : { "total_cloud_cover", "low_cloud_cover", "mid_cloud_cover", "high_cloud_cover" })
		{
			CMetVar& cloud_obs = station.getMetvar(cloud);
			varraysize bad_locs = npwhereOr(cloud_obs.getData(), float(9), float(10), "=", "=");

			station.setQc_flags(float(1), bad_locs, flag_col[c]);
			varraysize flag_locs = npwhere(station.getQc_flags(flag_col[c]), "=", float(1));
			print_flagged_obs_number(logfile, "Unobservable cloud", cloud, flag_locs.size());
			//copy flags into attribute
			cloud_obs.setFlags(flag_locs, float(1));
			c ++ ;
		}

	}
	void total_lt_max(CStation& station, int flag_col, std::ofstream& logfile)
	{
		CMetVar& total = station.getMetvar("total_cloud_cover");
		CMetVar& low = station.getMetvar("low_cloud_cover");
		CMetVar& mid = station.getMetvar("mid_cloud_cover");
		CMetVar& high = station.getMetvar("high_cloud_cover");

		varrayfloat maximum=low.getData();
		for (size_t i = 0; i < maximum.size(); i++)
		{
			if (maximum[i] < mid.getData()[i]) maximum[i] = mid.getData()[i];
			if (maximum[i] < high.getData()[i]) maximum[i] = high.getData()[i];
		}
	
		varraysize bad_locs = npwhere(total.getData(), "<", maximum);

		station.setQc_flags(float(1), bad_locs, flag_col);
		varraysize flag_locs = npwhere(station.getQc_flags(flag_col), "!", float(0));

		print_flagged_obs_number(logfile, "Total < Max cloud", "cloud", flag_locs.size());
		//copy flags into attribute

		total.setFlags(flag_locs, float(1));
		low.setFlags(flag_locs, float(1));
		mid.setFlags(flag_locs, float(1));
		high.setFlags(flag_locs, float(1));

	}
	void low_full(CStation& station, int flag_col, std::ofstream& logfile)
	{
		CMetVar& low = station.getMetvar("low_cloud_cover");
		CMetVar& mid = station.getMetvar("mid_cloud_cover");
		CMetVar& high = station.getMetvar("high_cloud_cover");

		varraysize low_full_locs = npwhere(low.getAllData(), '=', float(8));

		valarray<bool> dummy = mid.getAllData().m_mask[low_full_locs];
		varraysize	bad_mid = npwhere(dummy, "!", true);
		
		station.setQc_flags(float(1), low_full_locs[bad_mid], flag_col);

		dummy.resize(mid.getAllData().m_mask.size());
		dummy = high.getAllData().m_mask[low_full_locs];
		varraysize	bad_high = npwhere(dummy, "!", true);
		station.setQc_flags(float(1), low_full_locs[bad_high], flag_col);
		
		varraysize flag_locs = npwhere(station.getQc_flags(flag_col), "!", float(0));
		print_flagged_obs_number(logfile, "Low full cloud", "cloud", flag_locs.size());

		//copy flags into attribute
		mid.setFlags(flag_locs, float(1));
		high.setFlags(flag_locs, float(1));
	}

	void mid_full(CStation& station, int flag_col, std::ofstream& logfile)
	{
		CMetVar& mid = station.getMetvar("mid_cloud_cover");
		CMetVar& high = station.getMetvar("high_cloud_cover");

		varraysize mid_full_locs = npwhere(mid.getAllData(), '=', float(8));

		valarray<bool> dummy = high.getAllData().m_mask[mid_full_locs];

		varraysize bad_high = npwhere(dummy, "!", true);
		
		station.setQc_flags(float(1), mid_full_locs[bad_high], flag_col);


		varraysize flag_locs = npwhere(station.getQc_flags(flag_col), "!", float(0));
		print_flagged_obs_number(logfile, "Mid full cloud", "cloud", flag_locs.size());

		//copy flags into attribute
		high.setFlags(flag_locs, float(1));
	}

	void fix_cloud_base(CStation& station)
	{
		CMetVar& cloud_base = station.getMetvar("cloud_base");

		varraysize bad_cb = npwhere(cloud_base.getAllData(), '=', float(22000));

		//no flag set on purpose - just set to missing (unobservable)

		cloud_base.setData(bad_cb, Cast<float>(cloud_base.getMdi()),true);

	}
	void negative_cloud(CStation& station, int flag_col, std::ofstream& logfile)
	{
		//go through each cloud varaible and flag bad locations
		int c = 0;
		for (string cloud : { "total_cloud_cover", "low_cloud_cover", "mid_cloud_cover", "high_cloud_cover" })
		{

			CMetVar& cloud_obs = station.getMetvar(cloud);
			varraysize bad_locs = ma_masked_where(cloud_obs.getAllData(), "<", float(0));

			station.setQc_flags(float(1), bad_locs, flag_col);

			//copy flags into attribute
			cloud_obs.setFlags(bad_locs, float(1));
			c++;
		}
		varraysize flag_locs = npwhere(station.getQc_flags(flag_col), "!", float(0));
		print_flagged_obs_number(logfile, "Negative cloud", "cloud", flag_locs.size());

		
	}

	void ccc(CStation& station, std::vector<int> flag_col, std::ofstream& logfile)
	{
		if (flag_col.size() != 8)
		{
			cout << "insufficient flag columns given" << endl;
			return;
		}

		unobservable(station, { 33, 34, 35, 36 }, logfile);
		total_lt_max(station,flag_col[4], logfile);
		low_full(station, flag_col[5], logfile);
		mid_full(station, flag_col[6], logfile);
		fix_cloud_base(station);
		negative_cloud(station, flag_col[7], logfile);
		append_history(station, "Cloud - Logical Cross Check");
	}
}

