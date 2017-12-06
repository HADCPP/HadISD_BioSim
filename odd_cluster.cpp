#include "odd_cluster.h"



using namespace std;
using namespace UTILS;
using namespace PYTHON_FUNCTION;



COddCluster::COddCluster(float start, float end, int length,float location, float data_mdi, int last_data)
{
	COddCluster::m_start = start;
	COddCluster::m_end = end;
	COddCluster::m_length = length;
	COddCluster::m_data_mdi = data_mdi;
	COddCluster::m_last_data = last_data;
	COddCluster::m_locations.push_back(location);
}
COddCluster::~COddCluster()
{
}

string COddCluster::toString()
{
	return " Odd cluster, starting :" + to_string(m_start) + ", ending " + to_string(m_end) + ", length"+ to_string(m_length);
}

namespace INTERNAL_CHECKS
{
	/*
	Just a normal observation. (Not all inputs used here,
	but required for consistency)

	:param OddCluster cluster: cluster object
	:param int obs_type: type determinator of observation
	:param int time: timestamp
	:param array flags: flag array

	:returns:
	cluster - updated cluster information
	obs_type - updated observation type
	*/

	void occ_normal(COddCluster& cluster, int &obs_type, int time, const varrayfloat& flags)
	{
		cluster.m_last_data=time;
	}
	/*currently in a potential cluster.  (Not all inputs used here,
		but required for consistency)

		: param OddCluster cluster : cluster object
		: param int obs_type : type determinator of observation
		: param int time : timestamp
		: param array flags : flag array

		: returns :
		cluster - updated cluster information
		obs_type - updated observation type*/
	void occ_in_cluster(COddCluster& cluster, int &obs_type, int time, const varrayfloat& flags)
	{
		if (cluster.m_length == 6 || time - cluster.m_start > 24)
		{
			cluster.m_start = cluster.m_data_mdi;
			cluster.m_end = cluster.m_data_mdi;
			cluster.m_length = 0;
			cluster.m_locations.push_back(cluster.m_data_mdi);
			cluster.m_last_data = cluster.m_data_mdi;
			obs_type = 0;
		}
		else
		{ //in a cluster, less than 6hr and not over 24hr - increment
			cluster.m_length += 1;
			cluster.m_locations.push_back(time);
			cluster.m_end = time;
		}
	}

	/*
	There has been a gap in the data, check if long enough to start cluster
	(Not all inputs used here, but required for consistency)

	:param OddCluster cluster: cluster object
	:param int obs_type: type determinator of observation
	:param int time: timestamp
	:param array flags: flag array

	:returns:
	cluster - updated cluster information
	obs_type - updated observation type
	*/
	void occ_start_cluster(COddCluster& cluster, int &obs_type, int time, const varrayfloat& flags)
	{
		if (time - cluster.m_last_data >= 48)
		{
			//If gap of 48hr, then start cluster increments
			obs_type = 2;
			cluster.m_length += 1;
			cluster.m_start = time;
			cluster.m_end = time;
			cluster.m_locations.push_back(time);
		}
		else  //Gap in data not sufficiently large 
		{
			obs_type = 1;
			cluster.m_last_data = time;
			cluster.m_length = 0;
			cluster.m_start = cluster.m_data_mdi;
			cluster.m_end = cluster.m_data_mdi;
			cluster.m_locations.push_back(0);
		}
	}
	/*There has been a gap in the data after a cluster;
      check if long enough to mean cluster is sufficiently isolated.
      If so, flag else reset

    :param OddCluster cluster: cluster object
    :param int obs_type: type determinator of observation
    :param int time: timestamp
    :param array flags: flag array 
    
    :returns:
       cluster - updated cluster information
       obs_type - updated observation type
	*/
	void occ_after_cluster(COddCluster& cluster, int &obs_type, int time, const varrayfloat& flags)
	{
		if (time - cluster.m_end >= 48)
		{
			//isolated cluster with 48hr gap either side
			//flags[cluster.locations] = 1
			//as have had a 48hr gap, start a new cluster
			cluster.m_last_data = cluster.m_end;
			cluster.m_start = time;
			cluster.m_end = time;
			cluster.m_locations = { time };
			cluster.m_length = 1;
			obs_type = 2;
		}
		else if ((time - cluster.m_start) <= 24 && cluster.m_length < 6)
		{
			// have data, but cluster is small and within thresholds--> increment
			obs_type = 2;
			cluster.m_length += 1;
			cluster.m_locations.push_back(time);
			cluster.m_end = time;
		}
		else
		{
			//actually it is now normal data, so reset
			obs_type = 0;
			cluster.m_last_data = time;
			cluster.m_start = cluster.m_data_mdi;
			cluster.m_end = cluster.m_data_mdi;
			cluster.m_length = 0;
			cluster.m_locations.push_back(0);
		}
	}

	void occ(CStation& station, std::vector<std::string>& variable_list, std::vector<int>& flag_col,
		std::ofstream &logfile, bool second)
	{
		
		
		int v = 0;
		for (string variable : variable_list)
		{
			CMetVar& st_var = station.getMetvar(variable);
			CMaskedArray<float> filtered_data = UTILS::apply_filter_flags(st_var);

			varrayfloat var_flags = station.getQc_flags()[v];	
			unsigned int prev_flag_number = 0;
			if (second)
			{
				//count currently existing flags
				varrayfloat dummy = var_flags[var_flags != float(0)];
				prev_flag_number = dummy.size();
			}
			//using IDL copy as method to ensure reproducibility (initially)
			COddCluster oc_details = COddCluster::COddCluster(UTILS::Cast<float>(st_var.getMdi()), UTILS::Cast<float>(st_var.getMdi()), 0, UTILS::Cast<float>(st_var.getMdi()), UTILS::Cast<float>(st_var.getMdi()), -1);
			int obs_type = 1;
			for (int time : station.getMetvar("time").getData())
			{
				if (filtered_data.m_mask[time] == false)
				{
					//process observation point 
					switch (obs_type)
					{
					case 0: occ_normal(oc_details, obs_type, time, var_flags);
						break;
					case 1: occ_start_cluster(oc_details, obs_type, time, var_flags);
						break;
					case 2:occ_in_cluster(oc_details, obs_type, time, var_flags);
						break;
					case 3:occ_after_cluster(oc_details, obs_type, time, var_flags);
						break;
					}
				}
				else
				{
					//have missing data
					if (obs_type == 2)
						obs_type = 3;
					else if (obs_type == 0)
						obs_type = 1;
				}
				
			}
			station.setQc_flags(var_flags, v);
			varraysize flag_locs = PYTHON_FUNCTION::npwhere<float>(station.getQc_flags()[flag_col[v]], "!", float(0));
			UTILS::print_flagged_obs_number(logfile, "Odd Cluster", variable, flag_locs.size() - prev_flag_number);
			//copy flags into attribute
			st_var.setFlags(flag_locs, float(1));
			v++;
		}
		UTILS::append_history(station, "Isolated Odd Cluster Check");
	}
}