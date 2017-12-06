#include "records.h"



using namespace std;
using namespace UTILS;
using namespace PYTHON_FUNCTION;

namespace INTERNAL_CHECKS
{


	string krc_get_wmo_region(const string& wmo)
	{

		
		int wmoid = atoi(wmo.c_str());
		string region = "ROW";
		if(600000 <= wmoid && wmoid <= 699999)
			region = "Africa";
		if (200000 <= wmoid && wmoid <= 200999)
			region = "Asia";
		if (202000 <= wmoid && wmoid <= 219999)
			region = "Asia";
		if (230000 <= wmoid && wmoid <= 259999)
			region = "Asia";
		if (280000 <= wmoid && wmoid <= 329999)
			region = "Asia";
		if (350000 <= wmoid && wmoid <= 369999)
			region = "Asia";
		if (380000 <= wmoid && wmoid <= 399999)
			region = "Asia";
		if (403500 <= wmoid && wmoid <= 485999)
			region = "Asia";
		if (488000 <= wmoid && wmoid <= 499999)
			region = "Asia";
		if (500000 <= wmoid && wmoid <= 599999)
			region = "Asia";
		if (800000 <= wmoid && wmoid <= 889999)
			region = "South_America";
		if (700000 <= wmoid && wmoid <= 799999)
			region = "North_America";
		if (486000 <= wmoid && wmoid <= 487999)
			region = "Pacific";
		if (900000 <= wmoid && wmoid <= 989999)
			region = "Pacific";
		if (wmoid <= 199999)
			region = "Europe";
		if (201000 <= wmoid && wmoid <= 201999)
			region = "Europe";
		if (220000 <= wmoid && wmoid <= 229999)
			region = "Europe";
		if (260000 <= wmoid && wmoid <= 279999)
			region = "Europe";
		if (330000 <= wmoid && wmoid <= 349999)
			region = "Europe";
		if (370000 <= wmoid && wmoid <= 379999)
			region = "Europe";
		if (400000 <= wmoid && wmoid <= 403499)
			region = "Europe";
		if (890000 <= wmoid && wmoid <= 899999)
			region = "Antarctica";

		return region;
}
	
	void krc_set_flags(varraysize& locs, CStation& station, int col)
	{
		if (locs.size() > 0)
			station.setQc_flags(float(1), locs, col);
	}

	void krc(CStation& station, std::vector<std::string> var_list, std::vector<int>flag_col,std::ofstream &logfile)
	{
		//Mettre à jour les records mondiaux
		std::map< std::string, float> T_X, T_N, D_X, D_N, W_X, W_N, S_X, S_N;
		string reg[8] = { "Africa", "Asia", "South_America", "North_America", "Europe", "Pacific", "Antarctica", "ROW" };
		float record[8] =  { 55.0, 53.9, 48.9, 56.7, 48.0, 50.7, 15.0, 56.7 };
		float record1[8] = { -23.9, -67.8, -32.8, -63.0, -58.1, -23.0, -89.2, -89.2 };
		float record2[8] = { 55.0, 53.9, 48.9, 56.7, 48.0, 50.7, 15.0, 56.7 };
		float record3[8] = { -50., -100., -60., -100., -100., -50., -100., -100. };
		float record4[8] = { 113.2, 113.2, 113.2, 113.2, 113.2, 113.2, 113.2, 113.2 };
		float record5[8] = { 0., 0., 0., 0.,  0., 0., 0., 0. };
		float record6[8] = { 1083.3, 1083.3, 1083.3, 1083.3, 1083.3, 1083.3, 1083.3, 1083.3 };
		float record7[8] = { 870., 870., 870., 870., 870., 870., 870., 870. };

		for (int i = 0; i < 8; i++)
		{
			T_X[reg[i]] = record[i];
			T_N[reg[i]] = record1[i];
			D_X[reg[i]] = record2[i];
			D_N[reg[i]] = record3[i];
			W_X[reg[i]] = record4[i];
			W_N[reg[i]] = record5[i];
			S_X[reg[i]] = record6[i];
			S_N[reg[i]] = record7[i];

		}

		typedef map< string, map< string, float>> Tmap;
		Tmap maxes, mins;
		maxes["temperatures"] = T_X;
		maxes["dewpoints"] = D_X;
		maxes["windspeeds"] = W_X;
		maxes["slp"] = S_X;
		mins["temperatures"] = T_N;
		mins["dewpoints"] = D_N;
		mins["windspeeds"] = W_N;
		mins["slp"] = S_N;
		int v = 0;
		for (string variable : var_list)
		{
			CMetVar& st_var = station.getMetvar(variable);
			string st_region = krc_get_wmo_region(station.getWmoId());
			CMaskedArray<float> all_filtered = UTILS::apply_filter_flags(st_var);
			varraysize too_high = npwhere<float>(all_filtered.m_data, ">", maxes[variable][st_region]);
			krc_set_flags(too_high, station, flag_col[v]);
			//make sure that don't flag the missing values!
			varraysize too_low = npwhere(all_filtered, '<', mins[variable][st_region]);
			krc_set_flags(too_low, station, flag_col[v]);
			varraysize flag_locs = npwhere<float>(station.getQc_flags()[flag_col[v]], "!", float(0));
			UTILS::print_flagged_obs_number(logfile, "World Record", variable, flag_locs.size());
			st_var.setFlags(flag_locs, float(1));
			v++;
		}
		UTILS::append_history(station, "World Record Check");
	}
}