#pragma once

#include "python_function.h"
#include "Utilities.h"
#include "station.h"
#include "utils.h"

#include "netcdf.h"
#include "ncFile.h"
#include "ncDim.h"
#include "ncVar.h"
#include "ncException.h"
#include "ncVarAtt.h"

#include <vector>
#include <map>
#include <string>
#include <typeinfo.h>
#include <valarray>




namespace
{
	std::string LOG_OUTFILE_LOCS = "D:\\HadISD++\\Data\\Log\\";
	std::string NETCDF_DATA_LOCS = "D:\\HadISD++\\Data\\NetCDF_files\\";
	std::string CSV_OUTFILE_LOCS = "D:\\HadISD++\\Data\\PEI 2015-2016H\\";
	const int INTMDI = -999;
	const float FLTMDI = -1e30;
	const int NBVAR=69; // recuperer la taille des données dans le fichier netCDF
	
}
namespace NETCDFUTILS
{
	
	/*Write coordinates as variables

	: param str outfile : output netcdf file
	: param str short_name : netcdf short_name
	: param str standard_name : netcdf standard_name
	: param str long_name : netcdf long_name
	: param str units : netcdf units
	: param str axis : netcdf axis
	: param flt data : coordinate
	*/
	template<typename T>
	inline void  write_coordinates(netCDF::NcFile *ncFile, const std::string short_name, netCDF::NcDim dim, const std::string  standard_name, const std::string  long_name, const  std::string  units, const std::string  axis, double * data)
	{
		NcVar nc_var = (*ncFile).addVar(short_name, ncFloat, dim);
		nc_var.setCompression(false, true, 9);
		nc_var.putAtt("Standard_name", standard_name);
		nc_var.putAtt("long_name", long_name);
		nc_var.putAtt("units", units);
		nc_var.putAtt("axis", axis);
		if (short_name == "alt")
			nc_var.putAtt("positive ", "up");
		nc_var.putVar(data);

	}
	/*
	Write given attributes into ncdf variable

	:param object variable: netcdf Variable
	:param string long_name: long_name value for variable to be written
	:param string cell_methods: cell_methods value for variable to be written
	:param float/int missing_value: missing_value value for variable to be written
	:param string units: units value for variable to be written
	:param string axis: axis value for variable to be written
	:param float/int vmin: valid_min value for variable to be written
	:param float/int vmax: valid_max value for variable to be written
	:param string standard_name: standard_name value for variable to be written
	:param string coordinates: coordinates to associate to variable
	"""
	*/

	template<typename T>
	inline void  write_attributes(netCDF::NcVar nc_var, const std::string long_name, const std::string cell_method, T missing_value, 
		const  std::string  units, const std::string  axis, T vmin, T vmax, const std::string coordinates, const std::string standard_name = "")
	{
		nc_var.putAtt("long_name", long_name);
		nc_var.putAtt("cell_methods", cell_method);
		nc_var.putAtt("missing_value", to_string(missing_value));
		nc_var.putAtt("units", units);
		nc_var.putAtt("valid_min", to_string(vmin));
		nc_var.putAtt("valid_max", to_string(vmax));
		nc_var.putAtt("coordinates", coordinates);
		if (standard_name != "")
			nc_var.putAtt("standard_name", standard_name);

	}
	
	template<typename T>
	void fill(std::vector<T> *tab, int tab_size)
	{
		if (typeid(T).name() == typeid(int).name())
		{
			for (int i = 0; i < tab_size; i++)
				(*tab).push_back(INTMDI);
		}
		else
		{
			for(int i = 0; i < tab_size; i++)
				(*tab).push_back(FLTMDI);
		}

	}
	/* Write flag of attributes in netCDf files
	*/
	template<typename T>
	inline void WriteFlagAttributes(netCDF::NcVar nc_var, const std::string long_name, T missing_value, const std::string  axis)
	{
		nc_var.putAtt("long_name", long_name);
		nc_var.putAtt("missing_value", boost::lexical_cast<std::string>(missing_value));
		nc_var.putAtt("units", " 1");
		//nc_var.putAtt("axis", axis);
	}
	
	template<typename T>
	inline void CreateTab(std::vector<T>& v_tab, T *tab)
	{
		
		for (int i = 0; i < v_tab.size(); i++)
		{
			*(tab+i) =  v_tab[i];
		}
	}
	template<typename T>
	inline void CreateMat(std::valarray<std::valarray<T>>& v_tab, T **tab)
	{
		
		for (int i = 0; i < v_tab.size(); i++)
		{
			for (int j = 0; j < v_tab[0].size(); j++)
			{
				tab [i][j] = v_tab[i][j];
				
			}
		}
	}
	/*  Find the attribute whose name is in "chaine" 
	  of the netCDF variable var 
	*/
	template<typename T>
	inline T getAttribute(netCDF::NcVar var, const std::string chaine)
	{
		netCDF::NcVarAtt var_att = var.getAtt(chaine);
		T att;
		
	    var_att.getValues(att);
		return att;
		
	}
	template<>
	inline float getAttribute(netCDF::NcVar var, const std::string chaine)
	{
		netCDF::NcVarAtt var_att = var.getAtt(chaine);

			float *att = new float[4];
			var_att.getValues(att);
			return *att;

	}
	void MakeNetcdfFiles(const std::string fichier, std::string *date, CStation &stat);
	void read(std::string filename, CStation& stat, std::vector<std::string>& process_var, std::vector<std::string>&  opt_var_list ,
		bool read_input_station_id = true, bool read_qc_flags = true, bool read_flagged_obs = true);
	void write(const std::string filename1, CStation & stat, std::vector<std::string> var_list,
		std::vector<std::string>  opt_var_list, std::valarray<bool> compressed, bool write_QC_flags = true,
		bool write_flagged_obs = true);
}