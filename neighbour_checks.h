#pragma once

#include "Utilities.h"
#include "station.h"
#include "utils.h"
#include "Internal_checks.h"

#include <vector>
#include <map>
#include <ctime>
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include <errno.h>
#include <exception>
#include <cmath>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/chrono.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace
{
	std::map<const char*, int> FLAG_OUTLIER_DICT = { { "temperatures", 41 }, { "dewpoints", 42 }, { "slp", 43 } };
}


namespace NEIGHBOUR_CHECKS
{


	void get_distances_angles(CStation station, const std::vector<CStation>& station_info,
		boost::numeric::ublas::matrix<float>& distances, boost::numeric::ublas::matrix<float>& angles);
	/*
	Get the neighbours for the station using distance, angles and elevations
	Returns all neighbours within distance and elevation ranges

	:param float st_elev: station elevation
	:param array distances: distances of all other stations
	:param array bearings: bearings of all other stations
	:param array elevations: elevations of all other stations
	:param float sep_limit: separation limit (default 500km)
	:param float elev_limit: elevation limit (default 300m)
	:param int max_neighbours: maximum number to return (default 50)

	:returns: final set of neighbours as numbers in station list
	*/
	void get_all_neighbours(int station_loc,float st_elev, fvector& distances, fvector& bearings, fvector& elevations , ivector& neighbours, ivector& neighbour_quadrants,float sep_limit = 500, float elev_limit = 300, int max_neighbours = 20);


	/* 
		Format the time series to get allow for sensible correlations

		 - Process into 24h x N_days
		 - Obtain daily average and hence hourly anomalies from daily average (removes annual cycle)
		 - Obtain hourly average of anomalies to get double-anomalies (removes diurnal cycle)

		:param array timeseries: data to be processed in 1-D array
		:param int obs_per_day: number of observations per 24hr to get a daily average to process (default = 6)

		:returns: anomalies - timeseries with removed annual and diurnal cycles.
	*/
	CMaskedArray<float> hourly_daily_anomalies(CMaskedArray<float> timeseries, int obs_per_day = 6);

	void  select_neighbours(CStation station, std::string variable, const std::vector<CStation>& station_info, ivector neighbours,
		ivector neighbour_distances, ivector neighbour_quadrants, boost::gregorian::date start, boost::gregorian::date  end, std::ofstream& logfile);

	void neighbour_checks(CStation& station, const std::vector<CStation>& station_info, boost::gregorian::date start, boost::gregorian::date  end, std::ofstream&  logfile);
}