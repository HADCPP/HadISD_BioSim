#pragma once
#include "station.h"
#include "netCDFUtils.h"
#include "netcdf.h"
#include "ncFile.h"
#include "ncDim.h"
#include "ncVar.h"
#include "Utilities.h"
#include "utils.h"

// Tests
#include "duplicate_months.h"
#include "odd_cluster.h"
#include "frequent_values.h"
#include "Diurnal_cycle.h"
#include "records.h"
#include "streaks.h"
#include "spike.h"
#include "humidity.h"
#include "clouds.h"
#include "climatological.h"
#include "distributional_gap.h"
#include "variance.h"
#include "winds.h"

#include <vector>
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


struct test
{
	bool duplicate = false;
	bool odd = false;
	bool frequent = false;
	bool diurnal = false;
	bool gap = true;
	bool records = true;
	bool streaks = true;
	bool climatological = true;
	bool spike = true;
	bool humidity = true;
	bool cloud = true;
	bool variance = true;
	bool winds = true;
};
typedef struct test test;
namespace
{
	std::vector<std::string> process_var{ "temperatures", "dewpoints", "windspeeds", "winddirs", "slp", "total_cloud_cover", "low_cloud_cover", "mid_cloud_cover", "high_cloud_cover" };
	std::vector<std::string> carry_thru_vars{ "cloud_base", "precip1_depth", "precip1_period", "wind_gust", "past_sigwx1" };
	
}
namespace INTERNAL_CHECKS

{
	void internal_checks(std::vector<CStation> &station_info, test mytest , bool second,std::string *DATE);
}