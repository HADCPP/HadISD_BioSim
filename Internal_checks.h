#pragma once

#include "station.h"
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



namespace INTERNAL_CHECKS

{
	void internal_checks(CStation& station, test mytest, boost::gregorian::date  DATESTART, boost::gregorian::date  DATEEND, std::ofstream& logfile);
}