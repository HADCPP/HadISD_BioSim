#include "Internal_checks.h"


using namespace std;
using namespace boost;


namespace INTERNAL_CHECKS
{
	void internal_checks(CStation& station, test mytest, boost::gregorian::date  DATESTART, boost::gregorian::date  DATEEND, std::ofstream& logfile)
	{
		bool second = false;

		boost::posix_time::ptime process_start_time = boost::posix_time::second_clock::local_time();
			
		/* latitude and longitude check  */
		if (std::abs((station).getLat()) > 90.)
		{
			logfile << (station).getId() << "	Latitude check   " << DATESTART.year() << DATEEND.year() - 1 << "    Unphysical latitude " << endl;
			logfile.close();
			exit(1);
		}
		if (std::abs((station).getLon()) > 180.)
		{
			logfile << (station).getId() << "Longitude Check" << DATESTART.year() << DATEEND.year() - 1 << "  Unphysical longitude " << endl;
			logfile.close();
			exit(1);
		}
		// if running through the first time
		valarray<bool> match_to_compress;
		
			
		logfile << "Total CStation record size  " << station.getMetvar("time").getData().size() << endl;
		match_to_compress = UTILS::create_fulltimes(station, process_var, DATESTART, DATEEND, carry_thru_vars);

		//Initialiser CStation.qc_flags

		station.InitializeQcFlags(69, station.getMetvar("time").getData().size());

		for (string var : process_var)
		{
			CMetVar& st_var = station.getMetvar(var);
			st_var.setReportingStats(UTILS::monthly_reporting_statistics(st_var, DATESTART, DATEEND));
		}
		
		

		if (mytest.duplicate) //check on temperature ONLY
		{
			//Appel à la fonction duplicate_months de qc_tests
			vector<string> variable_list = { "temperatures"};
			dmc(station,variable_list, process_var,0, DATESTART, DATEEND,logfile);
		}
		if (mytest.odd)
		{
			vector<string> variable_list = { "temperatures","dewpoints","windspeeds","slp"};
			vector<int> flag_col = { 54, 55, 56, 57 };
			occ(station, variable_list, flag_col, logfile, second);
			UTILS::apply_windspeed_flags_to_winddir(station);
		}
		if (mytest.frequent)
		{
			fvc(station, { "temperatures","dewpoints","slp" }, { 1, 2, 3 }, DATESTART, DATEEND, logfile);
		}
		if (mytest.diurnal)
		{
			if (std::abs(station.getLat()) <= 60.)
				dcc(station, { "temperatures" }, process_var, { 4 }, logfile);
			else
				logfile << "   Diurnal Cycle Check not run as CStation latitude  " << station.getLat()<< "  >  60 ";
		}
		if (mytest.gap)
		{
				dgc(station, { "temperatures","dewpoints","slp"}, {5,6,7 },DATESTART,DATEEND, logfile);
		}
		if (mytest.records)
		{
			krc(station, { "temperatures", "dewpoints", "windspeeds", "slp" }, { 8, 9, 10, 11 }, logfile);
			UTILS::apply_windspeed_flags_to_winddir(station);
		}
		if (mytest.streaks)
		{
			rsc(station, { "temperatures", "dewpoints", "windspeeds", "slp", "winddirs" }, { { 12, 16, 20 }, { 13, 17, 21 }, { 14, 18, 22 }, { 15, 19, 23 }, { 66, 67, 68 } },
			DATESTART, DATEEND,logfile);
			UTILS::apply_windspeed_flags_to_winddir(station);

		}
		if (mytest.climatological)
		{
			coc( station, { "temperatures", "dewpoints" }, { 24, 25 }, DATESTART, DATEEND,logfile);
		}
			
		if (mytest.spike)
		{
			sc(station, {"temperatures", "dewpoints", "slp", "windspeeds"}, {27, 28, 29, 65}, DATESTART, DATEEND, logfile, second);
			UTILS::apply_windspeed_flags_to_winddir(station);
		}
		if (mytest.humidity)
		{
			hcc(station,{ 30, 31, 32}, DATESTART, DATEEND, logfile);
		}
		if (mytest.cloud)
		{
			//ccc(station, { 33,34, 35, 36, 37, 38, 39, 40 }, logfile);
		}
		if (mytest.variance)
		{
			evc(station, { "temperatures", "dewpoints", "slp", "windspeeds" }, { 58, 59, 60, 61 }, DATESTART, DATEEND, logfile);
			UTILS::apply_windspeed_flags_to_winddir(station); 
		}
		if (mytest.winds)
		{
			wdc(station, { 62, 63, 64 }, DATESTART, DATEEND, logfile);
				 
		}
		
		logfile << boost::gregorian::day_clock::local_day() << endl;
		logfile << "processing took " << posix_time::second_clock::local_time() - process_start_time << "  s" << endl;
		if (logfile)
			logfile.close();//Fermeture du fichier
		cout << "Internal Checks completed" << endl;
	} //end for CStation
		
	
}