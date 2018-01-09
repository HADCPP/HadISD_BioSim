#include "Internal_checks.h"


using namespace std;
using namespace boost;


namespace INTERNAL_CHECKS
{
	void internal_checks(CStation& station, test mytest, boost::gregorian::date  DATESTART, boost::gregorian::date  DATEEND, std::ofstream& logfile)
	{
		bool second = false;

		boost::posix_time::ptime process_start_time = boost::posix_time::second_clock::local_time();
			
		///////////////// LATITUDE AND LONGITUDE CHECK  ///////////////////
		if (std::abs((station).getLat()) > 90.)
		{
			logfile << (station).getId() << "	Latitude check   " << DATESTART.year() << DATEEND.year() - 1 << "    Unphysical latitude " << endl;
			logfile.close();
			return;
		}
		if (std::abs((station).getLon()) > 180.)
		{
			logfile << (station).getId() << "Longitude Check" << DATESTART.year() << DATEEND.year() - 1 << "  Unphysical longitude " << endl;
			logfile.close();
			return;
		}
		
		logfile << "Total CStation record size  " << station.getMetvar("time").getData().size() << endl;

		//expand the time axis of the variables ( create_fulltimes)

		//valarray<bool> match_to_compress = UTILS::create_fulltimes(station, process_var, DATESTART, DATEEND, carry_thru_vars);

		//Initialiser station.qc_flags

		station.InitializeQcFlags(69, station.getMetvar("time").getData().size());  //matrice de taille 69*nombre_obs

		// get reporting accuracies and frequencies.
		for (string var : process_var)
		{
			CMetVar& st_var = station.getMetvar(var);
			st_var.setReportingStats(UTILS::monthly_reporting_statistics(st_var, DATESTART, DATEEND));
		}
		/////////////////////////DEBUT DES TESTS INTERNES /////////////////////////////
		if (mytest.duplicate) //check on temperature ONLY
		{
			//  Test 2: duplicate_months check
			vector<string> variable_list = {"temperatures"};
			//dmc(station,variable_list, process_var,0, DATESTART, DATEEND,logfile); // flag à la ligne 0 de station.qc_flags
		}
		if (mytest.odd)
		{
			//  Test 3: odd cluster check
			vector<string> variable_list = { "temperatures","dewpoints","windspeeds","slp"};
			vector<int> flag_col = { 54, 55, 56, 57 };
			//occ(station, variable_list, flag_col, logfile, second);
			UTILS::apply_windspeed_flags_to_winddir(station);
		}
		if (mytest.frequent)
		{
			//  Test 4: frequent value check
			//fvc(station, { "temperatures","dewpoints","slp" }, { 1, 2, 3 }, DATESTART, DATEEND, logfile);
		}
		if (mytest.diurnal)
		{
			//  Test 5: diurnal cycle check
			if (std::abs(station.getLat()) <= 60.)
				second=false;//dcc(station, { "temperatures" }, process_var, { 4 }, logfile);
			else
				logfile << "   Diurnal Cycle Check not run as CStation latitude  " << station.getLat()<< "  >  60 ";
		}
		if (mytest.gap)
		{
			//  Test 6: distributional gap check
				//dgc(station, { "temperatures","dewpoints","slp"}, {5,6,7 },DATESTART,DATEEND, logfile);
		}
		if (mytest.records)
		{
			//  Test 7: known records check
			//krc(station, { "temperatures", "dewpoints", "windspeeds", "slp" }, { 8, 9, 10, 11 }, logfile);
			UTILS::apply_windspeed_flags_to_winddir(station);
		}
		if (mytest.streaks)
		{
			//  Test 8: repeated streaks/unusual spell frequency
			rsc(station, { "temperatures", "dewpoints", "windspeeds", "slp", "winddirs" }, { { 12, 16, 20 }, { 13, 17, 21 }, { 14, 18, 22 }, { 15, 19, 23 }, { 66, 67, 68 } },
			DATESTART, DATEEND,logfile);
			UTILS::apply_windspeed_flags_to_winddir(station);

		}
		if (mytest.climatological)
		{
			//  Test 9: climatological outlier check
			//Il faut effectuer ce test sur trois ans de données au moins
			//coc( station, { "temperatures", "dewpoints" }, { 24, 25 }, DATESTART, DATEEND,logfile);
		}
			
		if (mytest.spike)
		{
			//  Test 10: spike check
			sc(station, {"temperatures", "dewpoints", "slp", "windspeeds"}, {27, 28, 29, 65}, DATESTART, DATEEND, logfile, second);
			UTILS::apply_windspeed_flags_to_winddir(station);
		}
		if (mytest.humidity)
		{
			//  Test 11: temperature and dewpoint temperature cross - check
			
			//hcc(station,{ 30, 31, 32}, DATESTART, DATEEND, logfile);
		}
		if (mytest.cloud)
		{
			//ccc(station, { 33,34, 35, 36, 37, 38, 39, 40 }, logfile);
		}
		if (mytest.variance)
		{
			//  Test 13: unusual variance check
			evc(station, { "temperatures", "dewpoints", "slp", "windspeeds" }, { 58, 59, 60, 61 }, DATESTART, DATEEND, logfile);
			UTILS::apply_windspeed_flags_to_winddir(station); 
		}
		if (mytest.winds)
		{
			//wdc(station, { 62, 63, 64 }, DATESTART, DATEEND, logfile);
				 
		}
		
		logfile << boost::gregorian::day_clock::local_day() << endl;
		logfile << "processing took " << posix_time::second_clock::local_time() - process_start_time << "  s" << endl;
		if (logfile)
			logfile.close();
		cout << "Internal Checks completed" << endl;
	} //end for CStation
		
	
}