#include "Internal_checks.h"



using namespace std;
using namespace NETCDFUTILS;
using namespace netCDF;
using namespace boost;


namespace INTERNAL_CHECKS
{
	void internal_checks(vector<CStation>  &station_info, test mytest, bool second, string *DATE)
	{
		boost::gregorian::date  DATESTART = boost::gregorian::date_from_iso_string(DATE[0]);
		boost::gregorian::date  DATEEND = boost::gregorian::date_from_iso_string(DATE[1]);
		bool first = !second;
		
		int st = 0;
		for (CStation station: station_info )
		{
			
			cout << boost::gregorian::day_clock::local_day() << endl;  //Date du jour au format Www Mmm dd hh:mm:ss yyyy\n
			cout << "Station Number :  " << st + 1 << " / " << station_info.size() << endl;
			cout << "Station Identifier:  " << (station).getId() << endl;
			// Ouverture du fichier en mode ecriture lors du premier appel de la fonction ou en mode append la seconde fois
			ofstream logfile;
			stringstream sst;
			sst << LOG_OUTFILE_LOCS << (station).getId() << ".log";

			if (first)
			{
				try{logfile.open(sst.str().c_str());}
				catch (std::exception e)
				{
					cout << e.what() << endl;
				}
			}
			else
			{
				try{logfile.open(sst.str().c_str(), ofstream::out | ofstream::app);}
				catch (std::exception e)
				{
					cout << e.what() << endl;
				}
			}
			logfile << boost::gregorian::day_clock::local_day() << endl;
			logfile << " Internal Checks " << endl;
			logfile << " Station Identifier  :  " << (station).getId() << endl;

			boost::posix_time::ptime process_start_time = boost::posix_time::second_clock::local_time();
			
			/* latitude and longitude check  */
			if (std::abs((station).getLat()) > 90.)
			{
				logfile << (station).getId() << "	Latitude check   " << DATESTART.year() << DATEEND.year() - 1 << "    Unphysical latitude " << endl;
				logfile.close();
				continue;
			}
			if (std::abs((station).getLon()) > 180.)
			{
				logfile << (station).getId() << "Longitude Check" << DATESTART.year() << DATEEND.year() - 1 << "  Unphysical longitude " << endl;
				logfile.close();
				continue;
			}
			// if running through the first time
			valarray<bool> match_to_compress;
			if (first)
			{
				// tester si le fichier existe
				string filename = NETCDF_DATA_LOCS + (station).getId() + ".nc";
				boost::filesystem::path p{ filename};
				try
				{	boost::filesystem::exists(p);}
				catch ( std::exception&  e)
				{
					cout << e.what() << endl;
				}
				//read in data
				NETCDFUTILS::read(filename,station,process_var,carry_thru_vars);
				
				//lire dans le fichier netcdf
				logfile << "Total CStation record size  " << station.getMetvar("time").getData().size() << endl;
				match_to_compress = UTILS::create_fulltimes(station,process_var, DATESTART, DATEEND, carry_thru_vars);

				//Initialiser CStation.qc_flags
				
				station.InitializeQcFlags(69,station.getMetvar("time").getData().size());

				for (string var : process_var)
				{
					CMetVar& st_var = station.getMetvar(var);
					st_var.setReportingStats(UTILS::monthly_reporting_statistics(st_var, DATESTART, DATEEND));
				}
			}
			//or if second pass through?
			else if (second)
			{
				string filename = NETCDF_DATA_LOCS + (station).getId() + "_mask.nc";
				boost::filesystem::path p{ filename };
				try
				{
					boost::filesystem::exists(p);
				}
				catch (std::exception&  e)
				{
					cout << e.what() << endl;
				}

				NETCDFUTILS::read(filename, station, process_var, carry_thru_vars);
				match_to_compress = UTILS::create_fulltimes(station, process_var, DATESTART, DATEEND, carry_thru_vars);
			}
			if (mytest.duplicate) //check on temperature ONLY
			{
				//Appel à la fonction duplicate_months de qc_tests
				vector<string> variable_list = { "temperatures" };
				
				//dmc(station,variable_list, process_var,0, DATESTART, DATEEND,logfile);
			}
			if (mytest.odd)
			{
				vector<string> variable_list = { "temperatures","dewpoints","windspeeds","slp"};
				vector<int> flag_col = { 54, 55, 56, 57 };
				//occ(station, variable_list, flag_col, logfile, second);
				//UTILS::apply_windspeed_flags_to_winddir(station);
			}
			if (mytest.frequent)
			{
				
				//fvc(station, { "temperatures","dewpoints","slp" }, { 1, 2, 3 }, DATESTART, DATEEND, logfile);
			}
			if (mytest.diurnal)
			{
				/*if (std::abs(station.getLat()) <= 60.)
					dcc(station, { "temperatures" }, process_var, { 4 }, logfile);
				else
					logfile << "   Diurnal Cycle Check not run as CStation latitude  " << station.getLat()<< "  >  60 ";*/
			}
			if (mytest.gap)
			{
					//dgc(station, { "temperatures","dewpoints","slp"}, {5,6,7 },DATESTART,DATEEND, logfile);
			}
			if (mytest.records)
			{
				//krc(station, { "temperatures", "dewpoints", "windspeeds", "slp" }, { 8, 9, 10, 11 }, logfile);
				//UTILS::apply_windspeed_flags_to_winddir(station);
			}
			if (mytest.streaks)
			{
				//rsc(station, { "temperatures", "dewpoints", "windspeeds", "slp", "winddirs" }, { { 12, 16, 20 }, { 13, 17, 21 }, { 14, 18, 22 }, { 15, 19, 23 }, { 66, 67, 68 } },
				//DATESTART, DATEEND,logfile);
				//UTILS::apply_windspeed_flags_to_winddir(station);

			}
			if (mytest.climatological)
			{
				//coc( station, { "temperatures", "dewpoints" }, { 24, 25 }, DATESTART, DATEEND,logfile);
			}
			
			if (mytest.spike)
			{
				/*sc(station, {"temperatures", "dewpoints", "slp", "windspeeds"}, {27, 28, 29, 65}, DATESTART, DATEEND, logfile, second);
				UTILS::apply_windspeed_flags_to_winddir(station);*/
			}
			if (mytest.humidity)
			{
				//hcc(station,{ 30, 31, 32}, DATESTART, DATEEND, logfile);
			}
			if (mytest.cloud)
			{
				//ccc(station, { 33,34, 35, 36, 37, 38, 39, 40 }, logfile);
			}
			if (mytest.variance)
			{
				//evc(station, { "temperatures", "dewpoints", "slp", "windspeeds" }, { 58, 59, 60, 61 }, DATESTART, DATEEND, logfile);
				//UTILS::apply_windspeed_flags_to_winddir(station); 
			}
			if (mytest.winds)
			{
				wdc(station, { 62, 63, 64 }, DATESTART, DATEEND, logfile);
				 
			}
			//Write to file
			if (first)
			{
				string filename = NETCDF_DATA_LOCS + (station).getId() + "_internal.nc";
				NETCDFUTILS::write(filename, station, process_var, carry_thru_vars, match_to_compress);
			}
			else if (second)
			{
				string filename = NETCDF_DATA_LOCS + (station).getId() + "_internal2.nc";
				NETCDFUTILS::write(filename, station, process_var, carry_thru_vars, match_to_compress);
			}
			logfile << boost::gregorian::day_clock::local_day() << endl;
			logfile << "processing took " << posix_time::second_clock::local_time() - process_start_time << "  s" << endl;
			if (logfile)
				logfile.close();//Fermeture du fichier
			cout << "Internal Checks completed" << endl;
		} //end for CStation
		
	}
}