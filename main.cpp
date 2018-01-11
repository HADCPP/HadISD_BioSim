

#include "station.h"
#include "utils.h"
#include "Internal_checks.h"
#include "neighbour_checks.h"

//#include "netCDFUtils.h"

#include <boost/filesystem.hpp>
#include <boost/date_time.hpp>

#include <string>
#include <stdio.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace boost::filesystem;

const string DATES[2]={ "20100101", "20110101"};



int main(int arg, char * argv)
{

	//Lire les Ids et informations concernant les stations dans le fichier 
	
	
	string file = "D:\\HadISD_BioSim\\Data\\Meteo_data\\East-Canada\\Quebec 2002-2015.HourlyHdr.csv";
	
	vector<CStation> station_info;
	//Obtenir la liste des stations et mettre les infos dans des objets CStation
	DATA_READING::read_station_metadata(file, station_info);  // utils.cpp

	bool second = false;
	bool neighbour_checks = true;

	test internal_tests;

	

	boost::gregorian::date  DATESTART = boost::gregorian::date_from_iso_string(DATES[0]);
	boost::gregorian::date  DATEEND = boost::gregorian::date_from_iso_string(DATES[1]);

	cout << boost::gregorian::day_clock::local_day() << endl;  //Date du jour au format Www Mmm dd hh:mm:ss yyyy\n


	for (CStation& station : station_info)
	{
		////////////////////////////////******  INITIALISATION DES TESTS********************//////////////////////////////////

		internal_tests.climatological = false;
		internal_tests.diurnal = false;
		internal_tests.cloud = false;
		internal_tests.duplicate = false;
		internal_tests.frequent = false;
		internal_tests.gap = false;
		internal_tests.humidity = false;
		internal_tests.odd = false;
		internal_tests.records = false;
		internal_tests.spike = false;
		internal_tests.streaks = false;
		internal_tests.variance = false;
		internal_tests.winds = false;

		///////////////////////////////////////////////////////////////////////////////////////////////

		
		if(!DATA_READING::readData(station, DATESTART, DATEEND,internal_tests))
		{
			std::cout << "CSV file not found for station " << station.getName() << endl;
			continue;
		}
		else
		{
			// Ouverture du fichier log en mode ecriture 
			std::ofstream logfile;
			stringstream sst;
			sst << LOG_OUTFILE_LOCS << (station).getId() << ".log";

			try{ logfile.open(sst.str().c_str()); }
			catch (std::exception e)
			{
				std::cout << e.what() << endl;
			}

			logfile << boost::gregorian::day_clock::local_day() << endl;
			logfile << " Internal Checks " << endl;
			logfile << " Station Identifier  :  " << (station).getId() << endl;

		

			///////////////// TESTS INTERNES //////////////////////////////////////////////////////////////////
			INTERNAL_CHECKS::internal_checks(station, internal_tests, DATESTART, DATEEND, logfile);

		}
	}	
	for (CStation& station : station_info)
	{
		///////////////// TEST EXTERNE ////////////////////////////////////////////////////////////////
		std::ofstream logfile;
		logfile.open(LOG_OUTFILE_LOCS + (station).getId() + ".log", std::ofstream::app);

		//if(neighbour_checks) NEIGHBOUR_CHECKS::neighbour_checks(station, station_info, DATESTART, DATEEND, logfile);

		//////////////// SAVE CHANGES ////////////////////////////////////////////////////////////

		DATA_READING::writeCSV(station, DATESTART, DATEEND);
	}
	
	std::system("PAUSE");
	return 0;
	
}




	