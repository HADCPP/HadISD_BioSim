

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


const string DATES[2]={ "20020101", "20150101"};


int main(int arg, char * argv)
{

	//Lire les Ids et informations concernant les stations dans le fichier 
	
	
	string file = "D:\\HadISD_BioSim\\Data\\Meteo_data\\East-Canada\\Quebec 2002-2015.HourlyHdr.csv";
	
	vector<CStation> station_info;
	//Obtenir la liste des stations et mettre les infos dans des objets CStation
	DATA_READING::read_station_metadata(file, station_info);  // utils.cpp

	bool second = false;

	test internal_tests;

	////////////////////////////////******  INITIALISATION DES TESTS********************//////////////////////////////////

	internal_tests.climatological = true;
	internal_tests.diurnal = true;
	internal_tests.cloud = false;
	internal_tests.duplicate = true;
	internal_tests.frequent = false;
	internal_tests.gap = true;
	internal_tests.humidity = true;
	internal_tests.odd = true;
	internal_tests.records = true;
	internal_tests.spike = true;
	internal_tests.streaks = true;
	internal_tests.variance = true;
	internal_tests.winds = true;

	///////////////////////////////////////////////////////////////////////////////////////////////

	boost::gregorian::date  DATESTART = boost::gregorian::date_from_iso_string(DATES[0]);
	boost::gregorian::date  DATEEND = boost::gregorian::date_from_iso_string(DATES[1]);

	cout << boost::gregorian::day_clock::local_day() << endl;  //Date du jour au format Www Mmm dd hh:mm:ss yyyy\n


	for (CStation station : station_info)
	{

		if(!DATA_READING::readData(station, DATESTART, DATEEND))
		{
			std::cout << "CSV file not found for station " << station.getName() << endl;
			continue;
		}
		else
		{
			// Ouverture du fichier log en mode ecriture 
			ofstream logfile;
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

			///////////////// TEST EXTERNE ////////////////////////////////////////////////////////////////
			NEIGHBOUR_CHECKS::neighbour_checks(station, station_info, DATESTART, DATEEND, logfile);
		}


	}	
	std::system("PAUSE");
	return 0;
	
}




	