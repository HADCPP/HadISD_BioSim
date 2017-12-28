#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime> 
#include <errno.h>
#include <exception>
#include <chrono>
#include <cmath>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/date_time.hpp>
#include "duplicate_months.h"
#include "station.h"
#include "netCDFUtils.h"
#include "Internal_checks.h"



using namespace std;
using namespace boost;
using namespace NETCDFUTILS;
using namespace INTERNAL_CHECKS;



 string DATES[2]={ "20020101", "20150101"};

//string LOG_OUTFILE_LOCS = "P: \\Project\\wbsTools\\HadISD\\Weather\\Data\\";
//string CSV_OUTFILE_LOCS = "P: \\Project\\wbsTools\\HadISD\\Weather\\PEI 2015-2016H\\";



void read_file(const string file, vector<CStation>& station_info)
{
	/*    
		Fonction pour lire le contenu du fichier Quebec 2016-2017.HourlyHdr.csv et affecter ses éléments aux attributs
		de la classe CStation
	*/
	
		char_separator<char> sep(",");
		ifstream input;
		stringstream sst;
		sst << file;
		string chaine = sst.str();
		try
		{ 
			input.open(chaine.c_str());
		}
		catch (std::exception& e)
		{
			cout << e.what() << endl;
		}
		if (!input) exit(1);
		string ligne="";
		getline(input, ligne);
		/*while (!input.eof())
		{
			getline(input, ligne);
			int i = 0;
			string data[22] = { "0" };
			tokenizer<char_separator<char>> tokens(ligne, sep);
			for (tokenizer<char_separator<char>>::const_iterator t = tokens.begin(); t != tokens.end();t++)
			{

				(t->c_str() != "") ? data[i] = t->c_str() : data[i] = "";
				i++;
			}
			
			CStation station = CStation::CStation(data[0], data[1], atof(data[2].c_str()), atof(data[3].c_str()), 
							atof(data[4].c_str()), data[12].c_str());
						station_info.push_back(station);
		}*/
		
		
		char  delimiter = ',';
		while (!input.eof())
		{
			getline(input, ligne);
			istringstream iss(ligne);
			string token;
			int i = 0;
			string data[22] = { "0" };

			while (getline(iss, token, delimiter))
			{
				data[i] = token.c_str();
				i++;
			}
			CStation stat = CStation(data[0], data[1], atof(data[2].c_str()), atof(data[3].c_str()), atof(data[4].c_str()), data[12].c_str());
			station_info.push_back(stat);

		}
}

void ncdf(const vector<CStation>& station_info,test& internal_tests,bool second)
{
		
	for (CStation station : station_info)
	{
		
		string nom_file = (station).getName() + " [" + (station).getId() + "]";
		
		string fichier = CSV_OUTFILE_LOCS + nom_file + ".csv";
		boost::filesystem::path p{ fichier };
		if (!boost::filesystem::exists(p))
		{
			string nom = (station).getName();
			auto pos = std::string::npos;
			while ((pos = nom.find('é')) != std::string::npos || (pos = nom.find('è')) != std::string::npos || (pos = nom.find('ê')) != std::string::npos)
			{
					nom.replace(pos, 1, "e");
			}
			while ((pos = nom.find('à')) != std::string::npos)
			{
					nom.replace(pos, 1, "a");
			}
			while ((pos = nom.find('\'')) != std::string::npos)
			{
					nom.replace(pos, 1 , " ");
			}

			nom_file = nom + " [" + (station).getId() + "]";
			fichier = CSV_OUTFILE_LOCS + nom_file + ".csv";
			boost::filesystem::path p1{ fichier };
			if (boost::filesystem::exists(p1))
			{
				cout << p1.filename() << endl;
				cout << "\n Making NetCDF files from CSV files \n" << endl;
				cout << "Reading data from  " << p.parent_path() << endl;
				cout << "Writing data to   " << NETCDF_DATA_LOCS << endl;
				std::size_t pos = nom_file.find(' ');
				//nom_file = nom_file.substr(0,pos);
				if (!NETCDFUTILS::MakeNetcdfFiles(fichier, DATES, station, internal_tests)) // Creer le fichier netcdf
				{
					//Exécuter les tests
					INTERNAL_CHECKS::internal_checks(station, internal_tests, second, DATES);
				}
			}
			else continue;
		}
		else
		{
			cout << p.filename() << endl;
			cout << "\n Making NetCDF files from CSV files \n" << endl;
			cout << "Reading data from  " << p.parent_path() << endl;
			cout << "Writing data to   " << NETCDF_DATA_LOCS << endl;
			std::size_t pos = nom_file.find(' ');
			//nom_file = nom_file.substr(0,pos);
			if (!NETCDFUTILS::MakeNetcdfFiles(fichier, DATES, station, internal_tests)) // Creer le fichier netcdf
			{
				//Exécuter les tests
				INTERNAL_CHECKS::internal_checks(station, internal_tests, second, DATES);
			}
		} 
		
	}
	cout << "NetCDF files created" << endl;
}



int main(int arg, char * argv)
{

	//Lire les Ids et informations concernant les stations dans le fichier 
	//
	
	string file = "D:\\HadISD_BioSim\\Data\\Meteo_data\\East-Canada\\Quebec 2002-2015.HourlyHdr.csv";
	
	vector<CStation> station_info;
	//Obtenir la liste des stations et mettre les infos dans des objets CStation
	read_file(file,station_info);
	bool second = false;
	test internal_tests;
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
	if (!second)
	{//Creer des fichiers netcdf à partir de chaque csv file.
		ncdf(station_info,internal_tests,second);
	}
	
	std::system("PAUSE");
	return 0;
	
}




	