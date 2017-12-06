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



 string DATES[2]={ "20060101", "20100101"};

//string LOG_OUTFILE_LOCS = "P: \\Project\\wbsTools\\HadISD\\Weather\\Data\\";
//string CSV_OUTFILE_LOCS = "P: \\Project\\wbsTools\\HadISD\\Weather\\PEI 2015-2016H\\";



void read_file(string file, vector<CStation> &station_info)
{
	/*    
		Fonction pour lire le contenu du fichier candidate_stations.txt et affecter ses éléments aux attributs
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
		catch (std::exception e)
		{
			cout << e.what() << endl;
		}
		if (!input) exit(1);
		string ligne="";
		getline(input, ligne);
		while (!input.eof())
		{
			getline(input, ligne);
			int i = 0;
			string data[12] = { "0" };
			tokenizer<char_separator<char>> tokens(ligne, sep);
			for (const string& t : tokens)
			{
				data[i] = t.c_str();
				//cout << data[i] << endl;
				i++;
			}
			
			CStation station = CStation::CStation(data[0], data[1], atof(data[2].c_str()), atof(data[3].c_str()), 
							atof(data[4].c_str()), data[7].c_str());
						station_info.push_back(station);
		}
		/*for (size_t i = 0, len = station_info.size(); i < len; i++)
			cout <<station_info.at(i).toString()<<endl;*/
		
		/*char  delimiter = ';';
		while (!input.eof())
		{
			getline(input, ligne);
			cout << ligne << endl;
			istringstream iss(ligne);
			string token;
			int i = 0;
			string data[12] = { "0" };

			while (getline(iss, token, delimiter))
			{
				data[i] = token.c_str();
				i++;
			}
			CStation stat = CStation::CStation(data[0], atof(data[2].c_str()), atof(data[3].c_str()), atof(data[4].c_str()));
			station_info.push_back(stat);*/
		
		
}

void ncdf(vector<CStation> station_info)
{
		
	for (CStation stat : station_info)
	{
		
		string nom_file = (stat).getName() + " [" + (stat).getId() + "]";
		string fichier = CSV_OUTFILE_LOCS + nom_file + ".csv";
		boost::filesystem::path p{ fichier };
		boost::filesystem::exists(p);
		cout << p.filename() << endl;
		cout << "\n Making NetCDF files from CSV files \n" << endl; 
		cout << "Reading data from  " << p.parent_path() << endl;
		cout <<"Writing data to   "<<NETCDF_DATA_LOCS << endl;
		std::size_t pos = nom_file.find(' ');
		//nom_file = nom_file.substr(0,pos);
		MakeNetcdfFiles(fichier, DATES, stat); // Creer les fichiers netcdf
	}
	cout << "NetCDF files created" << endl;
}



int main(int arg, char * argv)
{


	//faire un changement
	//fichier d'entrée
	
	//char* file = "C:\\Users\\oali\\Documents\\Visual Studio 2013\\Projects\\duplicate_months\\duplicate_months\\candidate_stations.txt";
	string file = "D:\\HadISD++\\Data\\PEI2015-2016.HourlyHdr.csv";
			
	vector<CStation> station_info;
	//Obtenir la liste des stations
	read_file(file,station_info);
	/*for (size_t i = 0, len = station_info.size(); i < len; i++)
		cout <<station_info.at(i).toString()<<endl;*/
	//Creer des fichiers netcdf à partir de chaque csv file.
	//ncdf(station_info);
	bool second = false;
	test mytest;
	internal_checks(station_info, mytest, second, DATES);
	std::system("PAUSE");
	return 0;
	
}




	