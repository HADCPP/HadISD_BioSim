#include "station.h"


#include <string>
#include <iostream>
#include <algorithm>

using namespace std;



	CStation::CStation(string id,string name, double lat, double lon, double elev,string wmo)
	{
		m_id = id;
		m_lon = lon;
		m_lat = lat;
		m_elev = elev;
		m_name = name;
		m_wmoid = wmo;
	}
	CStation::~CStation()
	{
	}
	string CStation::toString()
	{
		return "CStation  :" + m_id + "  "+ m_name+  "+  Longitude: " + to_string(m_lon) + "  Latitude:  " + to_string(m_lat) + "  Elevation  :" + to_string(m_elev);

	}
	