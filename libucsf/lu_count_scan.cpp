/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_count_scan.cpp                                             *
*                                                                             *
*  Created    : August 4th 2005                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <stdexcept>
#include <lg_string.h>
#include <lg_io.h>
#include <lgen_error.h>
#include <lu_count_scan.h>
#include <lu_xml_data.h>
using std::string;
using std::getline;
using std::runtime_error;

CountScans::~CountScans () {}

XMLCountScans::XMLCountScans ( PPExpatCountScans* ppesc, const string& filename ) :
	ppesc ( ppesc )
{
	try {
		ppesc->parseXMLFromFile ( filename );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
}
XMLCountScans::~XMLCountScans ()
{
	delete ppesc;
}
int XMLCountScans::getCount () const
{
	return ppesc->getCount ();
}

MascotCountScans::MascotCountScans ( const string& filename )
{
	count = 0;
	string line;
	GenIFStream ist ( filename );
	while ( getline ( ist, line ) ) {
		if ( !line.compare ( 0, 8, "END IONS" ) ) {
			count++;
		}
	}
}

MS2CountScans::MS2CountScans ( const string& filename )
{
	count = 0;
	string line;
	GenIFStream ist ( filename );
	while ( getline ( ist, line ) ) {
		if ( !line.compare ( 0, 1, "S" ) ) {
			count++;
		}
	}
}

APLCountScans::APLCountScans ( const string& filename )
{
	count = 0;
	string line;
	GenIFStream ist ( filename );
	while ( getline ( ist, line ) ) {
		if ( !line.compare ( 0, 12, "peaklist end" ) ) {
			count++;
		}
	}
}

PPCountMSScans::PPCountMSScans ( const string& filename )
{
	count = 0;
	string line;
	GenIFStream ist ( filename );
	while ( getline ( ist, line ) ) {
		if ( !line.compare ( 0, 3, ">M1" ) ) {
			count++;
		}
	}
}

PPCountMSMSScans::PPCountMSMSScans ( const string& filename )
{
	count = 0;
	string line;
	GenIFStream ist ( filename );
	while ( getline ( ist, line ) ) {
		if ( !line.compare ( 0, 3, ">M2" ) ) {
			count++;
		}
	}
}
SpaceSeparatedSpectraCountScans::SpaceSeparatedSpectraCountScans ( const string& filename )
{
	count = 0;	// Spectra are separated by one or more blank lines which contain only space characters
	string line;
	GenIFStream ist ( filename );
	bool spectrumSearch = true;
	while ( getline ( ist, line ) ) {
		if ( spectrumSearch ) {
			if ( !genEmptyString ( line ) ) {
				count++;
				spectrumSearch = false;
			}
		}
		else {
			if ( genEmptyString ( line ) ) {
				spectrumSearch = true;
			}
		}
	}
}
