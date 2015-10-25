/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_r_plot.cpp                                                 *
*                                                                             *
*  Created    : March 8th 2006                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2006-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_stdlib.h>
#include <lgen_file.h>
#include <lu_getfil.h>
#include <lu_r_plot.h>
#include <lu_param_list.h>
using std::string;
using std::ostream;
using std::endl;

/*
struct PinkColormap {
	static double red [];
	static double green [];
	static double blue [];
	static int size;
	static StringVector getRColorMap ();
};
double PinkColormap::red [] = {
	0.1179,0.1959,0.2507,0.2955,0.3343,0.3691,0.4009,0.4303,0.4579,0.4839,0.5085,0.5320,0.5546,0.5762,0.5971,0.6172,
	0.6367,0.6557,0.6741,0.6920,0.7094,0.7265,0.7431,0.7594,0.7664,0.7732,0.7800,0.7868,0.7935,0.8001,0.8067,0.8133,
	0.8197,0.8262,0.8325,0.8389,0.8452,0.8514,0.8576,0.8637,0.8698,0.8759,0.8819,0.8879,0.8938,0.8997,0.9056,0.9114,
	0.9172,0.9230,0.9287,0.9344,0.9400,0.9456,0.9512,0.9567,0.9623,0.9677,0.9732,0.9786,0.9840,0.9894,0.9947,1.0000 };
double PinkColormap::green [] = { 
	0.0000,0.1029,0.1455,0.1782,0.2057,0.2300,0.2520,0.2722,0.2910,0.3086,0.3253,0.3412,0.3563,0.3709,0.3849,0.3984,
	0.4115,0.4241,0.4364,0.4484,0.4600,0.4714,0.4825,0.4933,0.5175,0.5407,0.5628,0.5842,0.6048,0.6247,0.6440,0.6627,
	0.6809,0.6986,0.7159,0.7328,0.7493,0.7655,0.7813,0.7968,0.8120,0.8270,0.8416,0.8560,0.8702,0.8842,0.8979,0.9114,
	0.9172,0.9230,0.9287,0.9344,0.9400,0.9456,0.9512,0.9567,0.9623,0.9677,0.9732,0.9786,0.9840,0.9894,0.9947,1.0000 };
double PinkColormap::blue [] = {
	0.0000,0.1029,0.1455,0.1782,0.2057,0.2300,0.2520,0.2722,0.2910,0.3086,0.3253,0.3412,0.3563,0.3709,0.3849,0.3984,
	0.4115,0.4241,0.4364,0.4484,0.4600,0.4714,0.4825,0.4933,0.5040,0.5143,0.5245,0.5345,0.5443,0.5540,0.5634,0.5727,
	0.5819,0.5909,0.5998,0.6086,0.6172,0.6257,0.6341,0.6424,0.6506,0.6587,0.6667,0.6746,0.6824,0.6901,0.6977,0.7052,
	0.7272,0.7485,0.7692,0.7893,0.8090,0.8282,0.8469,0.8653,0.8832,0.9008,0.9181,0.9351,0.9517,0.9681,0.9842,1.0000 };
int PinkColormap::size = 64;

StringVector PinkColormap::getRColorMap ()
{
	StringVector sv;
	for ( int i = 0 ; i < size ; i++ ) {
		char color [8];
		int r = red [i] * 255;
		int g = green [i] * 255;
		int b = blue [i] * 255;
		sprintf ( color, "#%.2X%.2X%.2X", r, g, b );
		sv.push_back ( color );
		std::cout << sv.back () << "<br />" << std::endl;
	}
	return sv;
}
*/
string RPlot::getRCommand ()
{
	return InfoParams::instance ().getRCommand ();
}

string RPlot::rCommand = getRCommand ();
bool RPlot::rFlag = !rCommand.empty ();
bool RPlot::keepRDataFile = InfoParams::instance ().getBoolValue ( "keep_r_data_file" );

RPlot::RPlot ( const string& rScriptName ) :
	rScriptName ( rScriptName )
{
	if ( rFlag ) {
		PPTempFile dataTempFile ( "data", ".txt" );
		dataFileFullPath = dataTempFile.getFullPath ();
	}
}
RPlot::~RPlot ()
{
	if ( rFlag && !keepRDataFile ) {
		genUnlink ( dataFileFullPath );
	}
}
void RPlot::printImage ( ostream& os, const string& imageType ) const
{
	if ( rFlag ) {
		if ( genFileSize ( dataFileFullPath ) != 0 ) {
			PPTempFile imageTempFile ( "image", imageType );
			string imageTempFileFullPath = imageTempFile.getFullPath ();
			string url = imageTempFile.getURL ();

			string command;
#ifdef VIS_C
			command += "\"";
			command += "\"" + rCommand + "\"";
#else
			command += rCommand;
#endif
			command += " --vanilla --slave --quiet --args ";
#ifdef VIS_C
			command += "\"" + dataFileFullPath + "\"";
#else
			command += dataFileFullPath;
#endif
			command += " ";
#ifdef VIS_C
			command += "\"" + imageTempFileFullPath + "\"";
#else
			command += imageTempFileFullPath;
#endif
			command += " ";
			command += "< ";
			command +=  rScriptName;
#ifdef VIS_C
			command +=  "\"";
#endif
			int ret = genSystem ( command, SystemCallBinDir::instance ().getSystemCallBinDir (), true );
			bool imageDone = ( ret == 0 );
			if ( imageType == ".pdf" ) {
				if ( imageDone ) {
					os << "<a href=\"";
					os << url;
					os << "\">";
					os << "PDF Image</a>";
				}
			}
			else {
				if ( imageDone ) {
					os << "<img src=\"";
					os << url;
					os << "\"";
					os << " ";
					os << "alt=\"Missing Image\"";
					os << " />";
				}
				else os << "<b>Image not plotted.</b>" << endl;
			}
		}
	}
}
void RPlot::printImageAndLink ( ostream& os ) const
{
	printImage ( os, ".png" );
	os << "<br />" << endl;
	//printImage ( os, ".pdf" );
	//os << "<br />" << endl;
}
