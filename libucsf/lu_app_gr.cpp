/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_app_gr.cpp                                                 *
*                                                                             *
*  Created    : October 20th 1997                                             *
*                                                                             *
*  Purpose    : Java graph applet interface functions.                        *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <iomanip>
#include <lg_string.h>
#include <lgen_file.h>
#include <lu_app_gr.h>
#include <lu_delim.h>
#include <lu_getfil.h>
#include <lu_param_list.h>
using std::string;
using std::ostream;
using std::endl;
using std::stringstream;
using std::fixed;
using std::setprecision;

static void writeParam ( ostream& os, const string& name, const string& value )
{
	os << "<param name=\"" << name << "\" value=\"" << genReplaceSubstrings ( value, "<", "&lt;" ) << "\" />" << endl;
}
static void writeParam ( ostream& os, const string& name, int index, const string& value )
{
	os << "<param name=\"" << name << index << "\" value=\"" << genReplaceSubstrings ( value, "<", "&lt;" ) << "\" />" << endl;
}
template <class T>
static void writeParam ( ostream& os, const string& name, const T& value )
{
	os << "<param name=\"" << name << "\" value=\"" << value << "\" />" << endl;
}
template <class T>
static void writeParam ( ostream& os, const string& name, int index, const T& value )
{
	os << "<param name=\"" << name << index << "\" value=\"" << value << "\" />" << endl;
}
static void writeParamDoubleFixed ( ostream& os, const string& name, double value, int precision )
{
	os << "<param name=\"" << name;
	os << "\" value=\"";
	genPrint ( os, value, precision );
	os << "\" />" << endl;
}
static void writeParamDoubleFixed ( ostream& os, const string& name, int index, double value, int precision )
{
	stringstream n;
	n << name << index;
	writeParamDoubleFixed ( os, n.str (), value, precision );
}
static void writeParamDoubleSigFig ( ostream& os, const string& name, double value, int precision )
{
	os << "<param name=\"" << name;
	os << "\" value=\"";
	if ( fabs ( value ) < 1e-10 ) os << fixed << setprecision ( 1 ) << 0.0;
	else genPrintSigFig ( os, value, precision );
	os << "\" />" << endl;
}
static void writeParamDoubleSigFig ( ostream& os, const string& name, int index, double value, int precision )
{
	stringstream n;
	n << name << index;
	writeParamDoubleSigFig ( os, n.str (), value, precision );
}
SpectrumGraph::SpectrumGraph ( const string& fileName )
{
	readParams ( fileName );
}
void SpectrumGraph::drawGraph ( ostream& os, GraphData& graphData, bool sorted, double xMin, double xMax, double heightMultiplier )
{
	os << "<applet";
	os << " ";
	os << "name=\"SpectrumGraph\"";
	os << " ";
	os << "code=\"SpectrumGraph.class\"";
	os << " ";
	os << "archive=\"";
	os << HTMLDir::instance ().getVirtualHTMLJavaDir ();
	os << "SpectrumGraph.jar\"";
	os << " ";
	os << "width=\"" << appletWidth << "\"";
	os << " ";
	int aHeight = appletHeight * heightMultiplier;
	os << "height=\"" << aHeight << "\"";
	os << ">" << endl;

//	writeParam ( os, "permissions", "sandbox" );
	if ( sorted == false ) writeParam ( os, "unsorted", "true" );
	writeParam ( os, "x_axis_label", xAxisLabel );
	writeParam ( os, "line_width", lineWidth );

	writeParam ( os, "applet_background_red", appletBackgroundColor.red );
	writeParam ( os, "applet_background_green", appletBackgroundColor.green );
	writeParam ( os, "applet_background_blue", appletBackgroundColor.blue );

	writeParam ( os, "axes_red", axesColor.red );
	writeParam ( os, "axes_green", axesColor.green );
	writeParam ( os, "axes_blue", axesColor.blue );

	writeParam ( os, "default_peak_red", defaultPeakColor.red );
	writeParam ( os, "default_peak_green", defaultPeakColor.green );
	writeParam ( os, "default_peak_blue", defaultPeakColor.blue );

	writeParam ( os, "num_application_colors", applicationColors.size () );
	for ( AppletColorVectorSizeType i = 0 ; i < applicationColors.size () ; i++ ) {
		stringstream red;
		red << "application" << i << "_red";
		writeParam ( os, red.str (), applicationColors [i].red );
		stringstream green;
		green << "application" << i << "_green";
		writeParam ( os, green.str (), applicationColors [i].green );
		stringstream blue;
		blue << "application" << i << "_blue";
		writeParam ( os, blue.str (), applicationColors [i].blue );
	}

	writeParam ( os, "default_font_family", defaultFont.family );
	writeParam ( os, "default_font_style", defaultFont.style );
	writeParam ( os, "default_font_points", defaultFont.points );

	writeParam ( os, "peak_label_font_family", peakLabelFont.family );
	writeParam ( os, "peak_label_font_style", peakLabelFont.style );
	writeParam ( os, "peak_label_font_points", peakLabelFont.points );

	graphData.draw ( os, xMin, xMax );

	os << "</applet>" << endl;
}
void SpectrumGraph::readParams ( const string& fileName )
{
	GenNameValueStream nvs ( MsparamsDir::instance ().getParamPath ( fileName ) );

	appletWidth = nvs.getIntValue ( "applet_width" );
	appletHeight = nvs.getIntValue ( "applet_height" );
	lineWidth = nvs.getIntValue ( "line_width" );

	appletBackgroundColor.red	= nvs.getIntValue ( "applet_background_color_red" );
	appletBackgroundColor.green	= nvs.getIntValue ( "applet_background_color_green" );
	appletBackgroundColor.blue	= nvs.getIntValue ( "applet_background_color_blue" );

	axesColor.red	= nvs.getIntValue ( "axes_color_red" );
	axesColor.green	= nvs.getIntValue ( "axes_color_green" );
	axesColor.blue	= nvs.getIntValue ( "axes_color_blue" );

	defaultPeakColor.red	= nvs.getIntValue ( "default_peak_color_red" );
	defaultPeakColor.green	= nvs.getIntValue ( "default_peak_color_green" );
	defaultPeakColor.blue	= nvs.getIntValue ( "default_peak_color_blue" );

	int numApplicationColors = nvs.getIntValue ( "number_application_colors" );
	applicationColors.resize ( numApplicationColors );
	for ( int i = 0 ; i < numApplicationColors ; i++ ) {
		string num = gen_itoa ( i+1 );
		applicationColors [i].red	= nvs.getIntValue ( "application_color_" + num + "_red" );
		applicationColors [i].green	= nvs.getIntValue ( "application_color_" + num + "_green" );
		applicationColors [i].blue	= nvs.getIntValue ( "application_color_" + num + "_blue" );
	}
	defaultFont.family	= nvs.getStringValue ( "default_font_family " );
	defaultFont.style	= nvs.getStringValue ( "default_font_style" );
	defaultFont.points	= nvs.getIntValue ( "default_font_points" );

	peakLabelFont.family	= nvs.getStringValue ( "peak_label_font_family " );
	peakLabelFont.style		= nvs.getStringValue ( "peak_label_font_style" );
	peakLabelFont.points	= nvs.getIntValue ( "peak_label_font_points" );

	xAxisLabel = nvs.getStringValue ( "x_axis_label" );
}
bool SpectrumGraph::drawGraphFlag = false;
void SpectrumGraph::setDrawGraph ( bool flag )
{
	drawGraphFlag = flag;
}
bool SpectrumGraph::getDrawGraph ()
{
	return drawGraphFlag;
}
int GraphData::xPrecision = 4;
int GraphData::yPrecision = 3;
GraphData::GraphData ( const XYData& xyData ) :
	XYData ( xyData )
{
}
void GraphData::draw ( ostream& os, double xMin, double xMax ) const
{
	writeParam ( os, "num_peaks", size () );
	writeRangeParams ( os, xMin, xMax );
	for ( first () ; isDone () ; next () ) {
		writeParamDoubleFixed ( os, "x", index (), x (), xPrecision );
		writeParamDoubleSigFig ( os, "y", index (), y (), yPrecision );
	}
}
void GraphData::writeXML ( ostream& os ) const
{
	os << "<graph_data>" << endl;
		ParameterList::printXML ( os, "num_peaks", size () );
		for ( first () ; isDone () ; next () ) {
			os << "<point";
			os << " ";
			os << "x=\""; 
			genPrint ( os, x (), xPrecision );
			os << "\" ";
			os << "y=\""; 
			genPrintSigFig ( os, y (), yPrecision );
			os << "\" />";
			os << endl;
		}
	os << "</graph_data>" << endl;
}
void GraphData::writeTabDelimitedText ( ostream& os ) const
{
	delimitedRowStart ( os );
		delimitedCell ( os, "Num Peaks:" );
		delimitedCell ( os, (int)size () );
	delimitedRowEnd ( os );
	delimitedRowStart ( os );
		delimitedHeader ( os, "X Value" );
		delimitedHeader ( os, "Y Value" );
	delimitedRowEnd ( os );
	for ( first () ; isDone () ; next () ) {
		delimitedRowStart ( os );
			delimitedCell ( os, x (), xPrecision );
			delimitedCellSigFig ( os, y (), yPrecision );
		delimitedRowEnd ( os );
	}
}
void GraphData::writeRangeParams ( ostream& os, double xMin, double xMax ) const
{
	if ( size () ) {
		xMin = xMin == 0.0 ? minX () : xMin;
		xMax = xMax == 0.0 ? maxX () : xMax;
		writeParamDoubleFixed ( os, "start_x", xMin, xPrecision );
		if ( xMin == xMax )
			writeParamDoubleFixed ( os, "end_x", xMax + 1, xPrecision );
		else
			writeParamDoubleFixed ( os, "end_x", xMax, xPrecision );
	}
}
FileGraphData::FileGraphData ( const XYData& xyData, bool retainZeroPeaks ) :
	GraphData ( xyData ),
	retainZeroPeaks ( retainZeroPeaks )
{
}
void FileGraphData::draw ( ostream& os, double xMin, double xMax ) const
{
	writeParam ( os, "url_flag", "true" );
	writeRangeParams ( os, xMin, xMax );
	PPTempFile tempFile ( "graph", ".txt" );
	writeParam ( os, "data_url", tempFile.getURL () );
	GenOFStream ofs ( tempFile.getAdjustedPath () );
	int num = 0;
	for ( first () ; isDone () ; next () ) {
		if ( retainZeroPeaks || y () != 0.0 ) {
			genPrint ( ofs, x (), xPrecision );
			ofs << " ";
			genPrintSigFig ( ofs, y (), yPrecision );
			ofs << endl;
			num++;
		}
	}
	writeParam ( os, "num_peaks", num );
}
ColouredGraphData::ColouredGraphData ( const XYData& xyData, int colour ) :
	GraphData ( xyData )
{
	for ( first () ; isDone () ; next () ) {
		colourList.push_back ( colour );
	}
}
void ColouredGraphData::add ( const double& x, const double& y, const int& colour )
{
	XYData::add ( x, y );
	colourList.push_back ( colour );
}
void ColouredGraphData::draw ( ostream& os, double xMin, double xMax ) const
{
	GraphData::draw ( os, xMin, xMax );
	if ( size () ) writeParam ( os, "labels", "true" );
	for ( first () ; isDone () ; next () ) {
		if ( colourList [index ()] ) {
			writeParam ( os, "c", index (), colourList [index ()] );
		}
	}
}
void LabelledGraphData::add ( const double& x, const double& y, const string& label, const int& colour )
{
	XYData::add ( x, y );
	colourList.push_back ( colour );
	labelList.push_back ( label );
}
void LabelledGraphData::draw ( ostream& os, double xMin, double xMax ) const
{
	GraphData::draw ( os, xMin, xMax );
	if ( size () ) writeParam ( os, "labels", "true" );
	for ( first () ; isDone () ; next () ) {
		if ( colourList [index ()] ) {
			writeParam ( os, "c", index (), colourList [index ()] );
		}
		if ( labelList [index ()] != "" ) {
			writeParam ( os, "l", index (), labelList [index ()] );
		}
	}
}
LabelledCatagorizedGraphData::LabelledCatagorizedGraphData ( const StringVector& categories, const BoolDeque& showCategory ) :
	categories ( categories ),
	showCategory ( showCategory )
{
}
void LabelledCatagorizedGraphData::add ( const double& x, const double& y, const string& label, const int& colour, const string& category )
{
	LabelledGraphData::add ( x, y, label, colour );
	catagoryList.push_back ( category );
}
void LabelledCatagorizedGraphData::draw ( ostream& os, double xMin, double xMax ) const
{
	LabelledGraphData::draw ( os, xMin, xMax );
	writeCategories ( os );
	defineCheckboxes ( os );
}
void LabelledCatagorizedGraphData::writeCategories ( ostream& os ) const
{
	if ( size () )	writeParam ( os, "categories", "true" );
	for ( first () ; isDone () ; next () ) {
		if ( catagoryList [index ()] != "" ) {
			writeParam ( os, "g", index (), catagoryList [index ()] );
		}
	}
}
void LabelledCatagorizedGraphData::defineCheckboxes ( ostream& os ) const
{
	int numCategories = categories.size ();
	writeParam ( os, "num_categories", numCategories );
	for ( int i = 0 ; i < numCategories ; i++ ) {
		string num = gen_itoa ( i );
		writeParam ( os, "cl"+num, categories [i] );
		if ( showCategory [i] ) writeParam ( os, "sc"+num, "true" );
	}
}
LabelledCatagorizedFileGraphData::LabelledCatagorizedFileGraphData ( const StringVector& categories, const BoolDeque& showCategory, const XYData& xyData, bool retainZeroPeaks ) :
	LabelledCatagorizedGraphData ( categories, showCategory ),
	retainZeroPeaks ( retainZeroPeaks )
{
	for ( xyData.first () ; xyData.isDone () ; xyData.next () ) {	// Write the URL data
		XYData::add ( xyData.x (), xyData.y () );
	}
	fileGraphSize = xyData.size ();		// Initial data size
}
void LabelledCatagorizedFileGraphData::draw ( ostream& os, double xMin, double xMax ) const
{
	writeParam ( os, "url_flag", "true" );					// Start URL
	PPTempFile tempFile ( "graph", ".txt" );
	writeParam ( os, "data_url", tempFile.getURL () );
	GenOFStream ofs ( tempFile.getAdjustedPath () );
	int num = 0;
	for ( first () ; isDone ( fileGraphSize ) ; next () ) {	// Write the URL data
		if ( retainZeroPeaks || y () != 0.0 ) {
			genPrint ( ofs, x (), xPrecision );
			ofs << " ";
			genPrintSigFig ( ofs, y (), yPrecision );
			ofs << endl;
			num++;
		}
	}
	writeRangeParams ( os, xMin, xMax );
	bool sz = false;
	for ( int i = 0 ; isDone () ; next (), i++ ) {
		sz = true;
		writeParamDoubleFixed ( os, "x", index (), x (), xPrecision );
		writeParamDoubleSigFig ( os, "y", index (), y (), yPrecision );
		if ( colourList [i] ) {
			writeParam ( os, "c", index (), 8 );
		}
		else {
			writeParam ( os, "c", index (), 9 );
		}
		if ( labelList [i] != "" ) {
			writeParam ( os, "l", index (), labelList [i] );
		}
		if ( catagoryList [i] != "" ) {
			writeParam ( os, "g", index (), catagoryList [i] );
		}
		num++;
	}
	if ( sz ) {
		writeParam ( os, "labels", "true" );
		writeParam ( os, "categories", "true" );
	}
	writeParam ( os, "num_peaks", num );					// End URL
	defineCheckboxes ( os );
}
