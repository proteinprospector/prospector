/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_data.cpp                                                   *
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
*  Copyright (2005-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef VIS_C
#include <stdexcept>
#endif
#include <iomanip>
#include <nr.h>
#include <lg_string.h>
#include <lg_time.h>
#include <lgen_error.h>
#include <lu_delim.h>
#ifdef BATCHTAG
#include <lu_xml_data.h>
#endif
#include <lu_charge.h>
#include <lu_mass_conv.h>
#include <lu_mgf.h>
#include <lu_t2d.h>
#include <lu_html.h>
#include <lu_cgi_val.h>
#include <lu_spec_id.h>
#include <lu_param_list.h>
#include <lu_getfil.h>
using std::string;
using std::cout;
using std::streampos;
using std::getline;
using std::istream;
using std::ostream;
using std::istringstream;
using std::ostringstream;
using std::ios;
using std::endl;
using std::runtime_error;

DataPoint::~DataPoint ()
{
}
string DataPoint::getSpotID () const
{
	ostringstream ost;
	ost << fraction << "-" << spot;
	return ost.str ();
}
string DataPoint::getSpotRunID () const
{
	ostringstream ost;
	ost << fraction << "-" << spot << "-" << run;
	return ost.str ();
}
string DataPoint::getSpecID () const
{
	ostringstream ost;
	ost << fraction << "-" << spot << "-" << run;
	return ost.str ();
}

double MSMSDataPoint::getParentMPlusH () const
{
	return mOverZToMPlusH ( precursorMZ, precursorCharge, true );
}
void MSMSDataPoint::calibrate ( const Tolerance* tol, double offset )
{
	precursorMZ -= tol->getCorrection ( precursorMZ, offset );
}
void MSMSDataPoint::setTolerance ( const Tolerance* tol )
{
	precursorTolerance = tol->getTolerance ( precursorMZ, precursorCharge );
}
void MSMSDataPoint::printSpecInfoDelimited ( ostream& os ) const
{
	ostringstream ost;
	ost << fraction << "-" << spot << "-" << run << "-" << spectrumNumber;
	if ( !msmsInfo.empty () ) {
		ost << "-" << msmsInfo;
	}
	delimitedCell ( os, ost.str () );
}
string MSMSDataPoint::getSpecID () const
{
	ostringstream ost;
	ost << fraction << "-" << spot << "-" << run << "-" << spectrumNumber;
	return ost.str ();
}
void MSMSDataPoint::putCGI ( std::ostream& os ) const
{
	printCGI ( os, "fraction", fraction );
	printCGI ( os, "spot_number", spot );
	printCGI ( os, "run", run );
	printCGI ( os, "spectrum_number", spectrumNumber );
	if ( !msmsInfo.empty () ) printCGI ( os, "msms_info", msmsInfo );
}
SpectrumRange::SpectrumRange ( const ParameterList* params ) :
	startFraction	( params->getIntValue ( "start_fraction" ) ),
	startSpectrum	( params->getIntValue ( "start_spectrum" ) ),
	endFraction		( params->getIntValue ( "end_fraction" ) ),
	endSpectrum		( params->getIntValue ( "end_spectrum" ) ),
	allFractions	( startFraction == 0 )
{
}
string DataReader::lastFilename;
string DataReader::lastVersion;
const int DataReader::MAX_CHARGE = 30;
DataReader::DataReader ( const string& s, bool fileFlag, bool intensityFlag ) :
	ifstr ( initFStream ( s, fileFlag ) ),
	isstr ( fileFlag ? "" : s ),
	istr ( initStream ( fileFlag ) ),
	intensityFlag ( intensityFlag ),
	num ( 0 ),
	endPos ( 0 )
{
	if ( fileFlag ) lastFilename = s;
}
DataReader::~DataReader ()
{
	delete ifstr;
}
istream& DataReader::initStream ( bool fileFlag )
{
	if ( fileFlag ) return *ifstr;
	else return isstr;
}
GenIFStream* DataReader::initFStream ( const string& s, bool fileFlag )
{
	if ( fileFlag ) {
		GenIFStream p ( s, ios::binary );
		string line;
		getline ( p, line );
		if ( line [line.length ()-1] == 13 ) return new GenIFStream ( s );	// DOS text file
		else return new GenIFStream ( s, ios::binary );						// UNIX text file
	}
	else return 0;
}
void DataReader::readParentData ( int charge, int isotopeOffset )
{
	string line;
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 ) {
			if ( line [0] != '#' ) {
				istringstream ist ( line );
				double mOZ;
				if ( ist >> mOZ ) {
					if ( charge == 0 ) {			// Read the charge from the file
						if ( ( ist >> charge ) == 0 ) charge = 1;
					}
					if ( isotopeOffset != 0 ) mOZ += ( isotopeOffset * 1.0029 / charge );
					msmsDataPoint.setPrecursor ( mOZ, charge, 100.0 );
					return;
				}
			}
		}
	}
}
void DataReader::skipParentData ()
{
	string line;
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 ) {
			if ( line [0] != '#' ) {
				istringstream ist ( line );
				double mOZ;
				if ( ist >> mOZ ) {
					return;
				}
			}
		}
	}
}
void DataReader::readData ()
{
	dataPoint->clear ();
	string line;

	int c;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( c == '>' || c == 'E' ) return;
		if ( !getline ( istr, line ) ) return;
		if ( line.length () != 0 ) {
			if ( isdigit ( line [0] ) ) readLine ( const_cast <char*> (line.c_str ()) );
			else {
				line = gen_strtrim ( line );	// Try trimming the line
				if ( isdigit ( line [0] ) ) readLine ( const_cast <char*> (line.c_str ()) );
			}
		}
	}
}
void DataReader::readLine ( char* ptr )
{
	double mass = atof ( ptr );
	while ( *ptr && isspace ( *ptr ) ) ptr++; // Skip leading white space
	while ( *ptr && !isspace ( *ptr ) ) ptr++; // Skip mass
	double intensity = 100.0;
	if ( intensityFlag ) {
		intensity = atof ( ptr );
		while ( *ptr && isspace ( *ptr ) ) ptr++; // Skip leading white space
		while ( *ptr && !isspace ( *ptr ) ) ptr++; // Skip intensity
	}
	while ( *ptr && isspace ( *ptr ) ) ptr++; // Skip leading white space
	int charge = 1;
	if ( *ptr ) {
		charge = atoi ( ptr );
		if ( charge > MAX_CHARGE ) ErrorHandler::genError ()->error ( "Maximum charge exceeded - possibly the Data Format parameter is incorrectly set.\n" );
	}
	dataPoint->addPeak ( mass, charge, intensity );
}
void DataReader::readMSData ( MSDataPointVector& msDataPointList )
{
	dataPoint = &msDataPoint;
	readData ();
	msDataPointList.push_back ( msDataPoint );
}
void DataReader::readMSMSData ( MSMSDataPointVector& msmsDataPointList, int charge, int isotopeOffset )
{
	readParentData ( charge, isotopeOffset );
	dataPoint = &msmsDataPoint;
	readData ();
	if ( msmsDataPoint.getPrecursorCharge () == 0 ) {
		int spectrumNumber = msmsDataPoint.getSpectrumNumber ();
		const IntVector& zRange = msmsDataPoint.getPrecursorChargeRange ();
		if ( zRange.empty () ) {
			for ( IntVectorSizeType i = 0 ; i < chargeRange.size () ; i++ ) {
				msmsDataPoint.setPrecursorCharge ( chargeRange [i] );
				int j = 10000 * chargeRange [i];
				msmsDataPoint.setSpectrumNumber ( j + spectrumNumber );
				msmsDataPointList.push_back ( msmsDataPoint );
			}
		}
		else {
			for ( IntVectorSizeType i = 0 ; i < zRange.size () ; i++ ) {
				msmsDataPoint.setPrecursorCharge ( zRange [i] );
				int j = 10000 * zRange [i];
				msmsDataPoint.setSpectrumNumber ( j + spectrumNumber );
				msmsDataPointList.push_back ( msmsDataPoint );
			}
		}
	}
	else
		msmsDataPointList.push_back ( msmsDataPoint );
}
bool DataReader::readMSMSData ( MSMSDataPointVector& msmsDataPointList, int charge, int isotopeOffset, double mOverZ )
{
	readParentData ( charge, isotopeOffset );
	if ( msmsDataPoint.getPrecursorMZ () != mOverZ ) return false;
	dataPoint = &msmsDataPoint;
	readData ();
	msmsDataPointList.push_back ( msmsDataPoint );
	return true;
}
void DataReader::readMSMSQuanData ( XYData& xyData, double startMass, double endMass )
{
	skipParentData ();
	dataPoint = &msmsDataPoint;
	readData ();
	DataFilePeakVector& dataPeaks = dataPoint->getDataPeaks ();
	for ( int i = 0 ; i < dataPeaks.size () ; i++ ) {
		double mz = dataPeaks [i].getMOverZ ();
		int z = dataPeaks [i].getCharge ();
		if ( z != 1 ) continue;
		if ( mz < startMass ) continue;
		if ( mz > endMass ) break;
		xyData.add ( mz, dataPeaks [i].getIntensity () );
	}
}
bool DataReader::writePeakList ( ostream& os, const SpecID& specID, const string& version )
{
	throw runtime_error ( "Writing peak lists in this format is currently not supported." );
}
bool DataReader::getMSMSQuanData ( XYData& xyData, const SpecID& specID, double startMass, double endMass, const string& version )
{
	throw runtime_error ( "MSMS quantitation doesn't work with this type of centroid file." );
}
void DataReader::rewind ()
{
	istr.seekg ( spos );
}
bool DataReader::isEOF ()
{
	int c = istr.peek ();
	return c == EOF;
}
void DataReader::getData ( MSMSDataPointVector& msmsDataPointList, const SpecID& specID, int chrg, const string& version, UpdatingJavascriptMessage* ujm )
{
	throw runtime_error ( "Writing peak lists in this format is currently not supported." );
}
#ifdef BATCHTAG
XMLDataReader::XMLDataReader ( const string& s, bool fileFlag, const string& dataFormat ) :
	DataReader ( s, fileFlag, true ),
	s ( s ),
	fileFlag ( fileFlag ),
	dataFormat ( dataFormat )
{
}
void XMLDataReader::getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const string& version, int off )
{
	PPExpat* ppe;
	if ( dataFormat == "mzXML" ) {
		ppe = new PPExpatMZXMLData ( msmsDataPointList, fraction, spectrumRange, specID, chrg, chargeRange );
	}
	else if ( dataFormat == "mzData" ) {
		ppe = new PPExpatMZDataData ( msmsDataPointList, fraction, spectrumRange, specID, chrg, chargeRange );
	}
	else {
	}
	try {
		if ( fileFlag )
			ppe->parseXMLFromFile ( s );
		else
			ppe->parseXMLFromString ( s );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
	delete ppe;
}
#endif
MS2DataReader::MS2DataReader ( const string& s, bool fileFlag ) :
	DataReader ( s, fileFlag, true ),
	readScan ( true ),
	readRT ( true )
{
}
bool MS2DataReader::parseHeader () const
{
	bool bullseye = false;
	int c;
	string line;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( c == 'S' ) {  //		Scan line - probably no header
			break;
		}
		getline ( istr, line );
		if ( line.length () != 0 ) {
			c = line [0];
			if ( c == 'H' ) {
				istringstream istr ( line );
				string dummy;
				istr >> dummy;				// Skip the H
				string s;
				istr >> s;
				if ( s == "FileGenerator" ) {
					string t;
					istr >> t;
					if ( t == "Bullseye" ) bullseye = true;
				}
				break;
			}
			else
				break;
		}
	}
	return bullseye;
}
void MS2DataReader::parseHeaderLines ( bool bullseye, string& spot, int& spectrumNumber, string& msmsInfo, DoubleVector& mVector, IntVector& zVector )
{
	string startScan;
	string endScan;
	string rt;
	int c;
	string line;
	int numS = 0;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( isdigit ( c ) ) {
			break;			// data section
		}
		if ( numS == 1 && c == 'S' ) {  //		2nd scan line, no data
			break;			// data section
		}
		getline ( istr, line );
		if ( line.length () != 0 ) {
			c = line [0];
			if ( c == 'S' ) {
				istringstream istr ( line );
				string dummy;
				istr >> dummy;
				istr >> startScan;
				istr >> endScan;
				numS++;
			}
			else if ( c == 'I' ) {
				istringstream istr ( line );
				string dummy;
				istr >> dummy;
				string temp;
				istr >> temp;
				if ( temp == "RTime" ) {
					istr >> rt;
				}
				//else if ( temp == "EZ" ) {
				//	string temp;
				//	istr >> temp;
				//	istr >> temp;
				//	istr >> temp;
				//	istr >> intVal;
				//	iVector.push_back ( intVal );
				//}
			}
			else if ( c == 'Z' ) {
				istringstream istr ( line );
				string dummy;
				istr >> dummy;
				int ch;
				istr >> ch;
				double m;
				istr >> m;
				zVector.push_back ( ch );
				if ( bullseye )	mVector.push_back ( mPlusHToMOverZ ( m, ch, true ) );
				else			mVector.push_back ( mToMOverZ ( m, ch, true ) );
			}
		}
	}
	spot = (!readRT || rt.empty ()) ? startScan : rt;
	if ( startScan == endScan )	msmsInfo = startScan;
	else						msmsInfo = startScan + "," + endScan;
	MapStringToIntIterator cur = spotNumber.find ( spot );
	if ( cur != spotNumber.end () ) cur->second++;
	else spotNumber [spot] = 1;
	cur = spotNumber.find ( spot );
	spectrumNumber = cur->second;
}
void MS2DataReader::getDataFromIndex ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, int index, int chrg, int off )
{
	string id = gen_itoa ( fraction ) + "-" + gen_itoa ( index ) + "-1-1";
	SpecID specID ( id );
	readRT = false;
	readScan = false;
	getData ( msmsDataPointList, fraction, spectrumRange, specID, chrg, "", off );
}
void MS2DataReader::getDataFromScanNumber ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const string& scanNumber, int chrg, int off )
{
	string id = gen_itoa ( fraction ) + "-" + scanNumber + "-1-1";
	SpecID specID ( id );
	readRT = false;
	getData ( msmsDataPointList, fraction, spectrumRange, specID, chrg, "", off );
}
void MS2DataReader::getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const string& version, int off )
{
	if ( off != 0 ) istr.seekg ( off );
	bool readPeaks = false;
	bool bullseye = parseHeader ();
	for ( ; ; ) {
		string spot = gen_itoa ( num + 1 );
		string msmsInfo;
		DoubleVector mVector;
		IntVector zVector;
		int spectrumNumber;
		parseHeaderLines ( bullseye, spot, spectrumNumber, msmsInfo, mVector, zVector );
		if ( readScan == false ) spot = gen_itoa ( num + 1 );
		num++;
		if ( spectrumRange.stopReadingFile ( fraction, num ) ) return;
		if ( spectrumRange.startReadingFile ( fraction, num ) ) {
			int run = 1;
			readPeaks = specID.isMSMSSpectrumToRead100 ( fraction, spot, run, spectrumNumber );
			if ( readPeaks ) {
				msmsDataPoint.setPointInfo ( fraction, spot, run, spectrumNumber, msmsInfo );
				dataPoint = &msmsDataPoint;
				double intensity = 100.0;
				msmsDataPoint.setPrecursorIntensity ( intensity );
				bool moreData = readMSMSData ();
				if ( specID.isFilterSpectrum () ) {		// Single spectrum
					int idx = specID.getSpecNum ();
					idx /= 100;
					msmsDataPoint.setPrecursorMZ ( mVector [idx] );
					if ( chrg != 0 )	msmsDataPoint.setPrecursorCharge ( chrg );
					else				msmsDataPoint.setPrecursorCharge ( zVector [idx] );
					msmsDataPointList.push_back ( msmsDataPoint );
					return;
				}
				else {
					for ( int i = 0 ; i < mVector.size () ; i++ ) {
						msmsDataPoint.setPrecursorMZ ( mVector [i] );
						msmsDataPoint.setPrecursorCharge ( zVector [i] );
						int j = 100 * i;
						msmsDataPoint.setSpectrumNumber ( j + spectrumNumber );
						msmsDataPointList.push_back ( msmsDataPoint );
					}
				}
				readPeaks = false;
				if ( moreData == false ) return;	// End of file
			}
			else {
				if ( skipMSMSData () == false ) break;	// End of file
			}
		}
		else {
			if ( skipMSMSData () == false ) break;	// End of file
		}
	}
}
bool MS2DataReader::readMSMSData ()
{
	dataPoint->clear ();
	string line;

	int c;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( isupper ( c ) ) return true;				// More scans
		if ( !getline ( istr, line ) ) return false;	// No more scans
		if ( line.length () != 0 ) {
			if ( isdigit ( line [0] ) ) readLine ( const_cast <char*> (line.c_str ()) );
			else {
				line = gen_strtrim ( line );	// Try trimming the line
				if ( isdigit ( line [0] ) ) readLine ( const_cast <char*> (line.c_str ()) );
			}
		}
	}
	return false;	// No more scans
}
bool MS2DataReader::skipMSMSData ()
{
	string line;
	int c;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( isupper ( c ) ) return true;				// More scans
		if ( !getline ( istr, line ) ) return false;	// No more scans
	}
	return false;	// No more scans
}
void MS2DataReader::writePeakList ( ostream& os )
{
	int c;
	int numS = 0;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( numS == 1 && c == 'S' ) {  //		2nd scan line, no data
			break;			// data section
		}
		string line;
		if ( !getline ( istr, line ) ) {
			break;
		}
		if ( line.length () != 0 ) {
			if ( line [line.length ()-1] == '\r' ) line = line.substr ( 0, line.length ()-1 );
			if ( line [0] == 'S' ) {
				numS++;
			}
			os << line << endl;
		}
	}
}
bool MS2DataReader::writePeakList ( ostream& os, const SpecID& specID, const string& version )
{
	string sID = specID.getPureSpecID100 ();
	MapStringToStreamposConstIterator iter = mapSpecID.find ( sID );
	if ( iter != mapSpecID.end () ) {			// Scan position already found
		streampos sp = iter->second;
		if ( sp > 0 ) {
			istr.clear ();
			istr.seekg ( sp );
			writePeakList ( os );
			mapSpecID [sID] = -1;			// To indicate that the scan has already been written
		}
		return true;
	}
	istr.seekg ( endPos );					// Scan not found, go to the furthest reached so far
	for ( ; ; ) {
		string spot = gen_itoa ( num + 1 );
		string msmsInfo;
		DoubleVector mVector;
		IntVector zVector;
		int spectrumNumber;
		spos = istr.tellg ();
		parseHeaderLines ( false, spot, spectrumNumber, msmsInfo, mVector, zVector );
		int run = 1;
		int fraction = specID.getFraction ();
		ostringstream ost;
		ost << fraction << '-' << spot << '-' << run << '-' << spectrumNumber;
		mapSpecID [ost.str ()] = spos;
		num++;
		if ( specID.isMSMSSpectrumToRead100 ( spot, run, spectrumNumber ) ) {
			istr.seekg ( spos );
			writePeakList ( os );
			endPos = istr.tellg ();
			mapSpecID [ost.str ()] = -1;
			return true;
		}
		if ( skipMSMSData () == false ) break;	// End of file
	}
	return false;
}
APLDataReader::APLDataReader ( const string& s, bool fileFlag ) :
	DataReader ( s, fileFlag, true )
{
}
//header=RawFile: OR9_20130607_EDG_13EDEG003_Ti-IMAC_Col1A Index: 26135 Precursor: 0 _multi_
//header=RawFile: OR9_20130607_EDG_13EDEG003_Ti-IMAC_Col2C Index: 38476 Silind: 264
void APLDataReader::parseHeaderLines ( string& spot, int& spectrumNumber, string& msmsInfo, DoubleVector& mVector, IntVector& zVector )
{
	int c;
	string line;
	int numS = 0;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( isdigit ( c ) ) {
			break;			// data section
		}
		if ( numS == 1 && c == 'p' ) {  //		2nd line, no data
			break;			// data section
		}
		getline ( istr, line );
		if ( line.length () != 0 ) {
			c = line [0];
			if ( c == 'm' ) {					// eg mz=1065.41113255839
				mVector.push_back ( atof ( line.substr ( 3 ).c_str () ) );
			}
			else if ( c == 'c' ) {				// charge=4
				zVector.push_back ( atoi ( line.substr ( 7 ).c_str () ) );
			}
			else if ( c == 'h' ) {
				int start = line.find ( "Index: " ) + 7;
				spot = gen_itoa ( atoi ( line.substr ( start ).c_str () ) );
			}
		}
	}
	msmsInfo = spot;
	MapStringToIntIterator cur = spotNumber.find ( spot );
	if ( cur != spotNumber.end () ) cur->second++;
	else spotNumber [spot] = 1;
	cur = spotNumber.find ( spot );
	spectrumNumber = cur->second;
}
void APLDataReader::getDataFromIndex ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, int index, int chrg, int off )
{
	string id = gen_itoa ( fraction ) + "-" + gen_itoa ( index ) + "-1-1";
	SpecID specID ( id );
	getData ( msmsDataPointList, fraction, spectrumRange, specID, chrg, "", off );
}
void APLDataReader::getDataFromScanNumber ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const string& scanNumber, int chrg, int off )
{
	string id = gen_itoa ( fraction ) + "-" + scanNumber + "-1-1";
	SpecID specID ( id );
	getData ( msmsDataPointList, fraction, spectrumRange, specID, chrg, "", off );
}
void APLDataReader::getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const string& version, int off )
{
	if ( off != 0 ) istr.seekg ( off );
	bool readPeaks = false;
	for ( ; ; ) {
		string spot = gen_itoa ( num + 1 );
		string msmsInfo;
		DoubleVector mVector;
		IntVector zVector;
		int spectrumNumber;
		parseHeaderLines ( spot, spectrumNumber, msmsInfo, mVector, zVector );
		num++;
		if ( spectrumRange.stopReadingFile ( fraction, num ) ) return;
		if ( spectrumRange.startReadingFile ( fraction, num ) ) {
			int run = 1;
			readPeaks = specID.isMSMSSpectrumToRead100 ( fraction, spot, run, spectrumNumber );
			if ( readPeaks ) {
				msmsDataPoint.setPointInfo ( fraction, spot, run, spectrumNumber, msmsInfo );
				dataPoint = &msmsDataPoint;
				double intensity = 100.0;
				msmsDataPoint.setPrecursorIntensity ( intensity );
				bool moreData = readMSMSData ();
				if ( specID.isFilterSpectrum () ) {		// Single spectrum
					int idx = specID.getSpecNum ();
					idx /= 100;
					msmsDataPoint.setPrecursorMZ ( mVector [idx] );
					if ( chrg != 0 )	msmsDataPoint.setPrecursorCharge ( chrg );
					else				msmsDataPoint.setPrecursorCharge ( zVector [idx] );
					msmsDataPointList.push_back ( msmsDataPoint );
					return;
				}
				else {
					for ( int i = 0 ; i < mVector.size () ; i++ ) {
						msmsDataPoint.setPrecursorMZ ( mVector [i] );
						msmsDataPoint.setPrecursorCharge ( zVector [i] );
						int j = 100 * i;
						msmsDataPoint.setSpectrumNumber ( j + spectrumNumber );
						msmsDataPointList.push_back ( msmsDataPoint );
					}
				}
				readPeaks = false;
				if ( moreData == false ) return;	// End of file
			}
			else {
				if ( skipMSMSData () == false ) break;	// End of file
			}
		}
		else {
			if ( skipMSMSData () == false ) break;	// End of file
		}
	}
}
bool APLDataReader::readMSMSData ()
{
	dataPoint->clear ();
	string line;

	int c;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( islower ( c ) ) return true;				// More scans
		if ( !getline ( istr, line ) ) return false;	// No more scans
		if ( line.length () != 0 ) {
			if ( isdigit ( line [0] ) ) readLine ( const_cast <char*> (line.c_str ()) );
			else {
				line = gen_strtrim ( line );	// Try trimming the line
				if ( isdigit ( line [0] ) ) readLine ( const_cast <char*> (line.c_str ()) );
			}
		}
	}
	return false;	// No more scans
}
bool APLDataReader::skipMSMSData ()
{
	string line;
	int c;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( islower ( c ) ) return true;				// More scans
		if ( !getline ( istr, line ) ) return false;	// No more scans
	}
	return false;	// No more scans
}
void APLDataReader::writePeakList ( ostream& os )
{
	int c;
	int numS = 0;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( numS == 1 && c == 'p' ) {  //		2nd scan line, no data
			break;			// data section
		}
		string line;
		if ( !getline ( istr, line ) ) {
			break;
		}
		if ( line.length () != 0 ) {
			if ( line [line.length ()-1] == '\r' ) line = line.substr ( 0, line.length ()-1 );
			if ( line [0] == 'p' ) {
				numS++;
			}
			os << line << endl;
		}
	}
}
bool APLDataReader::writePeakList ( ostream& os, const SpecID& specID, const string& version )
{
	string sID = specID.getPureSpecID100 ();
	MapStringToStreamposConstIterator iter = mapSpecID.find ( sID );
	if ( iter != mapSpecID.end () ) {			// Scan position already found
		streampos sp = iter->second;
		if ( sp > 0 ) {
			istr.clear ();
			istr.seekg ( sp );
			writePeakList ( os );
			mapSpecID [sID] = -1;			// To indicate that the scan has already been written
		}
		return true;
	}
	istr.seekg ( endPos );					// Scan not found, go to the furthest reached so far
	for ( ; ; ) {
		string spot = gen_itoa ( num + 1 );
		string msmsInfo;
		DoubleVector mVector;
		IntVector zVector;
		int spectrumNumber;
		spos = istr.tellg ();
		parseHeaderLines ( spot, spectrumNumber, msmsInfo, mVector, zVector );
		int run = 1;
		int fraction = specID.getFraction ();
		ostringstream ost;
		ost << fraction << '-' << spot << '-' << run << '-' << spectrumNumber;
		mapSpecID [ost.str ()] = spos;
		num++;
		if ( specID.isMSMSSpectrumToRead100 ( spot, run, spectrumNumber ) ) {
			istr.seekg ( spos );
			writePeakList ( os );
			endPos = istr.tellg ();
			mapSpecID [ost.str ()] = -1;
			return true;
		}
		if ( skipMSMSData () == false ) break;	// End of file
	}
	return false;
}
MGFDataReader::MGFDataReader ( const string& s, bool fileFlag ) :
	DataReader ( s, fileFlag, true )
{
	MGFInfo::instance ().reset ();
}
void MGFDataReader::initMaps ( int fraction, int chrg, const string& version, UpdatingJavascriptMessage* ujm )
// Used by BiblioSpec creation in Search Compare
{
	string line;
	bool readRT = version == "" || !Version::isOlderVersion ( version, "5.13.1" );
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			if ( !line.compare ( 0, 5, "BEGIN" ) ) spos = istr.tellg ();
			else {
				string spot = gen_itoa ( num + 1 );
				int run = 1;
				string msmsInfo;
				if ( MGFInfo::instance ().getTitleParams ( line, spot, run, msmsInfo ) ) {
					if ( msmsInfo == "s" ) {
						msmsInfo = MGFInfo::instance ().getMSMSInfoFromScans ( getScansFromScanLine () );
					}
					if ( readRT ) getRTInSeconds ( spot );		// The spot field is filled with the RTINSECONDS line if present
					MapStringToIntIterator cur = spotNumber.find ( spot );
					if ( cur != spotNumber.end () ) cur->second++;
					else spotNumber [spot] = 1;
					cur = spotNumber.find ( spot );
					int spectrumNumber = cur->second;
					ostringstream ost;
					num++;
					if ( num % 100 == 0 ) ujm->writeMessage ( cout, "Initializing peak list map " + gen_itoa ( num ) + " spectra processed" );
					ost << fraction << '-' << spot << '-' << run << '-' << spectrumNumber;
					mapSpecID [ost.str ()] = spos;
				}
			}
		}
	}
}
void MGFDataReader::getData ( MSMSDataPointVector& msmsDataPointList, const SpecID& specID, int chrg, const string& version, UpdatingJavascriptMessage* ujm )
// Used by BiblioSpec creation in Search Compare
{
	if ( mapSpecID.empty () ) initMaps ( specID.getFraction (), chrg, version, ujm );
	string sID = specID.getPureSpecID10000 ();
	MapStringToStreamposConstIterator iter = mapSpecID.find ( sID );
	if ( iter != mapSpecID.end () ) {	// Scan position already found
		streampos sp = iter->second;
		if ( sp > 0 ) {
			istr.clear ();
			istr.seekg ( sp );
			msmsDataPoint.setPointInfo ( 1, "", 1, 1, "" );
			readMSMSData ( msmsDataPointList, chrg, 0 );
			mapSpecID [sID] = -1;			// To indicate that the scan has already been written
			return;
		}
	}
}
void MGFDataReader::getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const string& version, int off )
// Used by MS-Product
{
	bool readRT = version == "" || !Version::isOlderVersion ( version, "5.13.1" );
	lastVersion = version;
	if ( off != 0 ) istr.seekg ( off );
	string line;
	spos = istr.tellg ();
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			if ( !line.compare ( 0, 5, "BEGIN" ) ) spos = istr.tellg ();
			else {
				string spot = gen_itoa ( num + 1 );
				int run = 1;
				string msmsInfo;
				if ( MGFInfo::instance ().getTitleParams ( line, spot, run, msmsInfo ) ) {
					if ( msmsInfo == "s" ) {
						msmsInfo = MGFInfo::instance ().getMSMSInfoFromScans ( getScansFromScanLine () );
					}
					if ( readRT ) getRTInSeconds ( spot );		// The spot field is filled with the RTINSECONDS line if present
					MapStringToIntIterator cur = spotNumber.find ( spot );
					if ( cur != spotNumber.end () ) cur->second++;
					else spotNumber [spot] = 1;
					cur = spotNumber.find ( spot );
					int spectrumNumber = cur->second;
					num++;
					if ( spectrumRange.stopReadingFile ( fraction, num ) ) return;
					if ( spectrumRange.startReadingFile ( fraction, num ) ) {
						if ( readSpectrumIfAppropriate ( msmsDataPointList, specID, fraction, spot, run, spectrumNumber, msmsInfo, spos, chrg ) ) {
							return;
						}
					}
				}
			}
		}
	}
}
void MGFDataReader::getDataFromTitle ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const string& title, int chrg, int off )
{
	if ( off != 0 ) istr.seekg ( off );
	string line;
	spos = istr.tellg ();
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			if ( !line.compare ( 0, 5, "BEGIN" ) ) spos = istr.tellg ();
			else {
				if ( !line.compare ( 0, 5, "TITLE" ) ) {
					num++;
					if ( !line.substr ( 6 ).compare ( title ) ) { 
						istr.seekg ( spos );
						msmsDataPoint.setPointInfo ( 1, "", 1, 1, "" );
						readMSMSData ( msmsDataPointList, chrg, 0 );
						return;
					}
				}
			}
		}
	}
}
void MGFDataReader::getDataFromIndex ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, int index, int chrg, int off )
{
	if ( off != 0 ) istr.seekg ( off );
	string line;
	spos = istr.tellg ();
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			if ( !line.compare ( 0, 5, "BEGIN" ) ) spos = istr.tellg ();
			else {
				if ( !line.compare ( 0, 5, "TITLE" ) ) {
					num++;
					if ( num == index ) { 
						istr.seekg ( spos );
						msmsDataPoint.setPointInfo ( 1, "", 1, 1, "" );
						readMSMSData ( msmsDataPointList, chrg, 0 );
						return;
					}
				}
			}
		}
	}
}
void MGFDataReader::getDataFromMOverZ ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, double mOverZ, int chrg, int off )
{
	if ( off != 0 ) istr.seekg ( off );
	string line;
	spos = istr.tellg ();
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			if ( !line.compare ( 0, 5, "BEGIN" ) ) spos = istr.tellg ();
			else {
				if ( !line.compare ( 0, 5, "TITLE" ) ) {
					num++;
					istr.seekg ( spos );
					msmsDataPoint.setPointInfo ( 1, "", 1, 1, "" );
					if ( readMSMSData ( msmsDataPointList, chrg, 0, mOverZ ) ) return;
				}
			}
		}
	}
}
void MGFDataReader::getDataFromScanNumber ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const string& scanNumber, int chrg, int off )
{
	if ( off != 0 ) istr.seekg ( off );
	string line;
	spos = istr.tellg ();
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			if ( !line.compare ( 0, 5, "BEGIN" ) ) spos = istr.tellg ();
			else {
				string spot = gen_itoa ( num + 1 );
				int run = 1;
				string msmsInfo;
				if ( MGFInfo::instance ().getTitleParams ( line, spot, run, msmsInfo ) ) {
					if ( msmsInfo == "s" ) {
						msmsInfo = MGFInfo::instance ().getMSMSInfoFromScans ( getScansFromScanLine () );
					}
					int c1 = scanNumber.find ( ',' );
					if ( c1 == string::npos ) {				// If no comma in scanNumber
						int c2 = msmsInfo.find ( ',' );
						if ( c2 != string::npos ) {			// If comma in msmsInfo
							msmsInfo = msmsInfo.substr ( 0, c2 );
						}
					}
					if ( msmsInfo == scanNumber ) {
						istr.seekg ( spos );
						msmsDataPoint.setPointInfo ( 1, "", 1, 1, "" );
						readMSMSData ( msmsDataPointList, chrg, 0 );
						return;
					}
				}
			}
		}
	}
}
void MGFDataReader::writePeakList ( ostream& os )
{
	int c;
	bool begin = false;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( c == '>' || c == 'E' ) break;
		if ( !begin ) {
			os << "BEGIN IONS" << endl;
			begin = true;
		}
		string line;
		if ( !getline ( istr, line ) ) {
			break;
		}
		if ( line.length () != 0 ) {
			if ( line [line.length ()-1] == '\r' ) line = line.substr ( 0, line.length ()-1 );
			os << line << endl;
		}
	}
	if ( begin ) os << "END IONS" << endl;
}
bool MGFDataReader::writePeakList ( ostream& os, const SpecID& specID, const string& version )
// Called by Search Compare - printMGF
// MSProductSearch - printBodyMGFSpectrum
{
	bool readRT = version == "" || !Version::isOlderVersion ( version, "5.13.1" );
	string sID = specID.getPureSpecID10000 ();
	MapStringToStreamposConstIterator iter = mapSpecID.find ( sID );
	if ( iter != mapSpecID.end () ) {	// Scan position already found
		streampos sp = iter->second;
		if ( sp > 0 ) {
			istr.clear ();
			istr.seekg ( sp );
			writePeakList ( os );
			mapSpecID [sID] = -1;			// To indicate that the scan has already been written
		}
		return true;
	}
	istr.seekg ( endPos );					// Scan not found, go to the furthest reached so far
	string line;
	spos = istr.tellg ();
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			if ( !line.compare ( 0, 5, "BEGIN" ) ) spos = istr.tellg ();
			else {
				string spot = gen_itoa ( num + 1 );
				int run = 1;
				string msmsInfo;
				if ( MGFInfo::instance ().getTitleParams ( line, spot, run, msmsInfo ) ) {
					if ( msmsInfo == "s" ) {
						msmsInfo = MGFInfo::instance ().getMSMSInfoFromScans ( getScansFromScanLine () );
					}
					if ( readRT ) getRTInSeconds ( spot );		// The spot field is filled with the RTINSECONDS line if present
					MapStringToIntIterator cur = spotNumber.find ( spot );
					if ( cur != spotNumber.end () ) cur->second++;
					else spotNumber [spot] = 1;
					cur = spotNumber.find ( spot );
					int spectrumNumber = cur->second;
					ostringstream ost;
					ost << specID.getFraction () << '-' << spot << '-' << run << '-' << spectrumNumber;
					mapSpecID [ost.str ()] = spos;
					num++;
					if ( specID.isMSMSSpectrumToRead10000 ( spot, run, spectrumNumber ) ) {
						istr.seekg ( spos );
						writePeakList ( os );
						endPos = istr.tellg ();
						mapSpecID [ost.str ()] = -1;
						return true;
					}
				}
			}
		}
	}
	return false;
}
bool MGFDataReader::getMSMSQuanData ( XYData& xyData, const SpecID& specID, double startMass, double endMass, const string& version )
// Used for iTRAQ quantitation using centroid files
{
	bool readRT = version == "" || !Version::isOlderVersion ( version, "5.13.1" );
	MapStringToStreamposConstIterator iter = mapSpecID.find ( specID.getSpecID () );
	if ( iter != mapSpecID.end () ) {			// Scan position already found
		istr.clear ();
		istr.seekg ( iter->second );
		readMSMSQuanData ( xyData, startMass, endMass );
		return true;
	}
	istr.seekg ( endPos );					// Scan not found, go to the furthest reached so far
	string line;
	spos = istr.tellg ();
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			if ( !line.compare ( 0, 5, "BEGIN" ) ) spos = istr.tellg ();
			else {
				string spot = gen_itoa ( num + 1 );
				int run = 1;
				string msmsInfo;
				if ( MGFInfo::instance ().getTitleParams ( line, spot, run, msmsInfo ) ) {
					if ( msmsInfo == "s" ) {
						msmsInfo = MGFInfo::instance ().getMSMSInfoFromScans ( getScansFromScanLine () );
					}
					if ( readRT ) getRTInSeconds ( spot );		// The spot field is filled with the RTINSECONDS line if present
					MapStringToIntIterator cur = spotNumber.find ( spot );
					if ( cur != spotNumber.end () ) cur->second++;
					else spotNumber [spot] = 1;
					cur = spotNumber.find ( spot );
					int spectrumNumber = cur->second;
					ostringstream ost;
					ost << specID.getFraction () << '-' << spot << '-' << run << '-' << spectrumNumber;
					mapSpecID [ost.str ()] = spos;
					num++;
					if ( specID.isMSMSSpectrumToRead10000 ( spot, run, spectrumNumber ) ) {
						istr.seekg ( spos );
						readMSMSQuanData ( xyData, startMass, endMass );
						endPos = istr.tellg ();
						return true;
					}
				}
			}
		}
	}
	return false;
}
void MGFDataReader::getRTInSeconds ( string& rt )
{
	string line;
	streampos spos2 = istr.tellg ();
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 ) {
			if ( !line.compare ( 0, 12, "RTINSECONDS=" ) ) {
				rt = secToMins ( line.substr ( 12 ) );
				break;
			}
			else if ( isdigit ( line [0] ) || line == "END IONS" || isdigit ( gen_strtrim ( line ) [0] ) ) {
				break;
			}
		}
	}
	istr.seekg ( spos2 );
}
IntVector MGFDataReader::getScansFromScanLine ()
{
	streampos spos2 = istr.tellg ();
	IntVector scans;
	bool range = false;
	string line;
	while ( getline ( istr, line ) ) {
		if ( !line.compare ( 0, 6, "SCANS=" ) ) {
			string s = line.substr ( 6 );
			string::size_type start = 0;
			string::size_type end = 0;
			for ( ; ; ) {
				end = s.find_first_of ( ",-\r", start );
				int sc = atoi ( s.substr ( start, end-start ).c_str () );
				scans.push_back ( sc );
				if ( range ) {
					range = false;
					for ( int i = scans.back () + 1 ; i <= sc ; i++ ) {
						scans.push_back ( sc );
					}
				}
				if ( end == string::npos ) break;
				if ( s [end] == '-' ) range = true;
				start = end + 1;
			}
			break;
		}
	}
	istr.seekg ( spos2 );
	return scans;
}
void MGFDataReader::readParentData ( int charge, int isotopeOffset )
{
	string line;
	double mOZ = 0.0;
	double intensity = 100.0;
	streampos spos2;
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 ) {
			if ( !line.compare ( 0, 8, "PEPMASS=" ) ) {
				string s = line.substr ( 8 );
				char* ptr = const_cast <char*> (s.c_str ());
				mOZ = atof ( ptr );
				while ( *ptr && isspace ( *ptr ) ) ptr++; // Skip leading white space
				while ( *ptr && !isspace ( *ptr ) ) ptr++; // Skip mass
				while ( *ptr && isspace ( *ptr ) ) ptr++; // Skip leading white space
				if ( *ptr ) intensity = atof ( ptr );
			}
			else if ( !line.compare ( 0, 7, "CHARGE=" ) ) {
				if ( charge == 0 ) {			// Read the charge from the file
					string s = line.substr ( 7 );
					char* ptr;
					double v = strtod ( s.c_str (), &ptr );
					int z = static_cast <int> ( v );
					char c = *ptr;
					if ( c == '+' ) ptr++;
					else if ( c == '-' ) {
						z = -z;
						ptr++;
					}
					while ( *ptr && isspace ( *ptr ) ) ptr++;	// Skip white space
					if ( *ptr == 0 ) {							// Single charge
						charge = z;
					}
					else {										// Multiple charges
						IntVector chargeRange;
						chargeRange.push_back ( z );
						do {
							char c = *ptr;
							if ( c == ',' || c == 'a' || isdigit ( c ) ) {
								if ( c == ',' ) ptr++;
								if ( c == 'a' ) ptr += 3;
								v = strtod ( ptr, &ptr );
								z = static_cast <int> ( v );
								c = *ptr;
								if ( c == '+' ) ptr++;
								else if ( c == '-' ) {
									z = -z;
									ptr++;
								}
							}
							while ( *ptr && isspace ( *ptr ) ) ptr++; // Skip white space
							chargeRange.push_back ( z );
						} while ( *ptr );
						msmsDataPoint.setPrecursorChargeRange ( chargeRange );
					}
				}
			}
			else if ( isdigit ( line [0] ) || line == "END IONS" || isdigit ( gen_strtrim ( line ) [0] ) ) {
				istr.seekg ( spos2 );
				if ( mOZ == 0.0 ) ErrorHandler::genError ()->error ( "m/z not specified in the data file.\n" );
				//if ( charge == 0 ) ErrorHandler::genError ()->error ( "Charge not specified in the data file.\n" );
				if ( isotopeOffset != 0 ) mOZ += ( isotopeOffset * 1.0029 / charge );
				msmsDataPoint.setPrecursor ( mOZ, charge, intensity );
				return;
			}
		}
		spos2 = istr.tellg ();
	}
}
void MGFDataReader::skipParentData ()
{
	string line;
	streampos spos2;
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 ) {
			if ( isdigit ( line [0] ) || line == "END IONS" || isdigit ( gen_strtrim ( line ) [0] ) ) {
				istr.seekg ( spos2 );
				return;
			}
		}
		spos2 = istr.tellg ();
	}
}
bool MGFDataReader::readSpectrumIfAppropriate ( MSMSDataPointVector& msmsDataPointList, const SpecID& specID, int fraction, const string& spot, int run, int spectrumNumber, const string& msmsInfo, streampos& spos2, int chrg )
{
	if ( specID.isMSMSSpectrumToRead10000 ( fraction, spot, run, spectrumNumber ) ) {
		istr.seekg ( spos2 );
		msmsDataPoint.setPointInfo ( fraction, spot, run, ( specID.getFraction () == -1 ) ? spectrumNumber : specID.getSpecNum (), msmsInfo );
		readMSMSData ( msmsDataPointList, chrg, 0 );
		if ( specID.isFilterSpectrum () ) return true;
	}
	return false;
}
PPDataReader::PPDataReader ( const string& s, bool fileFlag ) :
	DataReader ( s, fileFlag, true )
{
}
void PPDataReader::getData ( MSDataPointVector& msDataPointList, int fraction, const SpecID& specID )
{
	string line;
	bool readPeaks = false;		// Read peaks as well as title information
	IntVector possMSRuns;
	spos = istr.tellg ();
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			string spot;
			int run;
			if ( getMSTitleParams ( line, spot, run ) ) {
				readPeaks = specID.isMSSpectrumToRead ( fraction, spot, run );
				if ( readPeaks ) {
					msDataPoint.setPointInfo ( fraction, spot, run );
					readMSData ( msDataPointList );
					readPeaks = false;
					if ( specID.isFilterSpectrum () ) return;
				}
				else {
					if ( specID.isMSSpectrumToRead ( fraction, spot ) ) {
						possMSRuns.push_back ( run );
					}
				}
			}
		}
	}
#ifdef BATCHTAG
	if ( !possMSRuns.empty () ) {		// Deals with the case where the run number is different
		int idx = getT2DMSRunNumberIndex ( specID.getRun (), possMSRuns );
		SpecID sid = specID;
		sid.setRun ( possMSRuns [idx] );
		istr.clear ();
		istr.seekg ( spos );
		getData ( msDataPointList, fraction, sid );
	}
#endif
}
void PPDataReader::getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const string& version, int off )
{
	int num = 0;
	string line;
	bool readPeaks = false;		// Read peaks as well as title information
	while ( getline ( istr, line ) ) {
		if ( line.length () != 0 && line [0] != '#' ) {
			string spot;
			int run;
			int spectrumNumber;
			string msmsInfo;
			if ( getMSMSTitleParams ( line, spot, run, spectrumNumber, msmsInfo ) ) {
				num++;
				if ( spectrumRange.stopReadingFile ( fraction, num ) ) return;
				if ( spectrumRange.startReadingFile ( fraction, num ) ) {
					readPeaks = specID.isMSMSSpectrumToRead10000 ( fraction, spot, run, spectrumNumber );
					if ( readPeaks ) {
						msmsDataPoint.setPointInfo ( fraction, spot, run, spectrumNumber, msmsInfo );
						readMSMSData ( msmsDataPointList, chrg, 0 );	// Charge set in file
						readPeaks = false;
						if ( specID.isFilterSpectrum () ) return;
					}
				}
			}
		}
	}
}
string PPDataReader::getParameter ( const string& line, const string& param )
{
	int start = line.find ( param );
	if ( start == string::npos ) {
		return "";
	}
	else {
		start += param.length ();
		int end = line.find ( "$", start );
		return line.substr ( start, end - start );
	}
}

int PPDataReader::getIntValue ( const string& str, const string& tag )
{
	return atoi ( getParameter ( str, tag ).c_str () );
}
double PPDataReader::getDoubleValue ( const string& str, const string& tag )
{
	return atof ( getParameter ( str, tag ).c_str () );
}
string PPDataReader::getStringValue ( const string& str, const string& tag )
{
	return getParameter ( str, tag );
}
PPDataReader2::PPDataReader2 ( const string& s, bool fileFlag, const string& dataFormat ) :
	DataReader ( s, fileFlag, dataFormat == "PP M/Z Intensity Charge" )
{
}
void PPDataReader2::getData ( MSDataPointVector& msDataPointList, int fraction, const SpecID& specID )
{
	string line;
	int c;
	int i = 1;
	while ( ( c = istr.peek () ) != EOF ) {
		if ( c == '>' ) getline ( istr, line );
		if ( c == 'E' ) {		// This shouldn't really happen but did once.
			ErrorHandler::genError ()->error ( "The data is in an illegal format.\n" );
		}
		string s = gen_itoa ( i );
		if ( specID.isMSSpectrumToRead ( fraction, s, 1 ) ) {
			msDataPoint.setPointInfo ( fraction, s, 1 );
			readMSData ( msDataPointList );
			if ( specID.isFilterSpectrum () ) return;
		}
		else {
			while ( ( c = istr.peek () ) != EOF ) {		// Skip this spectrum
				if ( c == '>' || c == 'E' ) break;
				if ( !getline ( istr, line ) ) break;
			}
		}
		i++;
	}
}
void PPDataReader2::getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const string& version, int off )
{
	string line;
	int i = 1;
	while ( istr.peek () != EOF ) {
		msmsDataPoint.setPointInfo ( fraction, gen_itoa ( i ), 1, 1, "" );
		readMSMSData ( msmsDataPointList, chrg, 0 );	// Charge set in file
		i++;
	}
}
SpaceSeparatedSpectraDataReader::SpaceSeparatedSpectraDataReader ( const string& s, bool fileFlag ) :
	DataReader ( s, fileFlag, true )
{
	num = 1;
}
void SpaceSeparatedSpectraDataReader::readData ()	// OK
{
	dataPoint->clear ();
	string line;
	int c;
	while ( ( c = istr.peek () ) != EOF ) {	// Keep going until you come to a line with only blanks or the end of file.
		if ( !getline ( istr, line ) ) return;
		if ( genEmptyString ( line ) ) return;
		if ( isdigit ( line [0] ) ) readLine ( const_cast <char*> (line.c_str ()) );
	}
}
void SpaceSeparatedSpectraDataReader::getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const string& version, int off )
{
	string line;
	spos = istr.tellg ();
	bool readPeaks = false;
	while ( getline ( istr, line ) ) {
		if ( !genEmptyString ( line ) ) {
			if ( spectrumRange.stopReadingFile ( fraction, num ) ) return;
			if ( spectrumRange.startReadingFile ( fraction, num ) ) {
				string spot = gen_itoa ( num );
				readPeaks = specID.isMSMSSpectrumToRead10000 ( fraction, spot, 1, 1 );
				if ( readPeaks ) {
					istr.seekg ( spos );
					msmsDataPoint.setPointInfo ( fraction, spot, 1, 1, "" );
					readMSMSData ( msmsDataPointList, chrg, 0 );	// Charge set in file
					readPeaks = false;
					if ( specID.isFilterSpectrum () ) return;
				}
			}
			if ( !readPeaks ) {
				while ( getline ( istr, line ) ) {
					if ( genEmptyString ( line ) ) break;
				}
			}
			num++;
		}
		spos = istr.tellg ();
	}
}
bool SpaceSeparatedSpectraDataReader::getMSMSQuanData ( XYData& xyData, const SpecID& specID, double startMass, double endMass, const string& version )
{
	MapStringToStreamposConstIterator iter = mapSpecID.find ( specID.getSpecID () );
	if ( iter != mapSpecID.end () ) {			// Scan position already found
		istr.clear ();
		istr.seekg ( iter->second );
		readMSMSQuanData ( xyData, startMass, endMass );
		return true;
	}
	istr.seekg ( endPos );					// Scan not found, go to the furthest reached so far
	string line;
	spos = istr.tellg ();
	while ( getline ( istr, line ) ) {
		if ( !genEmptyString ( line ) ) {
			string spot = gen_itoa ( num );
			ostringstream ost;
			ost << specID.getFraction () << '-' << spot << "-1-1";
			mapSpecID [ost.str ()] = spos;
			num++;
			if ( specID.isMSMSSpectrumToRead10000 ( spot, 1, 1 ) ) {
				istr.seekg ( spos );
				readMSMSQuanData ( xyData, startMass, endMass );
				endPos = istr.tellg ();
				return true;
			}
			else {
				while ( getline ( istr, line ) ) {
					if ( genEmptyString ( line ) ) break;
				}
			}
		}
		spos = istr.tellg ();
	}
	return false;
}
PKLDataReader::PKLDataReader ( const string& s, bool fileFlag ) :
	SpaceSeparatedSpectraDataReader ( s, fileFlag )
{
}
void PKLDataReader::readParentData ( int charge, int isotopeOffset )
{
	string line;
	while ( getline ( istr, line ) ) {
		if ( !genEmptyString ( line ) ) {		// Skip any more blank lines
			istringstream ist ( line );
			double mOZ;
			double intensity;
			if ( ist >> mOZ ) {
				if ( ist >> intensity ) {
					if ( charge == 0 ) {			// Read the charge from the file
						if ( ( ist >> charge ) == 0 ) charge = 1;
					}
					if ( isotopeOffset != 0 ) mOZ += ( isotopeOffset * 1.0029 / charge );
					msmsDataPoint.setPrecursor ( mOZ, charge, intensity );
					return;
				}
			}
		}
	}
}
void PKLDataReader::skipParentData ()
{
	string line;
	while ( getline ( istr, line ) ) {
		if ( !genEmptyString ( line ) ) {		// Skip any more blank lines
			istringstream ist ( line );
			double mOZ;
			double intensity;
			if ( ist >> mOZ ) {
				if ( ist >> intensity ) {
					return;
				}
			}
		}
	}
}
DTADataReader::DTADataReader ( const string& s, bool fileFlag ) :
	SpaceSeparatedSpectraDataReader ( s, fileFlag )
{
}
void DTADataReader::readParentData ( int charge, int isotopeOffset )
{
	string line;
	while ( getline ( istr, line ) ) {
		if ( !genEmptyString ( line ) ) {		// Skip any more blank lines
			istringstream ist ( line );
			double mhPlus;
			if ( ist >> mhPlus ) {
				if ( charge == 0 ) {			// Read the charge from the file
					if ( ( ist >> charge ) == 0 ) charge = 1;
				}
				double mOZ = mPlusHToMOverZ ( mhPlus, charge, true );
				if ( isotopeOffset != 0 ) mOZ += ( isotopeOffset * 1.0029 / charge );
				msmsDataPoint.setPrecursor ( mOZ, charge, 100.0 );
				return;
			}
		}
	}
}
void DTADataReader::skipParentData ()
{
	string line;
	while ( getline ( istr, line ) ) {
		if ( !genEmptyString ( line ) ) {		// Skip any more blank lines
			istringstream ist ( line );
			double mhPlus;
			if ( ist >> mhPlus ) {
				return;
			}
		}
	}
}
