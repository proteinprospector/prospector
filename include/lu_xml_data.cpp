/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_xml_data.cpp                                               *
*                                                                             *
*  Created    : July 12th 2005                                                *
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
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#include <zlib.h>
#include <lg_io.h>
#include <lg_endian.h>
#include <lg_base64.h>
#include <lg_string.h>
#include <lu_aa_info.h>
#include <lu_mass_conv.h>
#include <lu_mass_elem.h>
#include <lu_spec_id.h>
#include <lu_xml_data.h>
using std::string;
using std::vector;
using std::ostringstream;

PPExpatSpecData::PPExpatSpecData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange ) :
	msmsDataPointList ( msmsDataPointList ),
	spectrumRange ( spectrumRange ),
	specID ( specID ),
	chrg ( chrg ),
	chargeRange ( chargeRange ),
	fraction ( fraction ),
	count ( 0 )
{
}
PPExpatSpecData::~PPExpatSpecData () {}
void PPExpatSpecData::makeSpectrum ()
{
	MapStringToIntIterator cur = spotNumber.find ( rt );
	if ( cur != spotNumber.end () ) cur->second++;
	else spotNumber [rt] = 1;
	cur = spotNumber.find ( rt );
	int spectrumNumber = cur->second;
	if ( specID.isMSMSSpectrumToRead10000 ( fraction, rt, 1, spectrumNumber ) ) {
		msmsDataPoint.setPointInfo ( fraction, rt, 1, ( specID.getFraction () == -1 ) ? spectrumNumber : specID.getSpecNum (), scan );
		int z;
		if ( chrg == -9999 )
			z = ( precursorCharge == 0 ) ? chrg : precursorCharge;
		else
			z = ( chrg == 0 ) ? precursorCharge : chrg;

		msmsDataPoint.setPrecursor ( precursorMz, z, precursorIntensity );
		if ( z == 0 ) {
			for ( IntVectorSizeType i = 0 ; i < chargeRange.size () ; i++ ) {
				msmsDataPoint.setPrecursorCharge ( chargeRange [i] );
				int j = 10000 * chargeRange [i];
				msmsDataPoint.setSpectrumNumber ( j + spectrumNumber );
				msmsDataPointList.push_back ( msmsDataPoint );
			}
		}
		else
			msmsDataPointList.push_back ( msmsDataPoint );
	}
}

static void makeScan32 ( const string& s1, bool bigEndian, MSMSDataPoint& msmsDataPoint )
{
	string s2;
	base64Decode ( s1, s2 );
	int siz = s2.length () / sizeof (float);
	unsigned int* ui = const_cast <unsigned int*> (reinterpret_cast<const unsigned int*> (s2.data ()));
	if ( bigEndian ) for ( int i = 0 ; i < siz ; i++ ) genEndianConvert ( ui [i] );
	msmsDataPoint.clear ();
	for ( int i = 0 ; i < siz ; ) {
		double mz = reinterpret_cast<float*> (ui) [i++];
		double inten = reinterpret_cast<float*> (ui) [i++];
		msmsDataPoint.addPeak ( mz, 1, inten );
	}
}
static void makeScan32Uncompress ( const string& s1, int numPeaks, bool bigEndian, MSMSDataPoint& msmsDataPoint )
{
	string s2;
	base64Decode ( s1, s2 );
	int numSrcBytes = s2.length ();
	size_t bufferLen = numPeaks*2;		// mz and intensity
	unsigned int* buffer = new unsigned int [bufferLen];
	unsigned long destLen = bufferLen * sizeof (unsigned int);
	const unsigned char* src = const_cast <const unsigned char*> (reinterpret_cast<const unsigned char*> (s2.data ()));
	uncompress ( (unsigned char*)&buffer[0], &destLen, src, numSrcBytes );
	if ( bigEndian ) for ( int i = 0 ; i < bufferLen ; i++ ) genEndianConvert ( buffer [i] );
	msmsDataPoint.clear ();
	for ( int i = 0 ; i < bufferLen ; ) {
		double mz = reinterpret_cast<float*> (buffer) [i++];
		double inten = reinterpret_cast<float*> (buffer) [i++];
		msmsDataPoint.addPeak ( mz, 1, inten );
	}
	delete [] buffer;
}
static void makeScan64 ( const string& s1, bool bigEndian, MSMSDataPoint& msmsDataPoint )
{
	string s2;
	base64Decode ( s1, s2 );
	int siz = s2.length () / sizeof (double);
	GENUINT64* ui = const_cast <GENUINT64*> (reinterpret_cast<const GENUINT64*> (s2.data ()));
	if ( bigEndian ) for ( int i = 0 ; i < siz ; i++ ) genEndianConvert ( ui [i] );
	msmsDataPoint.clear ();
	for ( int i = 0 ; i < siz ; ) {
		double mz = reinterpret_cast<double*> (ui) [i++];
		double inten = reinterpret_cast<double*> (ui) [i++];
		msmsDataPoint.addPeak ( mz, 1, inten );
	}
}
static void makeScan64Uncompress ( const string& s1, int numPeaks, bool bigEndian, MSMSDataPoint& msmsDataPoint )
{
	string s2;
	base64Decode ( s1, s2 );
	int numSrcBytes = s2.length ();
	size_t bufferLen = numPeaks*2;		// mz and intensity
	GENUINT64* buffer = new GENUINT64 [bufferLen];
	unsigned long destLen = bufferLen * sizeof (GENUINT64);
	const unsigned char* src = const_cast <const unsigned char*> (reinterpret_cast<const unsigned char*> (s2.data ()));
	uncompress ( (unsigned char*)&buffer[0], &destLen, src, numSrcBytes );
	if ( bigEndian ) for ( int i = 0 ; i < bufferLen ; i++ ) genEndianConvert ( buffer [i] );
	msmsDataPoint.clear ();
	for ( int i = 0 ; i < bufferLen ; ) {
		double mz = reinterpret_cast<double*> (buffer) [i++];
		double inten = reinterpret_cast<double*> (buffer) [i++];
		msmsDataPoint.addPeak ( mz, 1, inten );
	}
	delete [] buffer;
}
template<class T1>
static void makeVector32 ( const string& s1, bool bigEndian, vector <T1>& v )
{
	string s2;
	base64Decode ( s1, s2 );
	int siz = s2.length () / sizeof (float);
	unsigned int* ui = const_cast <unsigned int*> (reinterpret_cast<const unsigned int*> (s2.data ()));
	if ( bigEndian ) for ( int i = 0 ; i < siz ; i++ ) genEndianConvert ( ui [i] );
	for ( int i = 0 ; i < siz ; i++ ) {
		v.push_back ( reinterpret_cast<float*> (ui) [i] );
	}
}
template<class T1>
static void makeVector32Uncompress ( const string& s1, int numPeaks, bool bigEndian, vector <T1>& v )
{
	string s2;
	base64Decode ( s1, s2 );
	int numSrcBytes = s2.length ();
	if ( !numPeaks ) {							// Buffer length unknown
		numPeaks = numSrcBytes / sizeof (unsigned int);
		if ( numSrcBytes % sizeof (unsigned int) ) numPeaks++;
		numPeaks *= 2;
	}
	unsigned int* buffer = new unsigned int [numPeaks];
	unsigned long destLen = numPeaks * sizeof (unsigned int);
	const unsigned char* src = const_cast <const unsigned char*> (reinterpret_cast<const unsigned char*> (s2.data ()));
	uncompress ( (unsigned char*)&buffer[0], &destLen, src, numSrcBytes );
	numPeaks = destLen / sizeof (unsigned int);		// Actual numPeaks
	if ( bigEndian ) for ( int i = 0 ; i < numPeaks ; i++ ) genEndianConvert ( buffer [i] );
	for ( int i = 0 ; i < numPeaks ; i++ ) {
		v.push_back ( reinterpret_cast<float*> (buffer) [i] );
	}
	delete [] buffer;
}
template<class T1>
static void makeVector64 ( const string& s1, bool bigEndian, vector <T1>& v )
{
	string s2;
	base64Decode ( s1, s2 );
	int siz = s2.length () / sizeof (double);
	GENUINT64* ui = const_cast <GENUINT64*> (reinterpret_cast<const GENUINT64*> (s2.data ()));
	if ( bigEndian ) for ( int i = 0 ; i < siz ; i++ ) genEndianConvert ( ui [i] );
	for ( int i = 0 ; i < siz ; i++ ) {
		v.push_back ( reinterpret_cast<double*> (ui) [i] );
	}
}
template<class T1>
static void makeVector64Uncompress ( const string& s1, int numPeaks, bool bigEndian, vector <T1>& v )
{
	string s2;
	base64Decode ( s1, s2 );
	int numSrcBytes = s2.length ();
	if ( !numPeaks ) {							// Buffer length unknown
		numPeaks = numSrcBytes / sizeof (GENUINT64);
		if ( numSrcBytes % sizeof (GENUINT64) ) numPeaks++;
		numPeaks *= 2;
	}
	GENUINT64* buffer = new GENUINT64 [numPeaks];
	unsigned long destLen = numPeaks * sizeof (GENUINT64);
	const unsigned char* src = const_cast <const unsigned char*> (reinterpret_cast<const unsigned char*> (s2.data ()));
	uncompress ( (unsigned char*)&buffer[0], &destLen, src, numSrcBytes );
	numPeaks = destLen / sizeof (GENUINT64);		// Actual numPeaks
	if ( bigEndian ) for ( int i = 0 ; i < numPeaks ; i++ ) genEndianConvert ( buffer [i] );
	for ( int i = 0 ; i < numPeaks ; i++ ) {
		v.push_back ( reinterpret_cast<double*> (buffer) [i] );
	}
	delete [] buffer;
}
/*
mzXML version 2.0
-----------------
// MS scan
<scan
num="1"
msLevel="1"						// 1=MS, 2=MSMS
peaksCount="1807"
polarity="+"
scanType="Full"
retentionTime="PT0.2728S"       // ISO 8601 time period P=period T= date/time separator 0.2728S=0.2728 seconds
lowMz="300"
highMz="2000"
basePeakMz="370.901"
basePeakIntensity="28037.1"
totIonCurrent="852330">
	<peaks
	precision="32"
	byteOrder="network"
	pairOrder="m/z-int">
	Q5YxWURD5.......
	</peaks>
</scan>


// MSMS scan
<scan
num="451"
msLevel="2"
peaksCount="214"
polarity="+"
scanType="Full"
retentionTime="PT270.058S"
collisionEnergy="28"
lowMz="180"
highMz="1405"
basePeakMz="649.53"
basePeakIntensity="145.219"
totIonCurrent="1392.06">

  <precursorMz precursorIntensity="31354.2">696.0718384</precursorMz> 
  <peaks precision="32" byteOrder="network" pairOrder="m/z-int">Q03k....</peaks> 
</scan>
*/
/*
mzXML version 3.1
-----------------
// MS scan
<scan num="4"
	msLevel="1"
	peaksCount="104"
	polarity="+"
	scanType="Full"
	filterLine="FTMS + p ESI Full ms [400.00-2000.00]"
	retentionTime="PT2.6489S"
	lowMz="401.215"
	highMz="1946.74"
	basePeakMz="445.12"
	basePeakIntensity="16327.7"
	totIonCurrent="134848" >
   
	<peaks precision="32" byteOrder="network" contentType="m/z-int" compressionType="none" compressedLen="0" >Q8i...</peaks>
// MSMS scan
<scan
	num="5"
	msLevel="2"
	peaksCount="67"
	polarity="+"
	scanType="Full"
	filterLine="ITMS + c ESI d Full ms2 429.09@cid35.00 [105.00-440.00]"
	retentionTime="PT3.5119S"
	lowMz="147.172"
	highMz="420.017"
	basePeakMz="341.016"
	basePeakIntensity="3592.79"
	totIonCurrent="5212.12"
	collisionEnergy="35" >

	<precursorMz precursorIntensity="10609.2" precursorCharge="1" activationMethod="CID" >429.0887146</precursorMz>
	<peaks precision="32" byteOrder="network" contentType="m/z-int" compressionType="none" compressedLen="0" >QxMs...</peaks>
</scan>
*/
PPExpatMZXMLData::PPExpatMZXMLData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange ) :
	PPExpatSpecData ( msmsDataPointList, fraction, spectrumRange, specID, chrg, chargeRange ),
	peaksFlag ( false ),
	precursorMzFlag ( false )
{
	bigEndian = true;	// Currently all mzXML files are big-endian
	precision = 32;
	msLevel = 0;
	peaksCount = 0;
	compression = false;
}
PPExpatMZXMLData::~PPExpatMZXMLData () {}
void PPExpatMZXMLData::startElement ( const char* name, const char** attributes )
{
	if ( !strcmp ( name, "scan" ) ) {
		getAttributeValue ( attributes, "num", scan );
		getAttributeValue ( attributes, "msLevel", msLevel );
		getAttributeValue ( attributes, "peaksCount", peaksCount );
		if ( getAttributeValue ( attributes, "retentionTime", rt ) ) {
			string seconds = rt.substr ( 2, rt.length () - 3 );
			double dseconds = atof ( seconds.c_str () );
			rt = gen_ftoa ( dseconds / 60, "%.4f" );
		}
		else rt = scan;
		if ( msLevel == 2 ) count++;
	}
	else {
		if ( msLevel == 2 ) {
			if ( !spectrumRange.stopReadingFile ( fraction, count ) ) {
				if ( spectrumRange.startReadingFile ( fraction, count ) ) {
					if ( !strcmp ( name, "peaks" ) ) {
						getAttributeValue ( attributes, "precision", precision );
						peaksFlag = true;
						string compressionType;
						getAttributeValue ( attributes, "compressionType", compressionType );
						compression = (compressionType == "zlib");
					}
					else if ( !strcmp ( name, "precursorMz" ) ) {
						precursorMzFlag = true;
						precursorCharge = 0;
						precursorIntensity = 100.0;
						getAttributeValue ( attributes, "precursorIntensity", precursorIntensity );
						getAttributeValue ( attributes, "precursorCharge", precursorCharge );
					}
				}
			}
		}
	}
}
void PPExpatMZXMLData::endElement ( const char* name )
{
	if ( msLevel == 2 ) {
		if ( !strcmp ( name, "peaks" ) ) {
			if ( precision == 32 ) {
				if ( compression )	makeScan32Uncompress ( s, peaksCount, bigEndian, msmsDataPoint );	
				else				makeScan32 ( s, bigEndian, msmsDataPoint );
			}
			else {
				if ( compression )	makeScan64Uncompress ( s, peaksCount, bigEndian, msmsDataPoint );
				else				makeScan64 ( s, bigEndian, msmsDataPoint );
			}
			peaksFlag = false;
			s = "";
		}
		else if ( !strcmp ( name, "scan" ) ) {
			if ( !spectrumRange.stopReadingFile ( fraction, count ) ) {
				if ( spectrumRange.startReadingFile ( fraction, count ) ) {
					makeSpectrum ();
				}
			}
			msLevel = 1;
		}
		else if ( !strcmp ( name, "precursorMz" ) ) {
			precursorMzFlag = false;
			precursorMz = atof ( s.c_str () );
			s = "";
		}
	}
}
void PPExpatMZXMLData::characterDataHandler ( const char* str, int len )
{
	if ( peaksFlag || precursorMzFlag ) {
		if ( msLevel == 2 ) {
			s.append ( str, len );
		}
	}
}
/*
mzML version 1.1
----------------

http://www.peptideatlas.org/tmp/mzML1.1.0.html

// MS scan
<spectrum index="0" id="scan=19" defaultArrayLength="15">
	<referenceableParamGroupRef ref="CommonMS1SpectrumParams"/>
	<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1"/>
	<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" value=""/>
	<cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="400.38999999999999" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
	<cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="1795.5599999999999" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
	<cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="445.34699999999998" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
	<cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="120053" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of counts"/>
	<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="16675500"/>
	<scanList count="1">
		<cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>
		<scan instrumentConfigurationRef="LCQ_x0020_Deca">
			<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="5.8905000000000003" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
			<cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="+ c NSI Full ms [ 400.00-1800.00]"/>
			<cvParam cvRef="MS" accession="MS:1000616" name="preset scan configuration" value="3"/>
			<scanWindowList count="1">
				<scanWindow>
					<cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="400" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
					<cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="1800" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
				</scanWindow>
			</scanWindowList>
		</scan>
	</scanList>
	<binaryDataArrayList count="2">
		<binaryDataArray encodedLength="160" dataProcessingRef="CompassXtract_x0020_processing">
              <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
              <cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
              <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value="" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
              <binary>AAAAAAAAAAAAAAAAAADwPwAAAAAAAABAAAAAAAAACEAAAAAAAAAQQAAAAAAAABRAAAAAAAAAGEAAAAAAAAAcQAAAAAAAACBAAAAAAAAAIkAAAAAAAAAkQAAAAAAAACZAAAAAAAAAKEAAAAAAAAAqQAAAAAAAACxA</binary>
		</binaryDataArray>
		<binaryDataArray encodedLength="160" dataProcessingRef="CompassXtract_x0020_processing">
			<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
			<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
			<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of counts"/>
			<binary>AAAAAAAALkAAAAAAAAAsQAAAAAAAACpAAAAAAAAAKEAAAAAAAAAmQAAAAAAAACRAAAAAAAAAIkAAAAAAAAAgQAAAAAAAABxAAAAAAAAAGEAAAAAAAAAUQAAAAAAAABBAAAAAAAAACEAAAAAAAAAAQAAAAAAAAPA/</binary>
		</binaryDataArray>
	</binaryDataArrayList>
</spectrum>

// MSMS scan
<spectrum index="1" id="scan=20" defaultArrayLength="10">
	<referenceableParamGroupRef ref="CommonMS2SpectrumParams"/>
	<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2"/>
	<cvParam cvRef="MS" accession="MS:1000128" name="profile spectrum" value=""/>
	<cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="320.38999999999999" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
	<cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="1003.5599999999999" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
	<cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="456.34699999999998" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
	<cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="23433" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of counts"/>
	<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="16675500"/>
	<scanList count="1">
		<cvParam cvRef="MS" accession="MS:1000795" name="no combination" value=""/>
		<scan instrumentConfigurationRef="LCQ_x0020_Deca">
			<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="5.9904999999999999" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>
			<cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="+ c d Full ms2  445.35@cid35.00 [ 110.00-905.00]"/>
			<cvParam cvRef="MS" accession="MS:1000616" name="preset scan configuration" value="4"/>
			<scanWindowList count="1">
				<scanWindow>
					<cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="110" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
					<cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="905" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
				</scanWindow>
			</scanWindowList>
		</scan>
	</scanList>
	<precursorList count="1">
		<precursor spectrumRef="scan=19">
			<isolationWindow>
				<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="445.30000000000001" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
				<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="0.5" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
				<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.5" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
			</isolationWindow>
			<selectedIonList count="1">
				<selectedIon>
					<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="445.33999999999997" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
					<cvParam cvRef="MS" accession="MS:1000042" name="peak intensity" value="120053"/>
					<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2"/>
				</selectedIon>
			</selectedIonList>
			<activation>
				<cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" value=""/>
				<cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35" unitCvRef="UO" unitAccession="UO:0000266" unitName="electronvolt"/>
			</activation>
		</precursor>
	</precursorList>
	<binaryDataArrayList count="2">
		<binaryDataArray encodedLength="108" dataProcessingRef="CompassXtract_x0020_processing">
			<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
			<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
			<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" value="" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>
			<binary>AAAAAAAAAAAAAAAAAAAAQAAAAAAAABBAAAAAAAAAGEAAAAAAAAAgQAAAAAAAACRAAAAAAAAAKEAAAAAAAAAsQAAAAAAAADBAAAAAAAAAMkA=</binary>
		</binaryDataArray>
		<binaryDataArray encodedLength="108" dataProcessingRef="CompassXtract_x0020_processing">
			<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" value=""/>
			<cvParam cvRef="MS" accession="MS:1000576" name="no compression" value=""/>
			<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" value="" unitCvRef="MS" unitAccession="MS:1000131" unitName="number of counts"/>
			<binary>AAAAAAAANEAAAAAAAAAyQAAAAAAAADBAAAAAAAAALEAAAAAAAAAoQAAAAAAAACRAAAAAAAAAIEAAAAAAAAAYQAAAAAAAABBAAAAAAAAAAEA=</binary>
		</binaryDataArray>
	</binaryDataArrayList>
</spectrum>
*/
PPExpatMZMLData::PPExpatMZMLData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange ) :
	PPExpatSpecData ( msmsDataPointList, fraction, spectrumRange, specID, chrg, chargeRange ),
	scanFlag ( false ),
	paramGroupFlag ( false )
{
	bigEndian = false;	// All mzML files are little-endian
	msLevel = 0;
	selectedIonFlag = false;
	compression = false;
	mzArrayFlag = false;
	intenArrayFlag = false;
	precursorMz = 0.0;
	precursorIntensity = 100.0;
	precursorCharge = 0;
	scFlag = false;
	binaryFlag = false;
	binaryDataArray = false;
}
PPExpatMZMLData::~PPExpatMZMLData () {}
void PPExpatMZMLData::startElement ( const char* name, const char** attributes )
{
	if ( !scanFlag ) {
		if ( paramGroupFlag ) {
			if ( !strcmp ( name, "cvParam" ) ) {
				string n;
				string v;
				getAttributeValue ( attributes, "name", n );
				getAttributeValue ( attributes, "value", v );
				curParam [n] = v;
			}
		}
		else {
			if ( !strcmp ( name, "spectrum" ) ) {
				string id;
				getAttributeValue ( attributes, "id", id );
				int start = id.find ( "scan=" );
				if ( start != string::npos ) {
					scan = id.substr ( start + 5 );
				}
				else {
					scan = id;
				}
				scanFlag = true;
			}
			else if ( !strcmp ( name, "referenceableParamGroup" ) ) {
				getAttributeValue ( attributes, "id", paramGroupID );
				paramGroupFlag = true;
			}
		}
	}
	else {
		if ( msLevel == 0 ) {
			if ( !strcmp ( name, "referenceableParamGroupRef" ) ) {
				string ref;
				getAttributeValue ( attributes, "ref", ref );
				std::map <std::string, MapStringToString>::const_iterator cur = paramGroups.find ( ref );
				if ( cur != paramGroups.end () ) {
					MapStringToStringConstIterator cur2 = (*cur).second.find ( "ms level" );
					if ( cur2 != (*cur).second.end () ) msLevel = atoi ( (*cur2).second.c_str () );
				}
			}
			else {
				getCVAttributeValue ( name, attributes, "ms level", msLevel );
			}
			if ( msLevel == 2 ) count++;
		}
		else if ( msLevel == 2 ) {
			if ( !spectrumRange.stopReadingFile ( fraction, count ) ) {
				if ( spectrumRange.startReadingFile ( fraction, count ) ) {
					if ( !scFlag ) {
						if ( !strcmp ( name, "scan" ) ) scFlag = true;
					}
					else {
						string v;
						getCVAttributeValue ( name, attributes, "scan time", v );			// 1.0
						getCVAttributeValue ( name, attributes, "scan start time", v );		// 1.1
						if ( !v.empty () ) {
							rt = gen_ftoa ( atof ( v.c_str () ), "%.4f" );
						}
					}
					if ( !selectedIonFlag ) {
						if ( !strcmp ( name, "selectedIon" ) ) selectedIonFlag = true;
					}
					else {
						getCVAttributeValue ( name, attributes, "m/z", precursorMz );					// 1.0
						getCVAttributeValue ( name, attributes, "selected ion m/z", precursorMz );		// 1.1
						getCVAttributeValue ( name, attributes, "intensity", precursorIntensity );		// 1.0
						getCVAttributeValue ( name, attributes, "peak intensity", precursorIntensity );	// 1.1
						getCVAttributeValue ( name, attributes, "charge state", precursorCharge );
					}
					if ( !binaryDataArray ) {
						if ( !strcmp ( name, "binaryDataArray" ) ) {
							binaryDataArray = true;
						}
					}
					else {
						string dummy;
						if ( getCVAttributeName ( name, attributes, "zlib compression" ) ) {
							compression = true;
						}
						if ( getCVAttributeName ( name, attributes, "no compression" ) ) {
							compression = false;
						}
						if ( getCVAttributeName ( name, attributes, "m/z array" ) ) {
							mzArrayFlag = true;
						}
						if ( getCVAttributeName ( name, attributes, "intensity array" ) ) {
							intenArrayFlag = true;
						}
						if ( getCVAttributeName ( name, attributes, "32-bit float" ) ) {
							precision = 32;
						}
						if ( getCVAttributeName ( name, attributes, "64-bit float" ) ) {
							precision = 64;
						}
						if ( !strcmp ( name, "binary" ) ) {
							binaryFlag = true;
						}
					}
				}
			}
		}
	}
}
void PPExpatMZMLData::endElement ( const char* name )
{
	if ( binaryFlag ) {
		if ( !strcmp ( name, "binary" ) ) binaryFlag = false;
	}
	else if ( binaryDataArray ) {
		if ( !strcmp ( name, "binaryDataArray" ) ) {
			if ( mzArrayFlag ) {
				if ( compression ) {
					if ( precision == 32 ) makeVector32Uncompress ( s, 0, false, mzList );
					else makeVector64Uncompress ( s, 0, false, mzList );
				}
				else {
					if ( precision == 32 ) makeVector32 ( s, false, mzList );
					else makeVector64 ( s, false, mzList );
				}
			}
			else if ( intenArrayFlag ) {
				if ( compression ) {
					if ( precision == 32 ) makeVector32Uncompress ( s, 0, false, intensityList );
					else makeVector64Uncompress ( s, 0, false, intensityList );
				}
				else {
					if ( precision == 32 ) makeVector32 ( s, false, intensityList );
					else makeVector64 ( s, false, intensityList );
				}
			}
			s = "";
			binaryDataArray = false;
			mzArrayFlag = false;
			intenArrayFlag = false;
		}
	}
	else if ( selectedIonFlag ) {
		if ( !strcmp ( name, "selectedIon" ) ) selectedIonFlag = false;
	}
	else if ( scFlag ) {
		if ( !strcmp ( name, "scan" ) ) scFlag = false;
	}
	else if ( scanFlag ) {
		if ( !strcmp ( name, "spectrum" ) ) {
			if ( msLevel == 2 ) {
				if ( !spectrumRange.stopReadingFile ( fraction, count ) ) {
					if ( spectrumRange.startReadingFile ( fraction, count ) ) {
						makeSpectrum ();
					}
				}
			}
			scanFlag = false;
			msLevel = 0;
			precursorMz = 0.0;
			precursorIntensity = 100.0;
			precursorCharge = 0;
		}
	}
	else if ( paramGroupFlag ) {
		if ( !strcmp ( name, "referenceableParamGroup" ) ) {
			paramGroups [paramGroupID] = curParam;
			curParam.clear ();
			paramGroupFlag = false;
		}
	}
}
void PPExpatMZMLData::makeSpectrum ()
{
	msmsDataPoint.clear ();
	for ( DoubleVectorSizeType i = 0 ; i < mzList.size () ; i++ ) {
		if ( i < intensityList.size () ) {
			msmsDataPoint.addPeak ( mzList [i], 1, intensityList [i] );
		}
		else {
			msmsDataPoint.addPeak ( mzList [i], 1, 0.0 );
		}
	}
	PPExpatSpecData::makeSpectrum ();
	mzList.clear ();
	intensityList.clear ();
}
void PPExpatMZMLData::characterDataHandler ( const char* str, int len )
{
	if ( binaryFlag ) {
		s.append ( str, len );
	}
}
/*
psi - Protein Standards Institute
cv - controlled vocabulary

mzData version 1.04
-------------------
<spectrum id="21">
	<acqDesc>
		<acqSettings>
			<acqSpecification spectrumType="discrete" methodOfCombination="sum" count="1">
				<acquisition acqNumber="21" />
			</acqSpecification>
			<acqInstrument msLevel="2" mzRangeStart="100.000000" mzRangeStop="1390.000000">
				<cvParam cvLabel="psi" accession="" name="type" value="full" />
				<cvParam cvLabel="psi" accession="" name="polarity" value="+" />
				<cvParam cvLabel="psi" accession="" name="time.min" value="0.568500" />
			</acqInstrument>
		</acqSettings>
		<precursorList count="1">
			<precursor msLevel="2" spectrumRef="19">
				<ionSelection>
					<cvParam cvLabel="psi" accession="" name="mz" value="344.19" />
				</ionSelection>
				<activation>
					<cvParam cvLabel="psi" accession="" name="method" value="CID" />
					<cvParam cvLabel="psi" accession="" name="energy" value="28.00" />
				</activation>
			</precursor>
		</precursorList>
	</acqDesc>
	<mzArrayBinary>
		<data precision="32" endian="little" length="81">Ug0BQ1...</data>
	</mzArrayBinary>
	<intenArrayBinary>
		<data precision="32" endian="little" length="81">AIBlRA...</data> 
	</intenArrayBinary>
</spectrum>


mzData version 1.05
-------------------
<spectrum id="141">
	<spectrumDesc>
		<spectrumSettings>
			<acqSpecification spectrumType="discrete" methodOfCombination="sum" count="1">
				<acquisition acqNumber="141"/>
			</acqSpecification>
			<spectrumInstrument msLevel="2" mzRangeStart="205.000000" mzRangeStop="2000.000000">
				<cvParam cvLabel="psi" accession="PSI:1000036" name="ScanMode" value="MassScan"/>
				<cvParam cvLabel="psi" accession="PSI:1000037" name="Polarity" value="Positive"/>
				<cvParam cvLabel="psi" accession="PSI:1000038" name="TimeInMinutes" value="3.804667"/>
				<userParam name="ScanType" value="full"/>
			</spectrumInstrument>
		</spectrumSettings>
		<precursorList count="1">
			<precursor msLevel="1" spectrumRef="139">
				<ionSelection>
					<cvParam cvLabel="psi" accession="PSI:1000040" name="MassToChargeRatio" value="661.65"/>
					<cvParam cvLabel="psi" accession="PSI:1000041" name="ChargeState" value="3"/>
				</ionSelection>
				<activation>
					<cvParam cvLabel="psi" accession="PSI:1000044" name="Method" value="CID"/>
					<cvParam cvLabel="psi" accession="PSI:1000045" name="CollisionEnergy" value="28.00"/>
					<cvParam cvLabel="psi" accession="PSI:1000046" name="EnergyUnits" value="Percent"/>
				</activation>
			</precursor>
		</precursorList>
	</spectrumDesc>
	<mzArrayBinary>
		<data precision="32" endian="little" length="331">4h1mQ.....</data>
	</mzArrayBinary>
	<intenArrayBinary>
		<data precision="32" endian="little" length="331">AODYR.....</data>
	</intenArrayBinary>
</spectrum>
*/
PPExpatMZDataData::PPExpatMZDataData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange ) :
	PPExpatSpecData ( msmsDataPointList, fraction, spectrumRange, specID, chrg, chargeRange ),
	scanFlag ( false ),
	ionSelectionFlag ( false ),
	mzArrayFlag ( false ),
	intenArrayFlag ( false ),
	dataFlag ( false )
{
	bigEndian = false;
	precision = 32;
	msLevel = 0;
	precursorMz = 0.0;
	precursorIntensity = 100.0;
	precursorCharge = 0;
	spectrumFlag = false;
}
PPExpatMZDataData::~PPExpatMZDataData () {}
void PPExpatMZDataData::startElement ( const char* name, const char** attributes )
{
	if ( !spectrumFlag ) {
		if ( !strcmp ( name, "spectrum" ) ) {
			getAttributeValue ( attributes, "id", scan );
			spectrumFlag = true;
		}
	}
	else {
		if ( !scanFlag ) {
			if ( !strcmp ( name, "spectrumInstrument" ) ) {	// v1.0.5 (acqInstrument v1.04)
				getAttributeValue ( attributes, "msLevel", msLevel );
				scanFlag = true;
				if ( msLevel == 2 ) count++;
				return;
			}
		}
		else
			getCVAttributeValue ( name, attributes, "TimeInMinutes", rt );	//v1.0.5 (time.min v1.0.4)

		if ( msLevel == 2 ) {
			if ( !spectrumRange.stopReadingFile ( fraction, count ) ) {
				if ( spectrumRange.startReadingFile ( fraction, count ) ) {
					if ( !ionSelectionFlag ) {
						if ( !strcmp ( name, "ionSelection" ) ) ionSelectionFlag = true;
					}
					else {
						getCVAttributeValue ( name, attributes, "MassToChargeRatio", precursorMz );	// v1.0.5 (mz v1.0.4)
						getCVAttributeValue ( name, attributes, "Intensity", precursorIntensity );
						getCVAttributeValue ( name, attributes, "ChargeState", precursorCharge );

						getCVAttributeValue ( name, attributes, "selected ion m/z", precursorMz );	// alternative
						getCVAttributeValue ( name, attributes, "charge state", precursorCharge );

						getCVAttributeValue ( name, attributes, "Mass To Charge Ratio", precursorMz );	// alternative
						getCVAttributeValue ( name, attributes, "Charge State", precursorCharge );
					}
					if ( !mzArrayFlag ) {
						if ( !strcmp ( name, "mzArrayBinary" ) ) mzArrayFlag = true;
					}
					if ( !intenArrayFlag ) {
						if ( !strcmp ( name, "intenArrayBinary" ) ) intenArrayFlag = true;
					}
					if ( !dataFlag ) {
						if ( !strcmp ( name, "data" ) ) {
							dataFlag = true;
							getAttributeValue ( attributes, "precision", precision );
							string endian;
							getAttributeValue ( attributes, "endian", endian );
							bigEndian = (endian == "big");
						}
					}
				}
			}
		}
	}
}
void PPExpatMZDataData::endElement ( const char* name )
{
	if ( scanFlag ) {
		if ( !strcmp ( name, "spectrumInstrument" ) ) scanFlag = false;
	}
	else {
		if ( msLevel == 2 ) {
			if ( ionSelectionFlag ) {
				if ( !strcmp ( name, "ionSelection" ) ) ionSelectionFlag = false;
			}
			else if ( dataFlag ) {
				if ( !strcmp ( name, "data" ) ) {
					if ( mzArrayFlag ) {
						if ( precision == 32 ) makeVector32 ( s, bigEndian, mzList );
						else makeVector64 ( s, bigEndian, mzList );
					}
					else if ( intenArrayFlag ) {
						if ( precision == 32 ) makeVector32 ( s, bigEndian, intensityList );
						else makeVector64 ( s, bigEndian, intensityList );
					}
					s = "";
					dataFlag = false;
				}
			}
			else if ( mzArrayFlag ) {
				if ( !strcmp ( name, "mzArrayBinary" ) ) mzArrayFlag = false;
			}
			else if ( intenArrayFlag ) {
				if ( !strcmp ( name, "intenArrayBinary" ) ) intenArrayFlag = false;
			}
		}
	}
	if ( spectrumFlag ) {
		if ( !strcmp ( name, "spectrum" ) ) {
			if ( msLevel == 2 ) {
				if ( !spectrumRange.stopReadingFile ( fraction, count ) ) {
					if ( spectrumRange.startReadingFile ( fraction, count ) ) {
						makeSpectrum ();
						dPLMap [scan] = msmsDataPointList.size () - 1;
					}
					msLevel = 0;
					precursorMz = 0.0;
					precursorIntensity = 100.0;
					precursorCharge = 0;
				}
			}
			spectrumFlag = false;
		}
	}
}
void PPExpatMZDataData::makeSpectrum ()
{
	msmsDataPoint.clear ();
	for ( DoubleVectorSizeType i = 0 ; i < mzList.size () ; i++ ) {
		msmsDataPoint.addPeak ( mzList [i], 1, intensityList [i] );
	}
	PPExpatSpecData::makeSpectrum ();
	mzList.clear ();
	intensityList.clear ();
}
void PPExpatMZDataData::characterDataHandler ( const char* str, int len )
{
	if ( dataFlag && msLevel == 2 ) s.append ( str, len );
}

PPExpatMascotSearchResults::PPExpatMascotSearchResults ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange ) :
	PPExpatSpecData ( msmsDataPointList, fraction, spectrumRange, specID, chrg, chargeRange )
{
	queryFlag = false;
	queryMOverZFlag = false;
	queryChargeFlag = false;
	queryIntensityFlag = false;
	numValsFlag = false;
	stringIons1Flag = false;
	numVals = 0;
	precursorMz = 0.0;
	precursorIntensity = 100.0;
	precursorCharge = 0;
}
PPExpatMascotSearchResults::~PPExpatMascotSearchResults () {}

void PPExpatMascotSearchResults::startElement ( const char* name, const char** attributes )
{
	if ( !queryFlag ) {
		if ( !strcmp ( name, "query" ) ) {
			getAttributeValue ( attributes, "number", scan );
			rt = scan;
			queryFlag = true;
			count++;
		}
	}
	else {
		if ( !spectrumRange.stopReadingFile ( fraction, count ) ) {
			if ( spectrumRange.startReadingFile ( fraction, count ) ) {
				if ( !strcmp ( name, "query_moverz" ) )			queryMOverZFlag = true;
				else if ( !strcmp ( name, "query_charge" ) )	queryChargeFlag = true;
				else if ( !strcmp ( name, "query_intensity" ) )	queryIntensityFlag = true;
				else if ( !strcmp ( name, "NumVals" ) )			numValsFlag = true;
				else if ( !strcmp ( name, "StringIons1" ) )		stringIons1Flag = true;
			}
		}
	}
}
void PPExpatMascotSearchResults::endElement ( const char* name )
{
	if ( queryFlag ) {
		if ( !strcmp ( name, "query" ) ) {
			if ( !spectrumRange.stopReadingFile ( fraction, count ) ) {
				if ( spectrumRange.startReadingFile ( fraction, count ) ) {
					makeSpectrum ();
				}
				precursorMz = 0.0;
				precursorIntensity = 100.0;
				precursorCharge = 0;
				numVals = 0;
			}
			queryFlag = false;
		}
		if ( queryMOverZFlag && !strcmp ( name, "query_moverz" ) ) {
			queryMOverZFlag = false;
			if ( !s.empty () ) {
				precursorMz = atof ( s.c_str () );
				s = "";
			}
		}
		else if ( queryChargeFlag && !strcmp ( name, "query_charge" ) )	{
			queryChargeFlag = false;
			if ( !s.empty () ) {
				precursorCharge = atoi ( s.c_str () );
				s = "";
			}
		}
		else if ( queryIntensityFlag && !strcmp ( name, "query_intensity" ) ) {
			queryIntensityFlag = false;
			if ( !s.empty () ) {
				precursorIntensity = atof ( s.c_str () );
				s = "";
			}
		}
		else if ( numValsFlag && !strcmp ( name, "NumVals" ) ) {
			numValsFlag = false;
			if ( !s.empty () ) {
				numVals = atoi ( s.c_str () );
				s = "";
			}
		}
		else if ( stringIons1Flag && !strcmp ( name, "StringIons1" ) ) {
			stringIons1Flag = false;
			msmsDataPoint.clear ();
			int start = 0;
			for ( ; ; ) {
				string::size_type end = s.find ( ':', start );
				if ( end == string::npos ) break;
				double mz = atof ( s.substr ( start, end - start ).c_str () );
				start = end + 1;
				end = s.find ( ',', start );
				double inten = atof ( s.substr ( start, end - start ).c_str () );
				msmsDataPoint.addPeak ( mz, 1, inten );
				if ( end == string::npos ) break;
				start = end + 1;
			}
			s = "";
		}
	}
}
void PPExpatMascotSearchResults::characterDataHandler ( const char* str, int len )
{
	if ( queryMOverZFlag || queryChargeFlag || queryIntensityFlag || numValsFlag || stringIons1Flag ) {
		s.append ( str, len );
	}
}
PPExpatPepXML::PPExpatPepXML ( StringVectorVector& headerLines, StringVectorVector& rows, int& zColumnNumber, int& scanIDColumnNumber, int& peptideColumnNumber, int& variableModsColumnNumber ) :
	headerLines ( headerLines ),
	rows ( rows ),
	zColumnNumber ( zColumnNumber ),
	scanIDColumnNumber ( scanIDColumnNumber ),
	peptideColumnNumber ( peptideColumnNumber ),
	variableModsColumnNumber ( variableModsColumnNumber ),
	scoreHeaderSet ( false )
{
}
PPExpatPepXML::~PPExpatPepXML ()
{
}
void PPExpatPepXML::startElement ( const char* name, const char** attributes )
{
	if ( !strcmp ( name, "spectrum_query" ) ) {
		getAttributeValue ( attributes, "precursor_neutral_mass", precursorM );
		getAttributeValue ( attributes, "assumed_charge", z );
		getAttributeValue ( attributes, "index", index );
		getAttributeValue ( attributes, "start_scan", scanNumber );
		return;
	}
	if ( !strcmp ( name, "search_hit" ) ) {
		getAttributeValue ( attributes, "peptide", dbPeptide );
		mod = "";
		return;
	}
	if ( !strcmp ( name, "modification_info" ) ) {
		string modNTermMass;
		string modCTermMass;
		getAttributeValue ( attributes, "mod_nterm_mass", modNTermMass );
		getAttributeValue ( attributes, "mod_cterm_mass", modCTermMass );
		if ( !modNTermMass.empty () ) {
			static double hAdjust = formula_to_monoisotopic_mass ( "H" );
			double mass = atof ( modNTermMass.c_str () );
			mass -= hAdjust;
			mod += gen_ftoa ( mass, "%.4f" ) + "@N-term;";
		}
		if ( !modCTermMass.empty () ) {
			static double ohAdjust = formula_to_monoisotopic_mass ( "O H" );
			double mass = atof ( modCTermMass.c_str () );
			mass -= ohAdjust;
			mod += gen_ftoa ( mass, "%.4f" ) + "@C-term;";
		}
		return;
	}
	if ( !strcmp ( name, "mod_aminoacid_mass" ) ) {
		string position;
		double mass;
		getAttributeValue ( attributes, "position", position );
		getAttributeValue ( attributes, "mass", mass );
		mass -= AAInfo::getInfo ().getMonoisotopicMass ( dbPeptide[atoi(position.c_str ())-1] );
		mod += gen_ftoa ( mass, "%.4f" ) + "@" + position + ";";
		return;
	}
	if ( !strcmp ( name, "search_score" ) ) {
		double sc;
		getAttributeValue ( attributes, "value", sc );
		if ( !scoreHeaderSet ) {
			string scHeader;
			getAttributeValue ( attributes, "name", scHeader );
			scoreHeaders.push_back ( scHeader );
		}
		scores.push_back ( sc );
	}
}
void PPExpatPepXML::endElement ( const char* name )
{
	if ( !strcmp ( name, "search_hit" ) ) {
		StringVector cols;

		cols.push_back ( gen_ftoa ( mToMOverZ ( precursorM, z, true ), "%.4f" ) );

		cols.push_back ( gen_itoa ( z ) );
		zColumnNumber = cols.size ();

		cols.push_back ( gen_itoa ( index ) );

		cols.push_back ( gen_itoa ( scanNumber ) );
		scanIDColumnNumber = cols.size ();

		cols.push_back ( dbPeptide );
		peptideColumnNumber = cols.size ();

		if ( !mod.empty () ) mod = mod.substr ( 0, mod.length () - 1 ); // Delete trailing ";"
		cols.push_back ( mod );
		variableModsColumnNumber = cols.size ();

		for ( int i = 0 ; i < scores.size () ; i++ ) {
			ostringstream ostr;
			genPrintSigFig ( ostr, scores [i], 3 );
			cols.push_back ( ostr.str () );
		}
		scores.clear ();

		if ( rows.empty () ) {
			StringVector header;
			header.push_back ( "m/z" );
			header.push_back ( "z" );
			header.push_back ( "Scan ID" );
			header.push_back ( "Scan Number" );
			header.push_back ( "DB Peptide" );
			header.push_back ( "Mods" );
			for ( int j = 0 ; j < scoreHeaders.size () ; j++ ) {
				header.push_back ( scoreHeaders [j] );
			}
			scoreHeaderSet = true;
			headerLines.push_back ( header );
		}
		rows.push_back ( cols );
	}
}


PPExpatPrideXML::PPExpatPrideXML ( StringVector& constantHeader, SetString& variableHeader, StringVectorVector& rows, VectorMapStringToString& rowsVariable, MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange ) :
	PPExpatMZDataData ( msmsDataPointList, fraction, spectrumRange, specID, chrg, chargeRange ),
	constantHeader ( constantHeader ),
	variableHeader ( variableHeader ),
	rows ( rows ),
	rowsVariable ( rowsVariable )
{
	mzDataFlag = false;
	gelFreeIDFlag = false;
	peptideItemFlag = false;
	accessionFlag = false;
	startFlag = false;
	endFlag = false;
	spectrumReferenceFlag = false;
	sequenceFlag = false;
	modificationFlag = false;
	additionalFlag = false;
	modLocationFlag = false;
	modMonoDeltaFlag = false;
}
PPExpatPrideXML::~PPExpatPrideXML ()
{
}
void PPExpatPrideXML::startElement ( const char* name, const char** attributes )
{
	if ( mzDataFlag ) {
		PPExpatMZDataData::startElement ( name, attributes );
	}
	else {
		if ( !strcmp ( name, "mzData" ) ) {
			mzDataFlag = true;
		}
		else {
			if ( gelFreeIDFlag ) {
				if ( peptideItemFlag ) {
					if ( modificationFlag ) {
						if ( !strcmp ( name, "ModLocation" ) ) {
							modLocationFlag = true;
						}
						else if ( !strcmp ( name, "ModMonoDelta" ) ) {
							modMonoDeltaFlag = true;
						}
					}
					else if ( additionalFlag ) {
						string n;
						string v;
						if ( getCVParamNameAndValue ( name, attributes, n, v ) ) {
							att [n] = v;
						}
						else if ( getUserParamNameAndValue ( name, attributes, n, v ) ) {
							att [n] = v;
						}
					}
					else {
						if ( !strcmp ( name, "Start" ) ) {
							startFlag = true;
						}
						else if ( !strcmp ( name, "End" ) ) {
							endFlag = true;
						}
						else if ( !strcmp ( name, "SpectrumReference" ) ) {
							spectrumReferenceFlag = true;
						}
						else if ( !strcmp ( name, "Sequence" ) ) {
							sequenceFlag = true;
						}
						else if ( !strcmp ( name, "ModificationItem" ) ) {
							modificationFlag = true;
						}
						else if ( !strcmp ( name, "additional" ) ) {
							additionalFlag = true;
						}
					}
				}
				else {
					if ( !strcmp ( name, "PeptideItem" ) ) {
						peptideItemFlag = true;
					}
					else if ( !strcmp ( name, "Accession" ) ) {
						accessionFlag = true;
					}
				}
			}
			else {
				if ( !strcmp ( name, "GelFreeIdentification" ) ) {
					gelFreeIDFlag = true;
				}
			}
		}
	}
}
void PPExpatPrideXML::endElement ( const char* name )
{
	if ( additionalFlag ) {
		if ( !strcmp ( name, "additional" ) ) {
			additionalFlag = false;
		}
	}
	else if ( modificationFlag ) {
		if ( !strcmp ( name, "ModificationItem" ) ) {
			modificationFlag = false;
		}
		else if ( modLocationFlag ) {
			if ( !strcmp ( name, "ModLocation" ) ) {
				modLocation.push_back ( s );
				s = "";
				modLocationFlag = false;
			}
		}
		else if ( modMonoDeltaFlag ) {
			if ( !strcmp ( name, "ModMonoDelta" ) ) {
				modMonoDelta.push_back ( s );
				s = "";
				modMonoDeltaFlag = false;
			}
		}
	}
	else if ( sequenceFlag ) {
		if ( !strcmp ( name, "Sequence" ) ) {
			sequence = s;
			s = "";
			sequenceFlag = false;
		}
	}
	else if ( spectrumReferenceFlag ) {
		if ( !strcmp ( name, "SpectrumReference" ) ) {
			spectrumReference = s;
			s = "";
			spectrumReferenceFlag = false;
		}
	}
	else if ( endFlag ) {
		if ( !strcmp ( name, "End" ) ) {
			end = s;
			s = "";
			endFlag = false;
		}
	}
	else if ( startFlag ) {
		if ( !strcmp ( name, "Start" ) ) {
			start = s;
			s = "";
			startFlag = false;
		}
	}
	else if ( peptideItemFlag ) {
		if ( !strcmp ( name, "PeptideItem" ) ) {
			StringVector cols;
			cols.push_back ( accession );
			cols.push_back ( start );
			cols.push_back ( end );
			cols.push_back ( spectrumReference );
			cols.push_back ( sequence );
			string mod;
			for ( int i = 0 ; i < modLocation.size () ; i++ ) {
				double mass = atof ( modMonoDelta [i].c_str () );
				string loc = modLocation [i];
				if ( loc == "0" ) loc = "N-term";
				mod += gen_ftoa ( mass, "%.4f" ) + "@" + loc  + ";";
			}
			if ( !mod.empty () ) mod = mod.substr ( 0, mod.length () - 1 ); // Delete trailing ";"
			cols.push_back ( mod );
			MapStringToIntConstIterator cur = dPLMap.find ( spectrumReference );
			int idx = (*cur).second;
			cols.push_back ( gen_ftoa ( msmsDataPointList [idx].getPrecursorMZ (), "%.4f" ) );
			cols.push_back ( gen_itoa ( msmsDataPointList [idx].getPrecursorCharge () ) );
			peptideItemFlag = false;
			if ( rows.empty () ) {
				constantHeader.push_back ( "Accession" );
				constantHeader.push_back ( "Start" );
				constantHeader.push_back ( "End" );
				constantHeader.push_back ( "Scan ID" );
				constantHeader.push_back ( "DB Peptide" );
				constantHeader.push_back ( "Mods" );
				constantHeader.push_back ( "m/z" );
				constantHeader.push_back ( "z" );
			}
			for ( MapStringToStringConstIterator k = att.begin () ; k != att.end () ; k++ ) {
				variableHeader.insert ( (*k).first );
			}
			start = "";
			end = "";
			spectrumReference = "";
			sequence = "";
			modLocation.clear ();
			modMonoDelta.clear ();
			rows.push_back ( cols );
			rowsVariable.push_back ( att );
			att.clear ();
		}
	}
	else if ( accessionFlag ) {
		if ( !strcmp ( name, "Accession" ) ) {
			accession = s;
			s = "";
			accessionFlag = false;
		}
	}
	else if ( gelFreeIDFlag ) {
		if ( !strcmp ( name, "GelFreeIdentification" ) ) {
			accession = "";
			gelFreeIDFlag = false;
		}
	}
	else if ( mzDataFlag ) {
		if ( !strcmp ( name, "mzData" ) ) {
			mzDataFlag = false;
		}
		else {
			PPExpatMZDataData::endElement ( name );
		}
	}
}
void PPExpatPrideXML::characterDataHandler ( const char* str, int len )
{
	if ( mzDataFlag ) {
		PPExpatMZDataData::characterDataHandler ( str, len );
	}
	else {
		if ( startFlag || endFlag || spectrumReferenceFlag || sequenceFlag || modLocationFlag || modMonoDeltaFlag || accessionFlag ) {
			s.append ( str, len );
		}
	}
}
