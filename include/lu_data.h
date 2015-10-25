/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_data.h                                                     *
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

#ifndef __lu_data_h
#define __lu_data_h

#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <lgen_define.h>

class SpecID;
class Tolerance;
class ParameterList;
class GenIFStream;
class XYData;

class DataFilePeak {
	double mOverZ;
	double intensity;
	int charge;
	friend class sortDataFilePeaksByMOverZ;
public:
	DataFilePeak ( double mOverZ, int charge, double intensity ) :
		mOverZ ( mOverZ ), charge ( charge ), intensity ( intensity ) {}
	double getMOverZ () const { return mOverZ; }
	double getIntensity () const { return intensity; }
	void setIntensity ( double i ) { intensity = i; }
	void setZeroIntensity () { intensity = 0.0; }
	int getCharge () const { return charge; }
};

class sortDataFilePeaksByMOverZ {
public:
	bool operator () ( const DataFilePeak& a, const DataFilePeak& b ) const
	{
		return ( a.mOverZ < b.mOverZ );
	}
};

typedef std::vector <DataFilePeak> DataFilePeakVector;
typedef DataFilePeakVector::size_type DataFilePeakVectorSizeType;

class DataPoint {
protected:
	int fraction;			// lc fraction or spot plate
	std::string spot;		// spot number or elution time
	int run;				// different runs on same spot
	DataFilePeakVector dataPeaks;
	bool sorted;
public:
	DataPoint ()
	{
		sorted = true;
	}
	virtual ~DataPoint ();
	virtual void setPointInfo ( int f, const std::string& s, int r )
	{
		fraction = f;
		spot = s;
		run = r;
	}
	void addPeak ( double m, int c, double i )
	{
		if ( !dataPeaks.empty () ) {
			if ( m < dataPeaks.back ().getMOverZ () ) sorted = false;
		}
		dataPeaks.push_back ( DataFilePeak ( m, c, i ) );
	}
	void clear ()
	{
		sorted = true;
		dataPeaks.clear ();
	}
	DataFilePeakVector& getDataPeaks ()
	{
		if ( !sorted ) std::sort ( dataPeaks.begin (), dataPeaks.end (), sortDataFilePeaksByMOverZ () );
		return dataPeaks;
	}
	int size () const { return dataPeaks.size (); }
	std::string getSpotID () const;
	std::string getSpotRunID () const;
	virtual std::string getSpecID () const;
	int getFraction () const { return fraction; }
	std::string getSpot () const { return spot; };
};

class MSDataPoint : public DataPoint {
public:
	MSDataPoint () {}
};

class MSMSDataPoint : public DataPoint {
	int spectrumNumber;		// different MS/MS spectra on same spot
	std::string msmsInfo;	// info to construct MSMS raw data
	double precursorMZ;
	int precursorCharge;
	IntVector precursorChargeRange;
	double precursorIntensity;
	double precursorTolerance;
public:
	MSMSDataPoint () {}
	~MSMSDataPoint () {}
	virtual void setPointInfo ( int f, const std::string& s, int r, int n, const std::string& m )
	{
		DataPoint::setPointInfo ( f, s, r );
		spectrumNumber = n;
		msmsInfo = m;
	}
	void setPrecursor ( double m, int c, double i )
	{
		precursorMZ = m;
		precursorCharge = c;
		precursorIntensity = i;
	}
	void setPrecursorChargeRange ( const IntVector& cr )
	{
		precursorChargeRange = cr;
	}
	void setPrecursorCharge ( int c ) { precursorCharge = c; }
	void setPrecursorMZ ( double m ) { precursorMZ = m; }
	void setPrecursorIntensity ( double i ) { precursorIntensity = i; }
	void setSpectrumNumber ( int n ) { spectrumNumber = n; }
	double getParentMPlusH () const;
	double getPrecursorMZ () const { return precursorMZ; }
	int getPrecursorCharge () const { return precursorCharge; }
	IntVector getPrecursorChargeRange () const { return precursorChargeRange; }
	int getSpectrumNumber () const { return spectrumNumber; }
	double getPrecursorIntensity () const { return precursorIntensity; }
	double getPrecursorTolerance () const { return precursorTolerance; }
	void calibrate ( const Tolerance* tol, double offset );
	void setTolerance ( const Tolerance* tol );
	void printSpecInfoDelimited ( std::ostream& os ) const;
	std::string getSpecID () const;
	void putCGI ( std::ostream& os ) const;
	std::string getMSMSInfo () const { return msmsInfo; };
};

typedef std::vector <MSDataPoint> MSDataPointVector;
typedef MSDataPointVector::size_type MSDataPointVectorSizeType;

typedef std::vector <MSMSDataPoint> MSMSDataPointVector;
typedef MSMSDataPointVector::size_type MSMSDataPointVectorSizeType;

class SpectrumRange {
	int startFraction;
	int startSpectrum;
	int endFraction;
	int endSpectrum;
	bool allFractions;
public:
	SpectrumRange () : allFractions ( true ) {}
	SpectrumRange ( const ParameterList* params );
	bool readFile ( int fraction ) const { return allFractions || ( fraction >= startFraction && fraction <= endFraction ); }
	bool stopReadingFile ( int fraction, int spectrumNo ) const { return !allFractions && fraction == endFraction && spectrumNo > endSpectrum; }
	bool startReadingFile ( int fraction, int spectrumNo ) const { return allFractions || fraction != startFraction || spectrumNo >= startSpectrum; }
	bool stopReadingFile ( int spectrumNo ) const { return !allFractions || spectrumNo > endSpectrum; }
	bool startReadingFile ( int spectrumNo ) const { return allFractions || spectrumNo >= startSpectrum; }
	int getStartFraction () const { return startFraction; }
	int getStartSpectrum () const { return startSpectrum; }
	int getEndFraction () const { return endFraction; }
	int getEndSpectrum () const { return endSpectrum; }
	bool getAllFractions () const { return allFractions; }
};

class UpdatingJavascriptMessage;

class DataReader {
	bool intensityFlag;
	GenIFStream* ifstr;			// the input file stream
	std::istringstream isstr;	// the input string stream
	GenIFStream* initFStream ( const std::string& s, bool fileFlag );
	std::istream& initStream ( bool fileFlag );
	virtual void readParentData ( int charge, int isotopeOffset );
	virtual void skipParentData ();
	virtual void readData ();
	static const int MAX_CHARGE;
	static std::string lastFilename;
protected:
	int num;
	std::streampos spos;
	std::streampos endPos;
	IntVector chargeRange;
	MapStringToInt spotNumber;
	MapStringToStreampos mapSpecID;
	void readLine ( char* ptr );
	DataPoint* dataPoint;
	MSMSDataPoint msmsDataPoint;
	MSDataPoint msDataPoint;
	static std::string lastVersion;
	std::istream& istr;			// the input stream

	void readMSData ( MSDataPointVector& msDataPointList );
	void readMSMSData ( MSMSDataPointVector& msmsDataPointList, int charge, int isotopeOffset );
	bool readMSMSData ( MSMSDataPointVector& msmsDataPointList, int charge, int isotopeOffset, double mOverZ );
	void readMSMSQuanData ( XYData& xyData, double startMass, double endMass );
public:
	DataReader ( const std::string& s, bool fileFlag, bool intensityFlag );
	~DataReader ();
	virtual void getData ( MSMSDataPointVector& msmsDataPointList, const SpecID& specID, int chrg, const std::string& version, UpdatingJavascriptMessage* ujm );
	virtual void getData ( MSDataPointVector& msDataPointList, int fraction, const SpecID& specID ) = 0;
	virtual void getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const std::string& version, int off = 0 ) = 0;
	virtual void getDataFromTitle ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const std::string& title, int chrg, int off = 0 ) {} // STUB
	virtual void getDataFromIndex ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, int index, int chrg, int off = 0 ) {} // STUB
	virtual void getDataFromMOverZ ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, double mOverZ, int chrg, int off = 0 ) {} // STUB
	virtual void getDataFromScanNumber ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const std::string& scanNumber, int chrg, int off = 0 ) {} // STUB
	void setChargeRange ( const IntVector& iv ) { chargeRange = iv; }
	virtual void writePeakList ( std::ostream& os ) {}
	virtual bool writePeakList ( std::ostream& os, const SpecID& specID, const std::string& version );
	virtual bool getMSMSQuanData ( XYData& xyData, const SpecID& specID, double startMass, double endMass, const std::string& version );
	static std::string getLastFilename () { return lastFilename; }
	static std::string getLastVersion () { return lastVersion; }
	void rewind ();
	virtual bool isEOF ();
};

#ifdef BATCHTAG
class XMLDataReader : public DataReader {
	const std::string& s;
	bool fileFlag;
	std::string dataFormat;
public:
	XMLDataReader ( const std::string& s, bool fileFlag, const std::string& dataFormat );
	void getData ( MSDataPointVector& msDataPointList, int fraction, const SpecID& specID ) {}
	void getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const std::string& version, int off = 0 );
};
#endif

class MGFDataReader : public DataReader {
	void readParentData ( int charge, int isotopeOffset );
	void skipParentData ();
	bool readSpectrumIfAppropriate ( MSMSDataPointVector& msmsDataPointList, const SpecID& specID, int fraction, const std::string& spot, int run, int spectrumNumber, const std::string& msmsInfo, std::streampos& spos2, int chrg );
	bool getMSMSTitleParams ( const std::string& line, std::string& spot, int& run, std::string& msmsInfo );
	IntVector getScansFromScanLine ();
	void initMaps ( int fraction, int chrg, const std::string& version, UpdatingJavascriptMessage* ujm );
	void getRTInSeconds ( std::string& rt );
public:
	MGFDataReader ( const std::string& s, bool fileFlag );
	void getData ( MSMSDataPointVector& msmsDataPointList, const SpecID& specID, int chrg, const std::string& version, UpdatingJavascriptMessage* ujm );
	void getData ( MSDataPointVector& msDataPointList, int fraction, const SpecID& specID ) {}
	void getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const std::string& version, int off = 0 );
	void getDataFromTitle ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const std::string& title, int chrg, int off = 0 );
	void getDataFromIndex ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, int index, int chrg, int off = 0 );
	void getDataFromMOverZ ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, double mOverZ, int chrg, int off = 0 );
	void getDataFromScanNumber ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const std::string& scanNumber, int chrg, int off = 0 );
	void writePeakList ( std::ostream& os );
	bool writePeakList ( std::ostream& os, const SpecID& specID, const std::string& version );
	bool getMSMSQuanData ( XYData& xyData, const SpecID& specID, double startMass, double endMass, const std::string& version );
};

class MS2DataReader : public DataReader {
	bool readRT;
	bool readScan;
	bool parseHeader () const;
	void parseHeaderLines ( bool bullseye, std::string& spot, int& spectrumNumber, std::string& msmsInfo, DoubleVector& mVector, IntVector& zVector );
	bool readMSMSData ();
	bool skipMSMSData ();
public:
	MS2DataReader ( const std::string& s, bool fileFlag );
	void getData ( MSDataPointVector& msDataPointList, int fraction, const SpecID& specID ) {}
	void getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const std::string& version, int off = 0 );
	void getDataFromIndex ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, int index, int chrg, int off = 0 );
	void getDataFromScanNumber ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const std::string& scanNumber, int chrg, int off = 0 );
	void writePeakList ( std::ostream& os );
	bool writePeakList ( std::ostream& os, const SpecID& specID, const std::string& version );
};

class APLDataReader : public DataReader {
	void parseHeaderLines ( std::string& spot, int& spectrumNumber, std::string& msmsInfo, DoubleVector& mVector, IntVector& zVector );
	bool readMSMSData ();
	bool skipMSMSData ();
public:
	APLDataReader ( const std::string& s, bool fileFlag );
	void getData ( MSDataPointVector& msDataPointList, int fraction, const SpecID& specID ) {}
	void getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const std::string& version, int off = 0 );
	void getDataFromIndex ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, int index, int chrg, int off = 0 );
	void getDataFromScanNumber ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const std::string& scanNumber, int chrg, int off = 0 );
	void writePeakList ( std::ostream& os );
	bool writePeakList ( std::ostream& os, const SpecID& specID, const std::string& version );
};

class PPDataReader : public DataReader {
	static std::string getParameter ( const std::string& line, const std::string& param );
	static int getIntValue ( const std::string& str, const std::string& tag );
	static double getDoubleValue ( const std::string& str, const std::string& tag );
	static std::string getStringValue ( const std::string& str, const std::string& tag );
	MapStringToInt spotNumber;
	inline bool getMSTitleParams ( const std::string& line, std::string& spot, int& run );
	inline bool getMSMSTitleParams ( const std::string& line, std::string& sNum, int& run, int& spectrumNumber, std::string& msmsInfo );
public:
	PPDataReader ( const std::string& s, bool fileFlag );
	void getData ( MSDataPointVector& msDataPointList, int fraction, const SpecID& specID );
	void getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const std::string& version, int off = 0 );
};
inline bool PPDataReader::getMSTitleParams ( const std::string& line, std::string& spot, int& run )
{
	if ( line.substr ( 1, 2 ) == "M1" ) {
		spot = getStringValue ( line, "Spot:" );
		run = getIntValue ( line, "Run:" );
		return true;
	}
	else return false;
}
inline bool PPDataReader::getMSMSTitleParams ( const std::string& line, std::string& spot, int& run, int& spectrumNumber, std::string& msmsInfo )
{
	if ( line.substr ( 1, 2 ) == "M2" ) {
		spot = getStringValue ( line, "Spot:" );
		run = getIntValue ( line, "Run:" );
		msmsInfo = getStringValue ( line, "JobRunItem:" );
		std::ostringstream spotRun;
		spotRun << spot << "-" << run;
		MapStringToIntIterator cur = spotNumber.find ( spotRun.str () );
		if ( cur != spotNumber.end () )
			cur->second++;
		else
			spotNumber [spotRun.str ()] = 1;
		cur = spotNumber.find ( spotRun.str () );
		spectrumNumber = cur->second;
		return true;
	}
	else return false;
}

class PPDataReader2 : public DataReader {
public:
	PPDataReader2 ( const std::string& s, bool fileFlag, const std::string& dataFormat );
	void getData ( MSDataPointVector& msDataPointList, int fraction, const SpecID& specID );
	void getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const std::string& version, int off = 0 );
};

class SpaceSeparatedSpectraDataReader : public DataReader {
	void readData ();
public:
	SpaceSeparatedSpectraDataReader ( const std::string& s, bool fileFlag );
	void getData ( MSDataPointVector& msDataPointList, int fraction, const SpecID& specID ) {}
	void getData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const std::string& version, int off = 0 );
	bool getMSMSQuanData ( XYData& xyData, const SpecID& specID, double startMass, double endMass, const std::string& version );
};

class PKLDataReader : public SpaceSeparatedSpectraDataReader {
	void readParentData ( int charge, int isotopeOffset );
	void skipParentData ();
public:
	PKLDataReader ( const std::string& s, bool fileFlag );
};

class DTADataReader : public SpaceSeparatedSpectraDataReader {
	void readParentData ( int charge, int isotopeOffset );
	void skipParentData ();
public:
	DTADataReader ( const std::string& s, bool fileFlag );
};

#endif /* ! __lu_data_h */
