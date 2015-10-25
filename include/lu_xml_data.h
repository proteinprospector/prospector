/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_xml_data.h                                                 *
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

#ifndef __lu_xml_data_h
#define __lu_xml_data_h

#include <lgen_xml.h>
#include <lu_data.h>

class PPExpatSpecData : public PPExpat {
protected:
	MSMSDataPointVector& msmsDataPointList;
	const SpectrumRange& spectrumRange;
	const SpecID& specID;
	int chrg;
	IntVector chargeRange;
	MapStringToInt spotNumber;
	MSMSDataPoint msmsDataPoint;
	double precursorMz;
	double precursorIntensity;
	int precursorCharge;
	int fraction;
	bool bigEndian;	// Network or big endian byte order (ie Most Significant Byte first)
	int precision;
	int msLevel;
	std::string rt;
	std::string scan;
	int count;
protected:
	virtual void makeSpectrum ();
public:
	PPExpatSpecData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange );
	virtual ~PPExpatSpecData ();
};

class PPExpatCountScans : public PPExpat {
protected:
	int count;
public:
	PPExpatCountScans () : count ( 0 ) {}
	~PPExpatCountScans () {}
	int getCount () const { return count; }
};

class PPExpatMZXMLData : public PPExpatSpecData {
	bool peaksFlag;
	bool precursorMzFlag;
	int peaksCount;
	bool compression;

	void startElement ( const char* name, const char** attributes );
	void endElement ( const char* name );
	void characterDataHandler ( const char* str, int len );
public:
	PPExpatMZXMLData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange );
	~PPExpatMZXMLData ();
};

class PPExpatMZXMLCountScans : public PPExpatCountScans {
	void startElement ( const char* name, const char** attributes )
	{
		if ( !strcmp ( name, "scan" ) ) {
			int msLevel = 0;
			getAttributeValue ( attributes, "msLevel", msLevel );
			if ( msLevel == 2 ) count++;
		}
	}
public:
	PPExpatMZXMLCountScans () {}
	~PPExpatMZXMLCountScans () {}
};

class PPExpatMZMLData : public PPExpatSpecData {
	bool scanFlag;
	bool paramGroupFlag;
	bool selectedIonFlag;
	bool compression;
	bool mzArrayFlag;
	bool intenArrayFlag;
	bool chargeArrayFlag;
	bool binaryDataArray;
	bool binaryFlag;
	bool scFlag;

	std::string paramGroupID;
	std::map <std::string, MapStringToString> paramGroups;
	MapStringToString curParam;

	DoubleVector mzList;
	DoubleVector intensityList;
	DoubleVector chargeList;
	void startElement ( const char* name, const char** attributes );
	void endElement ( const char* name );
	void characterDataHandler ( const char* str, int len );
	void makeSpectrum ();
public:
	PPExpatMZMLData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange );
	~PPExpatMZMLData ();
};

class PPExpatMZMLCountScans : public PPExpatCountScans {
	void startElement ( const char* name, const char** attributes )
	{
		int msLevel = 0;
		if ( getCVAttributeValue ( name, attributes, "ms level", msLevel ) ) {
			if ( msLevel == 2 ) count++;
		}
	}
public:
	PPExpatMZMLCountScans () {}
	~PPExpatMZMLCountScans () {}
};

class PPExpatMZDataData : public PPExpatSpecData {
	bool spectrumFlag;
	bool scanFlag;
	bool ionSelectionFlag;
	bool mzArrayFlag;
	bool intenArrayFlag;
	bool dataFlag;

	DoubleVector mzList;
	DoubleVector intensityList;
	void makeSpectrum ();
protected:
	MapStringToInt dPLMap;
	void startElement ( const char* name, const char** attributes );
	void endElement ( const char* name );
	void characterDataHandler ( const char* str, int len );
public:
	PPExpatMZDataData ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange );
	~PPExpatMZDataData ();
};

class PPExpatMZDataCountScans : public PPExpatCountScans {
	void startElement ( const char* name, const char** attributes )
	{
		if ( !strcmp ( name, "spectrumInstrument" ) || !strcmp ( name, "acqInstrument" ) ) {
			int msLevel = 0;
			getAttributeValue ( attributes, "msLevel", msLevel );
			if ( msLevel == 2 ) count++;
		}
	}
public:
	PPExpatMZDataCountScans () {}
	~PPExpatMZDataCountScans () {}
};

class PPExpatMascotSearchResults : public PPExpatSpecData {
	bool queryFlag;
	bool queryMOverZFlag;
	bool queryChargeFlag;
	bool queryIntensityFlag;
	bool numValsFlag;
	bool stringIons1Flag;
	int numVals;
	void startElement ( const char* name, const char** attributes );
	void endElement ( const char* name );
	void characterDataHandler ( const char* str, int len );
public:
	PPExpatMascotSearchResults ( MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange );
	~PPExpatMascotSearchResults ();
};

class PPExpatMascotSearchResultsCountScans : public PPExpatCountScans {
	void startElement ( const char* name, const char** attributes )
	{
		if ( !strcmp ( name, "query" ) ) count++;
	}
public:
	PPExpatMascotSearchResultsCountScans () {}
	~PPExpatMascotSearchResultsCountScans () {}
};

class PPExpatPepXML : public PPExpat {

	std::string fractionName;

	StringVectorVector& headerLines;
	StringVectorVector& rows;

	int& zColumnNumber;
	int& fractionColumnNumber;
	int& scanIDColumnNumber;
	int& peptideColumnNumber;
	int& variableModsColumnNumber;

	int indexFlag;

	std::string dbPeptide;
	std::string mod;
	double precursorM;
	int z;
	int index;
	int scanNumber;
	bool precursorMOverZFlag;
	double mOverZ;
	StringVector scores;
	StringVector scoreHeaders;
	MapStringToInt scoreHeadersMap;
	int scoreHeadersSize;

	void startElement ( const char* name, const char** attributes );
	void endElement ( const char* name );
public:
	PPExpatPepXML ( const std::string& fractionName, StringVectorVector& headerLines, StringVectorVector& rows, int& zColumnNumber, int& fractionColumnNumber, int& scanIDColumnNumber, int& peptideColumnNumber, int& variableModsColumnNumber );
	~PPExpatPepXML ();
};

class PPExpatPrideXML : public PPExpatMZDataData {

	StringVector& constantHeader;
	SetString& variableHeader;
	StringVectorVector& rows;
	VectorMapStringToString& rowsVariable;

	bool mzDataFlag;
	bool gelFreeIDFlag;
	bool peptideItemFlag;
	bool accessionFlag;
	bool startFlag;
	bool endFlag;
	bool spectrumReferenceFlag;
	bool sequenceFlag;
	bool modificationFlag;
	bool additionalFlag;
	bool modLocationFlag;
	bool modMonoDeltaFlag;

	std::string accession;
	std::string start;
	std::string end;
	std::string spectrumReference;
	std::string sequence;
	StringVector modLocation;
	StringVector modMonoDelta;
	MapStringToString att;

	void startElement ( const char* name, const char** attributes );
	void endElement ( const char* name );
	void characterDataHandler ( const char* str, int len );
public:
	PPExpatPrideXML ( StringVector& constantHeader, SetString& variableHeader, StringVectorVector& rows, VectorMapStringToString& rowsVariable, MSMSDataPointVector& msmsDataPointList, int fraction, const SpectrumRange& spectrumRange, const SpecID& specID, int chrg, const IntVector& chargeRange );
	~PPExpatPrideXML ();
};

#endif /* ! __lu_xml_data_h */
