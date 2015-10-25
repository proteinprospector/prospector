/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_df_info.h                                                  *
*                                                                             *
*  Created    : August 8th 2005                                               *
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

#ifndef __lu_df_info_h
#define __lu_df_info_h

#include <string>
#include <lu_data.h>

class ParameterList;

class DataSetInfo {
protected:
	std::string dataFormat;
	std::string dataSource;
	std::string searchKey;
	std::string dataFilename;
	DataReader* dr;
	void parseData ( const ParameterList* params, bool noFilesOK, const SpecID& specID, bool msmsData, int off = 0 );
	void parseData ( const ParameterList* params, const std::string& title, int off = 0 );
	void parseData ( const ParameterList* params, int index, int off = 0 );
	void parseData ( const ParameterList* params, double mOverZ, int off = 0 );
	void parseData ( const ParameterList* params, const std::string& scanNumber, bool dummy, int off = 0 );
	void initDataReader ( const std::string& s, bool fileFlag );
	virtual void getData ( const std::string& s, bool fileFlag, int fraction, const SpecID& specID, const std::string& version, int off = 0 ) = 0;
	virtual void getDataFromTitle ( const std::string& s, bool fileFlag, int fraction, const std::string& title, int off = 0 ) {}
	virtual void getDataFromIndex ( const std::string& s, bool fileFlag, int fraction, int index, int off = 0 ) {}
	virtual void getDataFromMOverZ ( const std::string& s, bool fileFlag, int fraction, double mOverZ, int off = 0 ) {}
	virtual void getDataFromScanNumber ( const std::string& s, bool fileFlag, int fraction, const std::string& scanNumber, int off = 0 ) {}
public:
	DataSetInfo ( const std::string& dataFormat, const std::string& dataSource, const std::string& searchKey );
	DataSetInfo ( const ParameterList* params );
	DataSetInfo ();
	virtual ~DataSetInfo () = 0;

	virtual int getNumDataSets () const = 0;
};

class MSDataSetInfo : public DataSetInfo {
	MSDataPointVector dataPointList;
	void getData ( const std::string& s, bool fileFlag, int fraction, const SpecID& specID, const std::string& version, int off = 0 );
public:
	MSDataSetInfo ( const ParameterList* params );
	~MSDataSetInfo ();
	bool getNoDataFlag () const { return dataPointList.empty (); }
	int getNumDataSets () const { return dataPointList.size (); }
	MSDataPoint* getDataSet ( int i ) { return &dataPointList [i]; }
	int getDataSetSize ( int i ) { return dataPointList [i].size (); }
	std::string getSpotID ( int i ) const { return dataPointList [i].getSpotID (); }
	std::string getSpotRunID ( int i ) const { return dataPointList [i].getSpotRunID (); }
	std::string getSpecID ( int i ) const { return dataPointList [i].getSpecID (); }
};

class MSMSDataSetInfo : public DataSetInfo {
	SpectrumRange spectrumRange;
	int charge;
	IntVector chargeRange;
	MSMSDataPointVector dataPointList;
	std::string title;
	std::string index;
	std::string mOverZ;
	std::string scanNumber;
	int setCharge ( const ParameterList* params ) const;
	void calibrate ( const ParameterList* params );
	void getData ( const std::string& s, bool fileFlag, int fraction, const SpecID& specID, const std::string& version, int off = 0 );
	void getDataFromTitle ( const std::string& s, bool fileFlag, int fraction, const std::string& title, int off = 0 );
	void getDataFromIndex ( const std::string& s, bool fileFlag, int fraction, int index, int off = 0 );
	void getDataFromMOverZ ( const std::string& s, bool fileFlag, int fraction, double mOverZ, int off = 0 );
	void getDataFromScanNumber ( const std::string& s, bool fileFlag, int fraction, const std::string& scanNumber, int off = 0 );
public:
	MSMSDataSetInfo ( const ParameterList* params, bool noFilesOK );
	~MSMSDataSetInfo ();
	bool getNoDataFlag () const { return dataPointList.empty (); }
	int getNumDataSets () const { return dataPointList.size (); }
	MSMSDataPoint* getDataSet ( int i ) { return &dataPointList [i]; }
	void printSpecInfoDelimited ( std::ostream& os, int i ) const { dataPointList [i].printSpecInfoDelimited ( os ); }
	std::string getSpecID ( int i ) const { return dataPointList [i].getSpecID (); }
	std::string getSpot ( int i ) const { return dataPointList [i].getSpot (); }
	void putProductCGI ( std::ostream& os, int i ) const;
	static IntVector setChargeRange ( const ParameterList* params );
};
class MSMSQuanDataSetInfo : public DataSetInfo {
	virtual void getData ( const std::string& s, bool fileFlag, int fraction, const SpecID& specID, const std::string& version, int off = 0 ) {};	// stub
public:
	MSMSQuanDataSetInfo ( const std::string& filepath );
	~MSMSQuanDataSetInfo ();
	int getNumDataSets () const { return 0; }																// stub
	bool getMSMSQuanData ( XYData& xyData, const SpecID& specID, double startMass, double endMass, const std::string& version );
};

class UpdatingJavascriptMessage;

class MSMSPeakListDataSetInfo : public DataSetInfo {
	virtual void getData ( const std::string& s, bool fileFlag, int fraction, const SpecID& specID, const std::string& version, int off = 0 ) {};	// stub
public:
	MSMSPeakListDataSetInfo ( const std::string& filepath );
	~MSMSPeakListDataSetInfo ();
	int getNumDataSets () const { return 0; }																// stub
	bool writePeakList ( std::ostream& os, const SpecID& specID, const std::string& version );
	void getData ( MSMSDataPointVector& msmsDataPointList, const SpecID& specID, const std::string& version, UpdatingJavascriptMessage* ujm );
};
class MSMSPeakListDataFilterInfo : public DataSetInfo {
	virtual void getData ( const std::string& s, bool fileFlag, int fraction, const SpecID& specID, const std::string& version, int off = 0 ) {};	// stub
public:
	MSMSPeakListDataFilterInfo ( const std::string& filepath );
	~MSMSPeakListDataFilterInfo ();
	int getNumDataSets () const { return 0; }																// stub
	void writePeakList ( std::ostream& os );
	void readPeakList ( MSMSDataPointVector& msmsDataPointList, int index );
	void rewind ();
	bool isEOF ();
};

#endif /* ! __lu_df_info_h */
