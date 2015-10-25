/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_quan_multi.h                                               *
*                                                                             *
*  Created    : October 7th 2004                                              *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_quan_multi_h
#define __lu_quan_multi_h

#include <nr.h>
#include <lgen_xml.h>
#include <lu_parent.h>

class PurityCorrection {
	int numQuanPeaks;
	DoubleVectorVector matrix;
	double** a;
	double** b;
public:
	PurityCorrection ( const ParameterList* params );
	~PurityCorrection ();
	void loadA () const;
	void correction ( DoubleVector& dv ) const;
};

class QuantitationMulti : public QuantitationData {
protected:
	DoubleVector intensity;
	DoubleVector area;
	static DoubleVector iTRAQMass;
	static BoolDeque iTRAQPks;
	static int iTRAQRefPk;
	static int numQuanPks;
	static int numReporterPks;
	static const PurityCorrection* purityCorrection;
public:
	QuantitationMulti ();
	static void setPurityCorrection ( const PurityCorrection* pc )
		{ purityCorrection = pc; }
	static double getStartDisplayMass ()
	{
		return floor ( iTRAQMass [0] - 1.2 );
	}
	static double getEndDisplayMass ()
	{
		return floor ( iTRAQMass [numQuanPks-1] + 2.0 );
	}
	static int getColspan ( bool normal );
	static void printHTMLHeader ( std::ostream& os, const std::string& styleID, bool normal );
	static void printHTMLBlankLine ( std::ostream& os, const std::string& styleID, bool normal );
	static void printDelimitedHeader ( std::ostream& os, bool normal );
	static void printDelimitedBlankLine ( std::ostream& os, bool normal );
	static double getStartMass () { return getStartDisplayMass () + 0.5; }
	static double getEndMass () { return getEndDisplayMass () - 1.5; }
	bool outputQuanResults ( std::ostream& os, const std::string& searchName, int numRepeats, bool area ) const;
	DoubleVector getAreaRatios () const;
	DoubleVector getIntensityRatios () const;
	static void init ( const std::string& quanType );
	static std::string getMassString ( int i );
	static std::string getRatioString ( int i );
	static std::string getProteinRatioString ( int ind )
	{
		for ( int i = 0, j = 0 ; i < numQuanPks ; i++ ) {
			if ( i != iTRAQRefPk && iTRAQPks [i] ) {
				if ( j == ind ) return getRatioString ( i );
				j++;
			}
		}
		return "";
	}
};
	
class QuantitationMultiNormal : public QuantitationMulti {
	BoolDeque ok;
	DoubleVector mOverZ;
	DoubleVector snr;
	DoubleVector fwhm;
	DoubleVector resolution;

	void init ( const XYData& xyData, double resolution );
	void printHTMLMOverZ ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const;
public:
	QuantitationMultiNormal ( const XYData& xydata, double resolution );
	void printHTML ( std::ostream& os ) const;

	void printHTMLLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const;
	void printDelimitedLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber ) const;
};

class QuantitationMultiMassWindow : public QuantitationMulti {
	static double reporterIonWindow;
	void init ( const XYData& xyData );
	void printHTMLMOverZ ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const {}
public:
	QuantitationMultiMassWindow ( const XYData& xydata );
	void printHTML ( std::ostream& os ) const;

	void printHTMLLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber, const std::string& styleID ) const;
	void printDelimitedLine ( std::ostream& os, DoubleVectorVectorSizeType peakNumber ) const;
	static void setReporterIonWindow ( double riw )
		{ reporterIonWindow = riw; }
};

StringVector getQuanMSMSNames ();
StringVector getPurityNames ();
bool isQuanMSMS ( const std::string& quanType );
int getNumReporterIons ( const std::string& name );

struct QuanMSMSInfo {
	StringVector formulae;
	DoubleVector masses;
	BoolDeque quanPeaks;
	int numReporterIons;
	StringVector purityNames;
	QuanMSMSInfo () {}
	QuanMSMSInfo ( const StringVector& formulae, const DoubleVector& masses, const BoolDeque& quanPeaks, int numReporterIons, const StringVector& purityNames ) :
		formulae ( formulae ),
		masses ( masses ),
		quanPeaks ( quanPeaks ),
		numReporterIons ( numReporterIons ),
		purityNames ( purityNames ) {}
};
typedef std::map <std::string, QuanMSMSInfo> MapStringQuanMSMSInfo;
typedef MapStringQuanMSMSInfo::const_iterator MapStringQuanMSMSInfoConstIterator;

class QuanMSMSXMLData : public PPExpat {
	std::string quanName;
	bool quanNameFlag;
	StringVector formulae;
	DoubleVector masses;
	BoolDeque quanPeaks;
	int numReporterIons;
	MapStringQuanMSMSInfo quanMSMSInfo;
	void startElement ( const char* name, const char** attributes );
	void endElement ( const char* name );
	void characterDataHandler ( const char* str, int len );
	QuanMSMSXMLData ();
public:
	~QuanMSMSXMLData ();
	static QuanMSMSXMLData& instance ();
	QuanMSMSInfo getQuanMSMSInfo ( const std::string& name ) const;
	StringVector getQuanMSMSNames () const;
	StringVector getPurityNames () const;
	int getNumReporterIons ( const std::string& name ) const;
	int getNumQuanPeaks ( const std::string& name ) const;
	DoubleVectorVector getFormulaPurityCoefficients ( const std::string& name ) const;
};

#endif /* ! __lu_quan_multi_h */
