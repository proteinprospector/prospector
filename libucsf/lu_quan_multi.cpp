/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_quan_multi.cpp                                             *
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
#ifdef RAW_DATA
#include <stdexcept>
#include <lg_string.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lu_quan_multi.h>
#include <lu_app_gr.h>
#include <lu_delim.h>
#include <lu_table.h>
#include <lu_iso.h>
#include <lu_inst.h>
#include <lu_getfil.h>
#include <lu_mass_elem.h>
#include <lu_param_list.h>
#include <lu_mass.h>
using std::string;
using std::ostream;
using std::cout;
using std::ostringstream;
using std::make_pair;
using std::find;
using std::getline;
using std::istringstream;
using std::runtime_error;

namespace {
	string FROM_FORMULAE = "From formulae";
	string NO_CORRECTION = "No correction";
	string DEFAULT = "Default";
}

QuanMSMSXMLData::QuanMSMSXMLData ()
{
	try {
		parseXMLFromFile ( MsparamsDir::instance ().getParamPath ( "quan_msms.xml" ) );
	}
	catch ( runtime_error e ) {
		ErrorHandler::genError ()->error ( e );
	}
}
QuanMSMSXMLData::~QuanMSMSXMLData () {}
QuanMSMSXMLData& QuanMSMSXMLData::instance ()
{
	static QuanMSMSXMLData quanMSMSXMLData;
	return quanMSMSXMLData;
}
void QuanMSMSXMLData::startElement ( const char* name, const char** attributes )
{
	if ( !strcmp ( name, "quan_msms_type" ) ) {
		formulae.resize ( 0 );
		masses.resize ( 0 );
		quanPeaks.resize ( 0 );
		numReporterIons = 0;
		quanName = "";
		quanNameFlag = false;
	}
	if ( !strcmp ( name, "name" ) ) {
		quanName = "";
		quanNameFlag = true;
	}
	if ( !strcmp ( name, "reporter_ion" ) ) {
		string formula;
		bool flag = getAttributeValue ( attributes, "formula", formula );
		if ( flag ) {
			formulae.push_back ( formula );
			masses.push_back ( formula_to_monoisotopic_mass ( formula.c_str () ) - ELECTRON_REST_MASS );
		}
		else {
			formulae.push_back ( "" );
			double m;
			flag = getAttributeValue ( attributes, "mass", m );
			if ( flag ) {
				masses.push_back ( m );
			}
			else {
				ErrorHandler::genError ()->error ( "File quan_msms.xml has an invalid format.\n" );
			}
		}
		bool qPeak = true;
		getAttributeValue ( attributes, "quan_peak", qPeak );
		if ( qPeak ) numReporterIons++;
		quanPeaks.push_back ( qPeak );
	}
}
void QuanMSMSXMLData::endElement ( const char* name )
{
	if ( !strcmp ( name, "quan_msms_type" ) ) {
		string fName = MsparamsDir::instance ().getParamPath ( quanName + ".txt" );
		StringVector purityNames;
		if ( genFileExists ( fName ) ) {
			GenCommentedIFStream ifs ( fName );
			string line;
			while ( ifs.getUncommentedLine ( line ) ) {
				purityNames.push_back ( line );
				ParameterList params ( ifs );
			}
		}
		QuanMSMSInfo qmmi ( formulae, masses, quanPeaks, numReporterIons, purityNames );
		quanMSMSInfo [quanName] = qmmi;
	}
	if ( !strcmp ( name, "name" ) ) {
		quanNameFlag = false;
	}
}
void QuanMSMSXMLData::characterDataHandler ( const char* str, int len )
{
	if ( quanNameFlag ) quanName.append ( str, len );
}
QuanMSMSInfo QuanMSMSXMLData::getQuanMSMSInfo ( const string& name ) const
{
	MapStringQuanMSMSInfoConstIterator cur = quanMSMSInfo.find ( name );
	QuanMSMSInfo defaultVal; 
	if ( cur != quanMSMSInfo.end () )
		return (*cur).second;
	else
		return defaultVal;
}
int QuanMSMSXMLData::getNumReporterIons ( const string& name ) const
{
	return getQuanMSMSInfo ( name ).numReporterIons;
}
int QuanMSMSXMLData::getNumQuanPeaks ( const string& name ) const
{
	return getQuanMSMSInfo ( name ).quanPeaks.size ();
}
DoubleVectorVector QuanMSMSXMLData::getFormulaPurityCoefficients ( const string& name ) const
{
	DoubleVectorVector dvv;
	StringVector f = getQuanMSMSInfo ( name ).formulae;
	DoubleVector m = getQuanMSMSInfo ( name ).masses;
	for ( int i = 0 ; i < f.size () ; i++ ) {
		DoubleVector dv;
		if ( f [i] != "" ) {
			IsotopePeakStats ips ( f [i], 1 );
			int index = ips.getProbabilityIndexForMass ( m [i] );
			if ( index >= 2 ) dv.push_back ( ips.getProbability ( index - 2 ) );
			else dv.push_back ( 0.0 );
			if ( index >= 1 ) dv.push_back ( ips.getProbability ( index - 1 ) );
			else dv.push_back ( 0.0 );
			dv.push_back ( ips.getProbability ( index + 1 ) );
			dv.push_back ( ips.getProbability ( index + 2 ) );
		}
		else {									// No formulae so assume no correction
			for ( int j = 0 ; j < 4 ; j++ ) {
				dv.push_back ( 0.0 );
			}
		}
		dvv.push_back ( dv );
	}
	return dvv;
}
StringVector QuanMSMSXMLData::getQuanMSMSNames () const
{
	StringVector sv;
	for ( MapStringQuanMSMSInfoConstIterator i = quanMSMSInfo.begin () ; i != quanMSMSInfo.end () ; i++ ) {
		sv.push_back ( i->first );
	}
	return sv;
}
StringVector QuanMSMSXMLData::getPurityNames () const
{
	StringVector sv;
	sv.push_back ( FROM_FORMULAE );
	sv.push_back ( DEFAULT );
	sv.push_back ( NO_CORRECTION );
	for ( MapStringQuanMSMSInfoConstIterator i = quanMSMSInfo.begin () ; i != quanMSMSInfo.end () ; i++ ) {
		string n = i->first;
		StringVector pn = i->second.purityNames;
		for ( int j = 0 ; j < pn.size () ; j++ ) {
			if ( pn [j] != DEFAULT ) sv.push_back ( n + " " + pn [j] );
		}
	}
	return sv;
}
DoubleVector QuantitationMulti::iTRAQMass;
BoolDeque QuantitationMulti::iTRAQPks;
int QuantitationMulti::iTRAQRefPk = 0;
int QuantitationMulti::numQuanPks = 0;
int QuantitationMulti::numReporterPks = 0;
const PurityCorrection* QuantitationMulti::purityCorrection = 0;
QuantitationMulti::QuantitationMulti () :
	QuantitationData ( numQuanPks ),
	intensity ( numQuanPks, 0.0 ),
	area ( numQuanPks, 0.0 )
{
}
void QuantitationMulti::init ( const string& quanType )
{
	if ( numQuanPks == 0 ) {
		QuanMSMSInfo qmmi = QuanMSMSXMLData::instance ().getQuanMSMSInfo ( quanType );
		iTRAQRefPk = 0;
		iTRAQMass = qmmi.masses;
		iTRAQPks = qmmi.quanPeaks;
		numReporterPks = qmmi.numReporterIons;
		numQuanPks = iTRAQMass.size ();
	}
}
string QuantitationMulti::getMassString ( int i )
{
	return gen_itoa ( (int)iTRAQMass [i] );
}
string QuantitationMulti::getRatioString ( int i )
{
	return gen_itoa ( (int)iTRAQMass [iTRAQRefPk] ) + "/" + getMassString ( i );
}
int QuantitationMulti::getColspan ( bool normal )
{
	int colspan = 0;
	if ( normal ) {
		if ( reportPeakSNR )		colspan += numReporterPks;
		if ( reportPeakResolution )	colspan += numReporterPks;
		if ( reportPeakFWHM )		colspan += numReporterPks;
	}
	if ( reportPeakIntensity )					colspan += numReporterPks;
	if ( reportPeakArea )						colspan += numReporterPks;
	if ( reportActualLightHeavyIntensityRatio )	colspan += numReporterPks - 1;
	if ( reportActualLightHeavyAreaRatio )		colspan += numReporterPks - 1;
	return colspan;
}
void QuantitationMulti::printHTMLHeader ( ostream& os, const string& styleID, bool normal )
{
	if ( reportPeakIntensity )	for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) tableHeader ( os, "Int " + getMassString ( i ), styleID );
	if ( normal ) {
		if ( reportPeakSNR )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) tableHeader ( os, "SNR " + getMassString ( i ), styleID );
		if ( reportPeakResolution )	for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) tableHeader ( os, "Resolution " + getMassString ( i ), styleID );
		if ( reportPeakFWHM )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) tableHeader ( os, "FWHM " + getMassString ( i ), styleID );
	}
	if ( reportPeakArea )						for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) tableHeader ( os, "Area " + getMassString ( i ), styleID );
	if ( reportActualLightHeavyIntensityRatio )	for ( int i = 0 ; i < numQuanPks ; i++ ) if ( i != iTRAQRefPk && iTRAQPks [i] ) tableHeader ( os, getRatioString ( i ) + " Int", styleID );
	if ( reportActualLightHeavyAreaRatio )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( i != iTRAQRefPk && iTRAQPks [i] ) tableHeader ( os, getRatioString ( i ) + " Area", styleID );
}
void QuantitationMulti::printHTMLBlankLine ( ostream& os, const string& styleID, bool normal )
{
	if ( reportPeakIntensity )	tableEmptyNCells ( os, numReporterPks, styleID );
	if ( normal ) {
		if ( reportPeakSNR )		tableEmptyNCells ( os, numReporterPks, styleID );
		if ( reportPeakResolution )	tableEmptyNCells ( os, numReporterPks, styleID );
		if ( reportPeakFWHM )		tableEmptyNCells ( os, numReporterPks, styleID );
	}
	if ( reportPeakArea )		tableEmptyNCells ( os, numReporterPks, styleID );
	if ( reportActualLightHeavyIntensityRatio )	tableEmptyNCells ( os, numReporterPks-1, styleID );
	if ( reportActualLightHeavyAreaRatio )		tableEmptyNCells ( os, numReporterPks-1, styleID );
}
void QuantitationMulti::printDelimitedHeader ( ostream& os, bool normal )
{
	if ( reportPeakIntensity )	for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) delimitedHeader ( os, "Int " + getMassString ( i ) );
	if ( normal ) {
		if ( reportPeakSNR )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) delimitedHeader ( os, "SNR " + getMassString ( i ) );
		if ( reportPeakResolution )	for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) delimitedHeader ( os, "Resolution " + getMassString ( i ) );
		if ( reportPeakFWHM )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) delimitedHeader ( os, "FWHM " + getMassString ( i ) );
	}
	if ( reportPeakArea )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) delimitedHeader ( os, "Area " + getMassString ( i ) );
	if ( reportActualLightHeavyIntensityRatio )	for ( int i = 0 ; i < numQuanPks ; i++ ) if ( i != iTRAQRefPk && iTRAQPks [i] ) delimitedHeader ( os, getRatioString ( i ) + " Int" );
	if ( reportActualLightHeavyAreaRatio )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( i != iTRAQRefPk && iTRAQPks [i] ) delimitedHeader ( os, getRatioString ( i ) + " Area" );
}
void QuantitationMulti::printDelimitedBlankLine ( ostream& os, bool normal )
{
	if ( reportPeakIntensity )	delimitedEmptyNCells ( os, numReporterPks );
	if ( normal ) {
		if ( reportPeakSNR )		delimitedEmptyNCells ( os, numReporterPks );
		if ( reportPeakResolution )	delimitedEmptyNCells ( os, numReporterPks );
		if ( reportPeakFWHM )		delimitedEmptyNCells ( os, numReporterPks );
	}
	if ( reportPeakArea )		delimitedEmptyNCells ( os, numReporterPks );
	if ( reportActualLightHeavyIntensityRatio )	delimitedEmptyNCells ( os, numReporterPks-1 );
	if ( reportActualLightHeavyAreaRatio )		delimitedEmptyNCells ( os, numReporterPks-1 );
}
bool QuantitationMulti::outputQuanResults ( ostream& os, const string& searchName, int numRepeats, bool area ) const
{
	bool flag = false;
	for ( int i = 0 ; i < numQuanPks ; i++ ) {
		if ( i != iTRAQRefPk && iTRAQPks [i] ) {
			ostringstream s;
			s << searchName << " " << getRatioString ( i );
			if ( area ) {
				if ( outputQuanRatio ( os, s.str (), lightHeavyAreaRatio [0][i], lightHeavyAreaRatioType [0][i], numRepeats ) ) {
					flag = true;
				}
			}
			else {
				if ( outputQuanRatio ( os, s.str (), lightHeavyIntRatio [0][i], lightHeavyIntRatioType [0][i], numRepeats ) ) {
					flag = true;
				}
			}
		}
	}
	return flag;
}
DoubleVector QuantitationMulti::getAreaRatios () const
{
	DoubleVector dv;
	for ( int i = 0 ; i < numQuanPks ; i++ ) {
		if ( i != iTRAQRefPk && iTRAQPks [i] ) {
			dv.push_back ( getQuanRatio ( lightHeavyAreaRatio [0][i], lightHeavyAreaRatioType [0][i] ) );
		}
	}
	return dv;
}
DoubleVector QuantitationMulti::getIntensityRatios () const
{
	DoubleVector dv;
	for ( int i = 0 ; i < numQuanPks ; i++ ) {
		if ( i != iTRAQRefPk && iTRAQPks [i] ) {
			dv.push_back ( getQuanRatio ( lightHeavyIntRatio [0][i], lightHeavyIntRatioType [0][i] ) );
		}
	}
	return dv;
}
QuantitationMultiNormal::QuantitationMultiNormal ( const XYData& xydata, double res ) :
	ok ( numQuanPks, false ),
	mOverZ ( numQuanPks, 0.0 ),
	snr ( numQuanPks, 0.0 ),
	fwhm ( numQuanPks, 0.0 ),
	resolution ( numQuanPks, 0.0 )
{
	init ( xydata, res );
}
void QuantitationMultiNormal::init ( const XYData& xyData, double res )
{
	if ( xyData.size () ) {
		PeakFit::getNoise ( xyData, noiseMean, noiseStDev );

		double minimumMaxMass = floor ( iTRAQMass [numQuanPks-1] ) + 1.0;
		if ( xyData.maxX () < minimumMaxMass ) const_cast <XYData&> (xyData).add ( minimumMaxMass, 1 ); // Deal with the case where the data stops too close to the high mass peak
		ColouredGraphData* graphData = 0;
		bool graphFlag = PeakFit::getGraphs ();
		if ( graphFlag ) graphData = new ColouredGraphData ( xyData );
		DoubleVectorVectorVector coeff;
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			coeff.push_back ( PeakFit::getCoefficients ( xyData, graphData, iTRAQMass [i], 1, 1, res, noiseMean, noiseStDev, snrThreshold ) );
		}
		if ( graphFlag ) {
			PeakFit::drawGraph ( *graphData, false );
			delete graphData;
		}
		using nr::sqrtTwo;
		using nr::sqrtPi;

		for ( int j = 0 ; j < numQuanPks ; j++ ) {
			ok [j] = !coeff [j].empty ();
		}
		for ( int k = 0 ; k < numQuanPks ; k++ ) {
			if ( ok [k] ) {
				intensity [k] = coeff [k][0][0];
				mOverZ [k] = coeff [k][0][1];
				double width = coeff [k][0][2]; 
				snr [k] = intensity [k] / noiseStDev;
				fwhm [k] = 2.35 * width / sqrtTwo;
				resolution [k] = mOverZ [k] / fwhm [k];
				area [k] = intensity [k] * width * sqrtPi;
			}
		}
		if ( purityCorrection ) purityCorrection->correction ( intensity );
		if ( purityCorrection ) purityCorrection->correction ( area );
		for ( int x = 0 ; x < numQuanPks ; x++ ) {
			if ( ok [iTRAQRefPk] && ok [x] ) {
				if ( snr [iTRAQRefPk] < snrThreshold && snr [x] < snrThreshold ) {
					lightHeavyIntRatioType [0][x] = RATIO_NOT_CALCULATED;
					lightHeavyAreaRatioType [0][x] = RATIO_NOT_CALCULATED;
				}
				else if ( snr [iTRAQRefPk] < snrThreshold ) {
					double lowIntensity = noiseStDev * snrThreshold;
					if ( lowIntensity > intensityThreshold && intensity [x] > intensityThreshold )
						setLTRatio ( lightHeavyIntRatioType [0][x], lightHeavyIntRatio [0][x], lowIntensity, intensity [x] );
					else
						setRatio ( lightHeavyIntRatioType [0][x], lightHeavyIntRatio [0][x], intensityThreshold, lowIntensity, intensity [x] );
					double lowArea = lowIntensity * coeff [x][0][2] * sqrtPi;	// Use width of other peak
					if ( lowArea > areaThreshold && area [x] > areaThreshold )
						setLTRatio ( lightHeavyAreaRatioType [0][x], lightHeavyAreaRatio [0][x], lowArea, area [x] );
					else
						setRatio ( lightHeavyAreaRatioType [0][x], lightHeavyAreaRatio [0][x], areaThreshold, lowArea, area [x] );
				}
				else if ( snr [x] < snrThreshold ) {
					double highIntensity = noiseStDev * snrThreshold;
					if ( highIntensity > intensityThreshold && intensity [iTRAQRefPk] > intensityThreshold )
						setGTRatio ( lightHeavyIntRatioType [0][x], lightHeavyIntRatio [0][x], intensity [iTRAQRefPk], highIntensity );
					else
						setRatio ( lightHeavyIntRatioType [0][x], lightHeavyIntRatio [0][x], intensityThreshold, intensity [iTRAQRefPk], highIntensity );
					double highArea = highIntensity * coeff [iTRAQRefPk][0][2] * sqrtPi;	// Use width of other peak
					if ( highArea > areaThreshold && area [iTRAQRefPk] > areaThreshold )
						setGTRatio ( lightHeavyAreaRatioType [0][x], lightHeavyAreaRatio [0][x], area [iTRAQRefPk], highArea );
					else
						setRatio ( lightHeavyAreaRatioType [0][x], lightHeavyAreaRatio [0][x], areaThreshold, area [iTRAQRefPk], highArea );
				}
				else {
					setRatio ( lightHeavyIntRatioType [0][x], lightHeavyIntRatio [0][x], intensityThreshold, intensity [iTRAQRefPk], intensity [x] );
					setRatio ( lightHeavyAreaRatioType [0][x], lightHeavyAreaRatio [0][x], areaThreshold, area [iTRAQRefPk], area [x] );
				}
			}
		}
	}
}
void QuantitationMultiNormal::printHTML ( ostream& os ) const
{
	string styleID1 = "sc_stripe_1";
	tableStart ( os, true );
		tableRowStart ( os );
			for ( int i = 0 ; i < numQuanPks ; i++ ) {
				if ( iTRAQPks [i] ) tableHeader ( os, "m/z " + getMassString ( i ), styleID1 );
			}
			printHTMLNoiseHeader ( os, styleID1 );
		tableRowEnd ( os );

		tableRowStart ( os );
			printHTMLMOverZ ( os, 0, styleID1 );
			printHTMLNoiseLine ( os, styleID1 );
		tableRowEnd ( os );
	tableEnd ( os );

	tableStart ( os, true );
		tableRowStart ( os );
			printHTMLHeader ( os, styleID1, true );
		tableRowEnd ( os );

		tableRowStart ( os );
			printHTMLLine ( os, 0, styleID1 );
		tableRowEnd ( os );
	tableEnd ( os );
}
void QuantitationMultiNormal::printHTMLMOverZ ( ostream& os, DoubleVectorVectorSizeType peakNumber, const string& styleID ) const
{
	for ( int i = 0 ; i < numQuanPks ; i++ ) {
		if ( iTRAQPks [i] ) {
			if ( ok [i] ) tableCell ( os, mOverZ [i], 4, false, styleID );
			else tableEmptyCell ( os );
		}
	}
}
void QuantitationMultiNormal::printHTMLLine ( ostream& os, DoubleVectorVectorSizeType peakNumber, const string& styleID ) const
{
	if ( reportPeakIntensity ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( iTRAQPks [i] ) {
				if ( ok [i] ) tableCellSigFig ( os, intensity [i], 4, false, styleID );
				else tableEmptyCell ( os, styleID );
			}
		}
	}
	if ( reportPeakSNR ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( iTRAQPks [i] ) {
				if ( ok [i] ) tableCellSigFig ( os, snr [i], 4, false, styleID );
				else tableEmptyCell ( os, styleID );
			}
		}
	}
	if ( reportPeakResolution ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( iTRAQPks [i] ) {
				if ( ok [i] ) tableCellSigFig ( os, resolution [i], 0, false, styleID );
				else tableEmptyCell ( os, styleID );
			}
		}
	}
	if ( reportPeakFWHM ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( iTRAQPks [i] ) {
				if ( ok [i] ) tableCellSigFig ( os, fwhm [i], 4, false, styleID );
				else tableEmptyCell ( os, styleID );
			}
		}
	}
	if ( reportPeakArea ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( iTRAQPks [i] ) {
				if ( ok [i] ) tableCellSigFig ( os, area [i], 4, false, styleID );
				else tableEmptyCell ( os, styleID );
			}
		}
	}
	if ( reportActualLightHeavyIntensityRatio )	{
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( i != iTRAQRefPk && iTRAQPks [i] ) {
				if ( ok [i] ) printHTMLRatio ( os, lightHeavyIntRatio [0][i], lightHeavyIntRatioType [0][i], styleID );
				else tableEmptyCell ( os, styleID );
			}
		}
	}
	if ( reportActualLightHeavyAreaRatio ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( i != iTRAQRefPk && iTRAQPks [i] ) {
				if ( ok [i] ) printHTMLRatio ( os, lightHeavyAreaRatio [0][i], lightHeavyAreaRatioType [0][i], styleID );
				else tableEmptyCell ( os, styleID );
			}
		}
	}
}
void QuantitationMultiNormal::printDelimitedLine ( ostream& os, DoubleVectorVectorSizeType peakNumber ) const
{
	if ( reportPeakIntensity ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( iTRAQPks [i] ) {
				if ( ok [i] ) delimitedCellSigFig ( os, intensity [i], 4 );
				else delimitedEmptyCell ( os );
			}
		}
	}
	if ( reportPeakSNR ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( iTRAQPks [i] ) {
				if ( ok [i] ) delimitedCellSigFig ( os, snr [i], 4 );
				else delimitedEmptyCell ( os );
			}
		}
	}
	if ( reportPeakResolution ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( iTRAQPks [i] ) {
				if ( ok [i] ) delimitedCellSigFig ( os, resolution [i], 0 );
				else delimitedEmptyCell ( os );
			}
		}
	}
	if ( reportPeakFWHM ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( iTRAQPks [i] ) {
				if ( ok [i] ) delimitedCellSigFig ( os, fwhm [i], 4 );
				else delimitedEmptyCell ( os );
			}
		}
	}
	if ( reportPeakArea ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( iTRAQPks [i] ) {
				if ( ok [i] ) delimitedCellSigFig ( os, area [i], 4 );
				else delimitedEmptyCell ( os );
			}
		}
	}
	if ( reportActualLightHeavyIntensityRatio )	{
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( i != iTRAQRefPk && iTRAQPks [i] ) {
				if ( ok [i] ) printDelimitedRatio ( os, lightHeavyIntRatio [0][i], lightHeavyIntRatioType [0][i] );
				else delimitedEmptyCell ( os );
			}
		}
	}
	if ( reportActualLightHeavyAreaRatio ) {
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			if ( i != iTRAQRefPk && iTRAQPks [i] ) {
				if ( ok [i] ) printDelimitedRatio ( os, lightHeavyAreaRatio [0][i], lightHeavyAreaRatioType [0][i] );
				else delimitedEmptyCell ( os );
			}
		}
	}
}
double QuantitationMultiMassWindow::reporterIonWindow = 0.4;
QuantitationMultiMassWindow::QuantitationMultiMassWindow ( const XYData& xydata )
{
	init ( xydata );
}
void QuantitationMultiMassWindow::init ( const XYData& xyData )
{
	if ( xyData.size () ) {
		double quanTolerance = reporterIonWindow / 2.0;
		noiseMean = 0.0;
		noiseStDev = 0.0;

		GraphData graphData ( xyData );
		if ( PeakFit::getGraphs () ) {
			SpectrumGraph s ( "pr_graph.par.txt" );
			s.drawGraph ( cout, graphData );
		}
		for ( int i = 0 ; i < numQuanPks ; i++ ) {
			try {
				intensity [i] = xyData.getMaxYInTolRange ( iTRAQMass [i], quanTolerance );
			}
			catch ( lNrecEmptyXYDataRange) {}	// Not bothered already initialized to zero
			area [i] = xyData.getSumOfYInTolRange ( iTRAQMass [i], quanTolerance );
		}
		if ( purityCorrection ) purityCorrection->correction ( intensity );
		if ( purityCorrection ) purityCorrection->correction ( area );
		for ( int j = 0 ; j < numQuanPks ; j++ ) {
			setRatio ( lightHeavyIntRatioType [0][j], lightHeavyIntRatio [0][j], intensityThreshold, intensity [iTRAQRefPk], intensity [j] );
			setRatio ( lightHeavyAreaRatioType [0][j], lightHeavyAreaRatio [0][j], areaThreshold, area [iTRAQRefPk], area [j] );
		}
	}
}
void QuantitationMultiMassWindow::printHTML ( ostream& os ) const
{
	string styleID1 = "sc_stripe_1";
	tableStart ( os, true );
		tableRowStart ( os );
			printHTMLHeader ( os, styleID1, false );
		tableRowEnd ( os );

		tableRowStart ( os );
			printHTMLLine ( os, 0, styleID1 );
		tableRowEnd ( os );
	tableEnd ( os );
}
void QuantitationMultiMassWindow::printHTMLLine ( ostream& os, DoubleVectorVectorSizeType peakNumber, const string& styleID ) const
{
	if ( reportPeakIntensity )					for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) tableCellSigFig ( os, intensity [i], 4, false, styleID );
	if ( peakNumber == -1 ) {	// Joint report
		if ( reportPeakSNR )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) tableEmptyCell ( os, styleID );
		if ( reportPeakResolution )	for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) tableEmptyCell ( os, styleID );
		if ( reportPeakFWHM )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) tableEmptyCell ( os, styleID );
	}
	if ( reportPeakArea )						for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) tableCellSigFig ( os, area [i], 4, false, styleID );
	if ( reportActualLightHeavyIntensityRatio )	for ( int i = 0 ; i < numQuanPks ; i++ ) if ( i != iTRAQRefPk && iTRAQPks [i] ) printHTMLRatio ( os, lightHeavyIntRatio [0][i], lightHeavyIntRatioType [0][i], styleID );
	if ( reportActualLightHeavyAreaRatio )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( i != iTRAQRefPk && iTRAQPks [i] ) printHTMLRatio ( os, lightHeavyAreaRatio [0][i], lightHeavyAreaRatioType [0][i], styleID );
}
void QuantitationMultiMassWindow::printDelimitedLine ( ostream& os, DoubleVectorVectorSizeType peakNumber ) const
{
	if ( reportPeakIntensity )					for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) delimitedCellSigFig ( os, intensity [i], 4 );
	if ( peakNumber == -1 ) {	// Joint report
		if ( reportPeakSNR )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) delimitedEmptyCell ( os );
		if ( reportPeakResolution )	for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) delimitedEmptyCell ( os );
		if ( reportPeakFWHM )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) delimitedEmptyCell ( os );
	}
	if ( reportPeakArea )						for ( int i = 0 ; i < numQuanPks ; i++ ) if ( iTRAQPks [i] ) delimitedCellSigFig ( os, area [i], 4 );
	if ( reportActualLightHeavyIntensityRatio )	for ( int i = 0 ; i < numQuanPks ; i++ ) if ( i != iTRAQRefPk && iTRAQPks [i] ) printDelimitedRatio ( os, lightHeavyIntRatio [0][i], lightHeavyIntRatioType [0][i] );
	if ( reportActualLightHeavyAreaRatio )		for ( int i = 0 ; i < numQuanPks ; i++ ) if ( i != iTRAQRefPk && iTRAQPks [i] ) printDelimitedRatio ( os, lightHeavyAreaRatio [0][i], lightHeavyAreaRatioType [0][i] );
}
PurityCorrection::PurityCorrection ( const ParameterList* params )
{
	string quanType = params->getStringValue ( "quan_type", "" );
	if ( isQuanMSMS ( quanType ) ) {	// This is an MSMS quantitation type
		string purityName = params->getStringValue ( "purity_correction", "" );
		string purityFile = quanType + ".txt";
		string purityPath = MsparamsDir::instance ().getParamPath ( purityFile );
		if ( purityName == FROM_FORMULAE ) {
			numQuanPeaks = QuanMSMSXMLData::instance ().getNumQuanPeaks ( quanType );
			matrix = QuanMSMSXMLData::instance ().getFormulaPurityCoefficients ( quanType );
		}
		else if ( purityName == NO_CORRECTION ) {
			numQuanPeaks = QuanMSMSXMLData::instance ().getNumQuanPeaks ( quanType );
			matrix.resize ( numQuanPeaks );
			for ( int i = 0 ; i < numQuanPeaks ; i++ ) {
				for ( int j = 0 ; j < 4 ; j++ ) {
					matrix [i].push_back ( 0.0 );
				}
			}
		}
		else {
			if ( purityName != DEFAULT ) {
				if ( isPrefix ( purityName, quanType ) ) {
					purityName = purityName.substr ( quanType.length () + 1 );
				}
				else {
					ErrorHandler::genError ()->error ( "Invalid purity name.\n" );
				}
			}
			GenIFStream fromFile ( purityPath );
			string line;
			bool flag = false;
			while ( getline ( fromFile, line ) ) {
				if ( line == purityName ) {
					ParameterList p ( fromFile );
					StringVector sv = p.getNameList ();
					numQuanPeaks = p.size ();
					matrix.resize ( numQuanPeaks );
					for ( int i = 0 ; i < numQuanPeaks ; i++ ) {
						istringstream ist ( p.getStringValue ( sv [i] ) );
						double num;
						while ( ist >> num ) {
							matrix [i].push_back ( num / 100.0 );
						}
					}
					flag = true;
					break;
				}
			}
			if ( !flag ) {
				ErrorHandler::genError ()->error ( "Invalid purity name.\n" );
			}
		}
		a = nrmatrix ( 1, numQuanPeaks, 1, numQuanPeaks );
		loadA ();
		b = nrmatrix ( 1, numQuanPeaks, 1, 1 );
	}
	else numQuanPeaks = 0;
}
PurityCorrection::~PurityCorrection ()
{
	if ( numQuanPeaks ) {
		free_nrmatrix ( a, 1, numQuanPeaks, 1, numQuanPeaks );
		free_nrmatrix ( b, 1, numQuanPeaks, 1, 1 );
	}
}
void PurityCorrection::loadA () const
{
	for ( int i = 1 ; i <= numQuanPeaks ; i++ ) {
		for ( int j = 1 ; j <= numQuanPeaks ; j++ ) {
			a [i][j] = 0.0;
			if ( i == j+2 ) a [i][j] = matrix [j-1][3];
			if ( i == j+1 ) a [i][j] = matrix [j-1][2];
			if ( i == j ) a [i][j] = 1.0 - ( matrix [i-1][0] + matrix [i-1][1] + matrix [i-1][2] + matrix [i-1][3] );
			if ( i == j-1 ) a [i][j] = matrix [j-1][1];
			if ( i == j-2 ) a [i][j] = matrix [j-1][0];
		}
	}
}
void PurityCorrection::correction ( DoubleVector& dv ) const
{
	loadA ();			// a must be reloaded every time as it is overwritten by gaussj
	for ( int i = 0 ; i < numQuanPeaks ; i++ ) {
		b [i+1][1] = dv [i];
	}
	gaussj ( a, numQuanPeaks, b, 1 );
	for ( int j = 0 ; j < numQuanPeaks ; j++ ) {
		dv [j] = b [j+1][1];
	}
}
StringVector getQuanMSMSNames ()
{
	return QuanMSMSXMLData::instance ().getQuanMSMSNames ();
}
StringVector getPurityNames ()
{
	return QuanMSMSXMLData::instance ().getPurityNames ();
}
bool isQuanMSMS ( const string& quanType )
{
	StringVector qmmn = getQuanMSMSNames ();
	StringVectorIterator svi = find ( qmmn.begin (), qmmn.end (), quanType );
	return svi != qmmn.end ();	// This is an MSMS quantitation type
}
int getNumReporterIons ( const string& name )
{
	return QuanMSMSXMLData::instance ().getNumReporterIons ( name );
}
#endif
