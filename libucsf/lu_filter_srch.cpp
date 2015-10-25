/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_filter_srch.cpp                                            *
*                                                                             *
*  Created    : September 3rd 2012                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2012-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <iostream>
#include <lg_io.h>
#include <lgen_error.h>
#include <lgen_file.h>
#include <lgen_uncompress.h>
#include <lu_df_info.h>
#include <lu_charge.h>
#include <lu_filter_srch.h>
#include <lu_html.h>
#include <lu_mass.h>
#include <lu_mass_conv.h>
#include <lu_mass_elem.h>
#include <lu_repository.h>

using std::ostream;
using std::cout;
using std::string;
using std::endl;
using std::count;

MSFilterSearch::MSFilterSearch ( const MSFilterParameters& params ) :
	MSProgram ( params ),
	filterParams ( params )
{
	init_html ( cout, "MS-Filter Report" );
	processPeakListFile ();
	string keepOrRemove = params.getKeepOrRemove ();
	if ( keepOrRemove == "keep" ) {
		filterPeakLists ( true );
	}
	else if ( keepOrRemove == "remove" ) {
		filterPeakLists ( false );
	}
	else {
		filterPeakLists ( true, "-matching" );
		filterPeakLists ( false, "-non-matching" );
	}
	genUnlinkDirectory ( genDirectoryFromPath ( inPeakListFPath ) );
}
MSFilterSearch::~MSFilterSearch ()
{
}
void MSFilterSearch::printBodyHTML ( ostream& os )
{
}
void MSFilterSearch::processPeakListFile ()
{
	inPeakListFPath = filterParams.getInPeakListFpath ();
	if ( inPeakListFPath.empty () ) return;
	outPeakListFPath = filterParams.getOutPeakListFpath ();
	outPeakListURL = filterParams.getOutPeakListURL ();
	archiveName = filterParams.getArchiveName ();
	PPProject ppp ( inPeakListFPath );
	if ( ppp.initialised () ) {
		peakListFractionNames = ppp.getFractionNames ();
		peakListCentroidFiles = ppp.getCentroidFiles ();
	}
	else {
		string err = ppp.getErrMessage ();
		genUnlinkDirectory ( inPeakListFPath );
		if ( err.empty () )
			ErrorHandler::genError ()->error ( "Upload file has an invalid format.\n" );
		else
			ErrorHandler::genError ()->error ( err + "\n" );
	}
}
void MSFilterSearch::filterPeakLists ( bool keepFlag, const string& aSuffix )
{
	StringVector deletePaths;
	cout << "<p>" << endl;
	for ( int i = 0 ; i < peakListCentroidFiles.size () ; i++ ) {
		PairIntInt pii ( 0, 0 );
		string file = peakListCentroidFiles [i];
		string inFile = inPeakListFPath + SLASH + file;
		string outFile = outPeakListFPath + aSuffix + SLASH + file;
		deletePaths.push_back ( outFile );

		MSMSPeakListDataFilterInfo pkList ( inFile );
		GenOFStream ost ( outFile, std::ios_base::out | std::ios_base::app );

		for ( int index = 1 ; ; index++ ) {
			MSMSDataPointVector msmsDataPointList;
			pkList.readPeakList ( msmsDataPointList, index );
			bool flag = true;
			if ( !msmsDataPointList.empty () ) {
				pii.second++;
				Tolerance* parentTol = filterParams.getParentMassTolerance ();
				double systematicError = filterParams.getSystematicError ();
				msmsDataPointList [0].calibrate ( parentTol, systematicError );
				flag = checkMPlusHRange ( msmsDataPointList [0] );
				if ( flag ) {
					flag = checkCharge ( msmsDataPointList [0] );
				}
				if ( flag ) {
					flag = checkNeutralLoss ( msmsDataPointList [0] );
				}
				if ( flag ) {
					flag = checkFragmentMZs ( msmsDataPointList [0] );
				}
				if ( ( flag && keepFlag ) || ( !flag && !keepFlag ) ) {
					pii.first++;
					pkList.rewind ();
					pkList.writePeakList ( ost );
				}
			}
			if ( pkList.isEOF () ) {
				cout << file << " " << pii.first << "/" << pii.second << " retained<br />" << endl;
				break;
			}
		}
	}
	cout << "</p>" << endl;
	bool flag = gen7zaCreate ( outPeakListFPath + aSuffix + SLASH + archiveName + aSuffix, outPeakListFPath + aSuffix + SLASH + "*", "zip" );
	genUnlink ( deletePaths );
	if ( flag ) {
		cout << "<a href=\"" << outPeakListURL + aSuffix + "/" + archiveName + aSuffix + ".zip" << "\">";
		cout << "Archive file (" << archiveName + aSuffix << ".zip)</a><br />" << endl;
	}
}
bool MSFilterSearch::checkMPlusHRange ( MSMSDataPoint& mmdp )
{
	bool fullMPlusHRange = filterParams.getFullMPlusHRange ();
	if ( !fullMPlusHRange ) {
		double lowMPlusH = filterParams.getLowMPlusH ();
		double highMPlusH = filterParams.getHighMPlusH ();
		double parentMPlusH = mmdp.getParentMPlusH ();
		if ( parentMPlusH < lowMPlusH || parentMPlusH > highMPlusH ) return false;
	}
	return true;
}
bool MSFilterSearch::checkCharge ( MSMSDataPoint& mmdp )
{
	bool allCharges = filterParams.getAllCharges ();
	if ( !allCharges ) {
		StringVector chargeFilter = filterParams.getChargeFilter ();
		SetInt si;
		bool tenAndAbove = false;
		for ( StringVectorSizeType i = 0 ; i < chargeFilter.size () ; i++ ) {
			if ( chargeFilter [i] == "10 and above" )
				tenAndAbove = true;
			else
				si.insert ( atoi ( chargeFilter [i].c_str () ) );
		}
		int z = mmdp.getPrecursorCharge ();
		if ( tenAndAbove && z >= 10 ) return true;
		return si.find ( z ) != si.end ();  
	}
	return true;
}
bool MSFilterSearch::checkNeutralLoss ( MSMSDataPoint& mmdp )
{
	string lossFormula = filterParams.getLossFormula ();			// Phospho "H3 P O4" Met ox "S O C H4"
	if ( !lossFormula.empty () ) {
		double precursorMZ = mmdp.getPrecursorMZ ();
		double precursorZ = mmdp.getPrecursorCharge ();
		int protonLoss = 0;
		if ( instInf->getChargeReducedFragmentation () ) protonLoss = 1;
		double protonMass = formula_to_monoisotopic_mass ( "H" ) - ELECTRON_REST_MASS;
		//double testMass = ( ( precursorMZ * precursorZ )/* - ( protonLoss * protonMass )*/ ) / ( precursorZ - protonLoss );
		double testMass = ( precursorMZ * precursorZ ) / ( precursorZ - protonLoss );
		double lossMass = formula_to_monoisotopic_mass ( lossFormula.c_str () );
		testMass -= ( lossMass / ( precursorZ - protonLoss ) );
		return checkNeutralLoss ( &mmdp, testMass );
	}
	return true;
}
bool MSFilterSearch::checkNeutralLoss ( MSMSDataPoint* mmdp, double testMass )
{
	PeakContainer::setMultiChargeAssign ( false );
	MSMSPeakFilterOptions* mmpfo = filterParams.getMSMSPeakFilterOptions ();
	Tolerance* parentTol = filterParams.getParentMassTolerance ();
	Tolerance* fragTol = filterParams.getProductMassTolerance ();
	bool monoisotopicFlag = filterParams.getMonoisotopicFlag ();
	bool averageParentMonoFragments = filterParams.getAverageParentMonoFragments ();
	Peak parentPeak ( mmdp->getPrecursorMZ (), mmdp->getPrecursorTolerance (), mmdp->getPrecursorCharge (), mmdp->getPrecursorIntensity (), getAdductMass ( monoisotopicFlag ), averageParentMonoFragments );
	PeakContainer peaks ( mmdp, mmpfo, &parentPeak, parentTol, fragTol, monoisotopicFlag, averageParentMonoFragments );
	for ( int i = peaks.size () ; i-- ; ) {
		if ( peaks [i]->isMatch ( testMass ) ) {
			return true;
		}
		if ( peaks [i]->isLowerMatch ( testMass ) ) break;
	}
	return false;
}
bool MSFilterSearch::checkFragmentMZ ( MSMSDataPoint& mmdp )
{
	double fragmentMZ = filterParams.getFragmentMZ ();
	if ( fragmentMZ ) {
		PeakContainer::setMultiChargeAssign ( false );
		MSMSPeakFilterOptions* mmpfo = filterParams.getMSMSPeakFilterOptions ();
		Tolerance* parentTol = filterParams.getParentMassTolerance ();
		Tolerance* fragTol = filterParams.getProductMassTolerance ();
		bool monoisotopicFlag = filterParams.getMonoisotopicFlag ();
		bool averageParentMonoFragments = filterParams.getAverageParentMonoFragments ();
		Peak parentPeak ( mmdp.getPrecursorMZ (), mmdp.getPrecursorTolerance (), mmdp.getPrecursorCharge (), mmdp.getPrecursorIntensity (), getAdductMass ( monoisotopicFlag ), averageParentMonoFragments );
		PeakContainer peaks ( &mmdp, mmpfo, &parentPeak, parentTol, fragTol, monoisotopicFlag, averageParentMonoFragments );
		for ( int i = 0 ; i < peaks.size () ; i++ ) {
			if ( peaks [i]->isMatch ( fragmentMZ ) ) {
				return true;
			}
			if ( peaks [i]->isUpperMatch ( fragmentMZ ) ) break;
		}
		return false;
	}
	return true;
}
bool MSFilterSearch::checkFragmentMZs ( MSMSDataPoint& mmdp )
{
	const DoubleVector& fragmentMZs = filterParams.getFragmentMZs ();
	int numFragments = fragmentMZs.size ();
	CharVector massMatched ( numFragments, 0 );
	if ( numFragments ) {
		PeakContainer::setMultiChargeAssign ( false );
		MSMSPeakFilterOptions* mmpfo = filterParams.getMSMSPeakFilterOptions ();
		Tolerance* parentTol = filterParams.getParentMassTolerance ();
		Tolerance* fragTol = filterParams.getProductMassTolerance ();
		bool monoisotopicFlag = filterParams.getMonoisotopicFlag ();
		bool averageParentMonoFragments = filterParams.getAverageParentMonoFragments ();
		Peak parentPeak ( mmdp.getPrecursorMZ (), mmdp.getPrecursorTolerance (), mmdp.getPrecursorCharge (), mmdp.getPrecursorIntensity (), getAdductMass ( monoisotopicFlag ), averageParentMonoFragments );
		PeakContainer peaks ( &mmdp, mmpfo, &parentPeak, parentTol, fragTol, monoisotopicFlag, averageParentMonoFragments );
		double lowMass = peaks.getMinMassMinusTol ();
		double highMass = peaks.getMaxMassPlusTol ();
		int numPeaks = peaks.size ();
		for ( int i = 0 ; i < numFragments ; i++ ) {
			double fragmentMass = fragmentMZs [i];
			int j = 0;
			if ( fragmentMass > highMass ) break;
			if ( fragmentMass > lowMass ) {
				for ( ; j < numPeaks && peaks [j]->isLowerMatch ( fragmentMass ) ; j++ ) {
					if ( !massMatched [i] ) {
						if ( peaks [j]->isMatch ( fragmentMass ) ) {
							massMatched [i] = 1;
						}
					}
				}
			}
		}
		return count ( massMatched.begin (), massMatched.end (), 1 ) >= filterParams.getMinMatches ();
	}
	return true;
}
