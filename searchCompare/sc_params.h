/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_params.h                                                   *
*                                                                             *
*  Created    : April 9th 2003                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __sc_params_h
#define __sc_params_h

#include <map>
#include <string>
#include <lg_string.h>
#include <lu_mass.h>
#include <lu_mass_elem.h>
#include <lu_quan_multi.h>

class ParameterList;

class SearchCompareParams {
	std::string saveFormat;
	std::string outputDirectory;
	std::string outputFilename;
	std::string reportType;
	std::string reportHomologousProteins;
	StringVector databaseType;
	bool multiSample;
	StringVector idFilterList;
	SetInt idFilterSet;
	std::map <std::string, double> minBestDiscScore;
	bool bestDiscOnly;
	bool discScoreGraph;
	double modificationScoreThreshold;
	double maxProteinEValue;
	double maxPeptideEValue;
	double minPPProteinScore;
	double minPPPeptideScore;
	double xlMinLowScore;
	double xlMinScoreDiff;
	double xlMaxLowExpectation;
	bool remove;
	std::string peptideFilter;
	std::string reportHitsType;
	std::string mergeOption;
	std::string sortType;
	std::string sortType2;
	bool unmatchedSpectra;
	StringVector filenames;
	double resolution;
	bool reportMMods;
	bool reportCoverage;
	AAInitInfo aaii;
	IsotopePurity isotopePurity;
	PurityCorrection purityCorrection;
	double reporterIonWindow;
	bool checkboxes;
	std::vector <std::map <std::string, VectorPairStringString> > peptidesToDelete;

	MapStringToStringVector accessionNumbers;

	void setReportItems ( const ParameterList* params );
	void initPeptidesToDelete ( const StringVector& s );

	StringVector projectNames;
	StringVector resultsNames;
	StringVector resultsFullPaths;
public:
	SearchCompareParams ( const ParameterList* params );
	char getSaveFormat () const
	{
		if ( saveFormat == "HTML" ) return 'H';
		else if ( saveFormat == "mzIdentML" ) return 'M';
		else if ( isPrefix ( saveFormat, "pepXML" ) ) return 'P';
		else if ( saveFormat == "Filtered Peak Lists" ) return 'F';
		else if ( saveFormat == "MS-Viewer Files" ) return 'V';
		else if ( isPrefix ( saveFormat, "BiblioSpec" ) ) return 'B';
		else return 'T';
	}
	bool getNormalizedBiblioSpec () const
	{
		return saveFormat == "BiblioSpec (Normalized RT)";
	}
	std::string getOutputDirectory () const { return outputDirectory; }
	std::string getOutputFilename () const { return outputFilename; }
	std::string getReportType () const { return reportType; }
	std::string getReportHomologousProteins () const { return reportHomologousProteins; }
	bool getRemove () const { return remove; }
	double getMinBestDiscScore ( const std::string& instrument ) const
	{
		MapStringToDoubleConstIterator cur = minBestDiscScore.find ( instrument );
		if ( cur != minBestDiscScore.end () ) {
			return (*cur).second;
		}
		else return 0.0;
	}
	double getMaxPeptideEValue () const { return maxPeptideEValue; }
	bool getBestDiscOnly () const { return bestDiscOnly; }
	bool getDiscScoreGraph () const { return discScoreGraph && idFilterSet.empty (); }
	double getModificationScoreThreshold () const { return modificationScoreThreshold; }
	double getMaxProteinEValue () const { return maxProteinEValue; }
	double getMinPPProteinScore () const { return minPPProteinScore; }
	double getMinPPPeptideScore () const { return minPPPeptideScore; }
	double getXLMinLowScore () const { return xlMinLowScore; }
	double getXLMinScoreDiff () const { return xlMinScoreDiff; }
	double getXLMaxLowExpectation () const { return xlMaxLowExpectation; }
	std::string getPeptideFilter () const { return peptideFilter; }
	std::string getReportHitsType () const { return reportHitsType; }
	std::string getMergeOption () const { return mergeOption; }
	std::string getSortType () const { return sortType; }
	std::string getSortType2 () const { return sortType2; }
	bool getUnmatchedSpectra () const { return unmatchedSpectra; }
	StringVector getFilenames () const { return filenames; }
	StringVector getIDFilterList () const { return idFilterList; }
	SetInt getIDFilterSet () const { return idFilterSet; }
	bool getMultiSample () const { return multiSample; }
	bool getReportMMods () const { return reportMMods; }
	bool getReportCoverage () const { return reportCoverage; }
	bool getCheckboxes () const { return checkboxes; }
	MapStringToStringVector getAccessionNumbers () const { return accessionNumbers; }
	bool getRemovePeptides () const { return !peptidesToDelete.empty (); }
	VectorPairStringString getRemovePeptides ( int fileNum, const std::string& id ) const;
	MapStringConstModPtr getConstMods () const { return aaii.getConstMods (); }

	int getNumCommandLineSearches () const { return projectNames.size (); }
	StringVector getCommandLineProjectNames () const { return projectNames; }
	StringVector getCommandLineResultsNames () const { return resultsNames; }
	StringVector getCommandLineResultsFullPaths () const { return resultsFullPaths; }
};

#endif /* ! __sc_params_h */
