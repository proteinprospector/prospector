/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_viewer_par.h                                               *
*                                                                             *
*  Created    : March 8th 2011                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2011-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_viewer_par_h
#define __lu_viewer_par_h

#include <lgen_file.h>
#include <lu_prog_par.h>
#include <lu_tol.h>
#include <lu_mass.h>

class MSViewerFilter {
	IntVector filterColumn;
	StringVector filterType;
	StringVectorVector filterValue;
	static int getColumnNumber ( const std::string& name );
public:
	MSViewerFilter ();
	void initParameters ( const ParameterList* params );
	IntVector getFilterColumn () const { return filterColumn; }
	StringVector getFilterType () const { return filterType; }
	StringVectorVector getFilterValue () const { return filterValue; }
	void printHTML ( std::ostream& os ) const;
};
class MSViewerSortOrder {
	IntVector sortLevel;
	StringVector sortOrderType;
	StringVector sortOrderDirection;
	static int getColumnNumber ( const std::string& name );
public:
	MSViewerSortOrder ();
	void initParameters ( const ParameterList* params );
	IntVector getSortLevel () const { return sortLevel; }
	StringVector getSortOrderType () const { return sortOrderType; }
	StringVector getSortOrderDirection () const { return sortOrderDirection; }
	void printHTML ( std::ostream& os ) const;
};
class MSViewerRemoveReplicate {
	IntVector replicateTest;
	static int getColumnNumber ( const std::string& name );
public:
	MSViewerRemoveReplicate ();
	void initParameters ( const ParameterList* params );
	IntVector getReplicateTest () const { return replicateTest; }
	void printHTML ( std::ostream& os ) const;
};

class MSViewerParameters : public MSProgramParameters {
	const ParameterList* pList;
	bool pListFlag;

	bool saveSettings;
	std::string searchKey;
	bool deleteFlag;
	std::string viewerOutputType;
	std::string resultsFileFormat;
	std::string instrumentFilter;
	double probabilityLimit;

	std::string peakListFpath;
	std::string resultsFpath;
	std::string tempFileFullPath;
	std::string resName;

	ToleranceInfo* parentTolerance;
	ToleranceInfo* fragmentTolerance;
	std::string instrumentName;
	std::string linkSearchType;
	std::string bridgeFormula;

	int fractionColumnNumber;
	std::string spectrumIdentifier;
	int scanIDColumnNumber;

	int peptideColumnNumber;
	int zColumnNumber;
	std::string modificationReporting;
	int constantModColumnNumber;
	int variableModColumnNumber;
	int allModColumnNumber;

	std::string separator;
	int numTitleLines;
	int numHeaderLines;
	int page;
	std::string rowsPerPage;
	AAInitInfo* aaInitInfo;
	StringVector removeColumn;
	MSViewerSortOrder msvSortOrder;
	MSViewerFilter msvFilter;
	MSViewerRemoveReplicate msvRemoveReplicate;
	bool commandLine;
	int numPepXML;
	int numMSF;
	static int getColumnNumber ( const std::string& name );
	static std::string getColumnSeparator ( const std::string& name );
	void init ( const ParameterList* params );
	void scriptConversion ();
	void initPeakList ( const ParameterList* params );
	void initResults ( const ParameterList* params );
	void initParams ( const ParameterList* params );
public:
	MSViewerParameters ( const ParameterList* params );
	~MSViewerParameters ();

	std::string getParamsSearchKey ( const ParameterList* params ) const;
	bool getSaveSettings () const { return saveSettings; }
	std::string getSearchKey () const { return searchKey; }
	bool getDeleteFlag ()	const { return deleteFlag; }
	bool getTabDelimitedText ()	const { return viewerOutputType == "Tab delimited text"; }
	bool getTabDelimitedTextWithURL ()	const { return viewerOutputType == "Tab delimited text with URL"; }
	bool getFilteredViewerFiles ()	const { return viewerOutputType == "Filtered Viewer files"; }
	bool getViewerFiles ()			const { return viewerOutputType == "Viewer files"; }
	bool getProspector ()			const { return resultsFileFormat == "Protein Prospector Tab Delimited"; }
	bool getProspectorXL ()			const { return resultsFileFormat == "Protein Prospector Crosslinked Peptides Tab Delimited"; }
	bool getMaxQuant ()				const { return resultsFileFormat == "MaxQuant"; }
	bool getOther ()				const { return resultsFileFormat == "Other"; }
	bool getAuto () const
	{
		return resultsFileFormat == "PRIDE XML" || resultsFileFormat == "pepXML" || resultsFileFormat == "BiblioSpec" || resultsFileFormat == "Thermo MSF" || resultsFileFormat == "NIST MSP" || resultsFileFormat == "SpectraST sptxt";
	}
	std::string getScript () const;
	void getScriptParameters ( int& numTitleLines, int& numHeaderLines, std::string& delimiter ) const;
	void getScriptParameters2 ( std::string& spectrumIdentifier, std::string& fraction, std::string& scanID, std::string& peptide, std::string& charge, std::string& modifications ) const;
	std::string getPeakListFpath () const
	{
		if ( genIsFullPath ( peakListFpath ) )	return peakListFpath;
		else {
			if ( peakListFpath.empty () ) return "";
			else return getViewerRepositoryPath ( searchKey ) + SLASH + peakListFpath;
		}
	}
	std::string getResultsFpath () const
	{
		if ( genIsFullPath ( resultsFpath ) )	return resultsFpath;
		else {
			if ( resultsFpath.empty () )	return "";
			else return getViewerRepositoryPath ( searchKey ) + SLASH + resultsFpath;
		}
	}
	std::string getResName () const { return resName; }
	std::string getKey () const;

	std::string getSeparator () const { return separator; }
	bool getCSV () const { return separator == ","; }
	int getNumTitleLines () const { return numTitleLines; }
	int getNumHeaderLines () const { return numHeaderLines; }
	int getPage () const { return page; }
	std::string getRowsPerPage () const { return rowsPerPage; }

	int getFractionColumnNumber () const { return fractionColumnNumber; }
	std::string getSpectrumIdentifier () const { return spectrumIdentifier; }
	int getScanIDColumnNumber () const { return scanIDColumnNumber; }

	int getPeptideColumnNumber () const { return peptideColumnNumber; }
	int getZColumnNumber () const { return zColumnNumber; }
	std::string getModificationReporting () const { return modificationReporting; }
	int getConstantModColumnNumber () const { return constantModColumnNumber; }
	int getVariableModColumnNumber () const { return variableModColumnNumber; }
	int getAllModColumnNumber () const { return allModColumnNumber; }

	BoolDeque getRemoveColumn ( int numCols ) const;
	IntVector getSortLevel () const { return msvSortOrder.getSortLevel (); }
	StringVector getSortOrderType () const { return msvSortOrder.getSortOrderType (); }
	StringVector getSortOrderDirection () const { return msvSortOrder.getSortOrderDirection (); }

	IntVector getFilterColumn () const { return msvFilter.getFilterColumn (); }
	StringVector getFilterType () const { return msvFilter.getFilterType (); }
	StringVectorVector getFilterValue () const { return msvFilter.getFilterValue (); }

	IntVector getReplicateTest () const { return msvRemoveReplicate.getReplicateTest (); }

	const Tolerance* getParentTolerance () const { return parentTolerance->getTolerance (); }
	const Tolerance* getFragmentTolerance () const { return fragmentTolerance->getTolerance (); }
	std::string getInstrumentName () const { return instrumentName; }
	std::string getLinkSearchType () const { return linkSearchType; }
	std::string getBridgeFormula () const { return bridgeFormula; }
	MapCharToString getConstMods () const { return aaInitInfo->getConstModMap (); }

	bool getPListFlag () const { return pListFlag; }
	const ParameterList* getParameterList () const { return pList; }
	static std::string getViewerRepositoryContainerPath ();
	static std::string getViewerRepositoryContainerPath ( const std::string& searchKey );
	static std::string getViewerRepositoryPath ( const std::string& searchKey );
	void deleteDataSet () const;
	void saveDataSet () const;

	void setPeakListFpath ( const std::string& f );
	void setResultsFpath ( const std::string& f );
	void setResultsFileFormat ( const std::string& f );
	void setSpectrumIdentifier ( const std::string& f );

	void setScanIDColumnNumber ( int n );
	void setPeptideColumnNumber ( int n );
	void setVariableModColumnNumber ( int n );
	void setZColumnNumber ( int n );
	void setFractionColumnNumber ( int n );

	void setSeparator ( const std::string& s );
	void setNumHeaderLines ( int n );
	void setPage ( int n );
	void setConstMod ( const StringVector& s );
	bool getCommandLine () const { return commandLine; }
	int getNumPepXML () const { return numPepXML; }
	int getNumMSF () const { return numMSF; }
	void processAPLFiles ();
};

#endif /* ! __lu_viewer_par_h */
