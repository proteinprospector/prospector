/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_viewer_srch.h                                              *
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

#ifndef __lu_viewer_srch_h
#define __lu_viewer_srch_h

#include <lu_viewer_par.h>
#include <lu_program.h>


class MSViewerSearch : public MSProgram {
	const MSViewerParameters& viewerParams;
	int spectrumIdentifier;
	std::string modificationReporting;

	StringVector fractionName;
	StringVector scanID;
	StringVector peptide;
	StringVector peptide2;
	StringVector linker;
	IntVector charge;
	StringVector constantMods;
	StringVector variableMods;
	StringVector allMods;

	std::string resultsFPath;

	std::string peakListFPath;
	StringVector peakListFractionNames;
	StringVector peakListCentroidFiles;
	StringVector titleLines;
	StringVectorVector headerLines;
	StringVectorVector rows;
	int extraTitleLines;
	int numCols;

	int fractionColumnNumber;
	int scanIDColumnNumber;
	int ppSpectrumColumnNumber;
	int peptideColumnNumber;
	int peptide2ColumnNumber;
	int linkerColumnNumber;
	int zColumnNumber;
	int constantModsColumnNumber;
	int variableModsColumnNumber;
	int allModsColumnNumber;

	MapCharToString constMods;
	void makeHeaderLines ( const StringVector& constantHeader, const SetString& variableHeader );
	void makeRowLines ( const SetString& variableHeader, const VectorMapStringToString rowsVariable );
	void writeResultsFile ( const std::string& resFile ) const;
	void writeMGFFile ( const std::string& dir, const std::string& shortFName, MSMSDataPointVector& msmsDataPointList, const IntVector& scans ) const;
	void processPeakListFile ();
	void processResultsFile ();
	void processPepXMLResultsFile ( bool multi );
	void processPrideXMLResultsFile ();
	void setParams ( bool peakList, const std::string& separator, const std::string& suffix );
	void processBLibLibrary ();
	void processLibrary ( const std::string& type );
	void processMSF ( bool multi );
	std::string getConstModsString ( const std::string& dbPeptide );
	std::string convertTandemMod ( const std::string& tMod, int startAA );
	void getMenuColumnNumbers ();
	int getColumnNumber ( const std::string& id ) const;
	int getOptionalColumnNumber ( const std::string& id ) const;
	void getColumnNumbers ();
	void getProspectorColumnNumbers ();
	void getProspectorXLColumnNumbers ();
	void setColumnNumbers ();
	void sortRows ();
	void filterRows ();
	void removeReplicates ();
	void setSpecialColumns ();
	void setFractionName ();
	void setScanID ();
	void setPeptide ();
	void setPeptide2 ();
	void setLinker ();
	void setCharge ();
	void setConstantMods ();
	void setVariableMods ();
	void setAllMods ();
	void processPeakListAndResults ();
	void urlDelimitedCell ( std::ostream& os, std::map <std::string, MSProductLink*>& productLink, int m );
	std::string getSpecIDStringFromScanID ( const std::string& scanID );
public:
	MSViewerSearch ( const MSViewerParameters& params );
	~MSViewerSearch ();
	void printBodyHTML ( std::ostream& os );
	void printSaveSettingsHTML ( std::ostream& os );
	void printMSViewerRepositoryHMTL ( std::ostream& os );
	void printMSViewerReportHMTL ( std::ostream& os );
	void printMSViewerReportTabDelimited ( std::ostream& os, bool url, bool printNumbers );
	void printMSViewerReportTabDelimitedText ( std::ostream& os, bool printNumbers );
	void printMSViewerReportTabDelimitedTextWithURL ( std::ostream& os, bool printNumbers );
	void printMSViewerFilteredPeakList ( std::ostream& os );
	void printMSViewerReportViewerFiles ( std::ostream& os );
};

#endif /* ! __lu_viewer_srch_h */
