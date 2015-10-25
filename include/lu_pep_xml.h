/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pep_xml.h                                                  *
*                                                                             *
*  Created    :                                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2010-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_pep_xml_h
#define __lu_pep_xml_h

#include <lu_xml.h>

void printPepXMLHeader ( std::ostream& os );

class PepXMLMSMSPipelineAnalysis : public XMLOutputContainerItem {
	void init ( const std::string& name, const std::string& summaryXML );
public:
	PepXMLMSMSPipelineAnalysis ( const std::string& name, const std::string& summaryXML );
	PepXMLMSMSPipelineAnalysis ( const VectorXMLOutputItemPtr& subItems, const std::string& name, const std::string& summaryXML );
};

class PepXMLAnalysisSummary : public XMLOutputContainerItem {
	void init ( const std::string& analysis, const std::string& t );
public:
	PepXMLAnalysisSummary ( const std::string& analysis, const std::string& t );
	PepXMLAnalysisSummary ( const VectorXMLOutputItemPtr& subItems, const std::string& analysis, const std::string& t );
};

class PepXMLMSMSRunSummary : public XMLOutputContainerItem {
	void init ( const std::string& baseName, const std::string& rawDataType, const std::string& rawData, const std::string& msManufacturer, const std::string& msModel, const std::string& msIonization, const std::string& msMassAnalyzer, const std::string& msDetector );
public:
	PepXMLMSMSRunSummary ( const std::string& baseName, const std::string& rawDataType, const std::string& rawData, const std::string& msManufacturer = "", const std::string& msModel = "", const std::string& msIonization = "", const std::string& msMassAnalyzer = "", const std::string& msDetector = "" );
	PepXMLMSMSRunSummary ( const VectorXMLOutputItemPtr& subItems, const std::string& baseName, const std::string& rawDataType, const std::string& rawData, const std::string& msManufacturer = "", const std::string& msModel = "", const std::string& msIonization = "", const std::string& msMassAnalyzer = "", const std::string& msDetector = "" );
	void setBaseName ( const std::string& baseName );
};

class PepXMLSampleEnzyme : public XMLOutputContainerItem {
public:
	PepXMLSampleEnzyme ( const VectorXMLOutputItemPtr& subItems, const std::string& name, const std::string& description, const std::string& fidelity, int independent = -1 );
};

class PepXMLSpecificity : public XMLOutputAttrItem {
public:
	PepXMLSpecificity ( const std::string& sense, int minSpacing, const std::string& cut, const std::string& noCut );
};

class PepXMLSearchSummary : public XMLOutputContainerItem {
	void init ( const std::string& baseName, const std::string& searchEngine, const std::string& precursorMassType, const std::string& fragmentMassType, int searchID, const std::string& outDataType, const std::string& outData );
public:
	PepXMLSearchSummary ( const std::string& baseName, const std::string& searchEngine, const std::string& precursorMassType, const std::string& fragmentMassType, int searchID, const std::string& outDataType = "", const std::string& outData = "" );
	PepXMLSearchSummary ( const VectorXMLOutputItemPtr& subItems, const std::string& baseName, const std::string& searchEngine, const std::string& precursorMassType, const std::string& fragmentMassType, int searchID, const std::string& outDataType = "", const std::string& outData = "" );
	void setBaseName ( const std::string& baseName );
};

class PepXMLSearchDatabase : public XMLOutputAttrItem {
public:
	PepXMLSearchDatabase ( const std::string& localPath, const std::string& type, int sizeInDBEntries, const std::string& url = "", const std::string& databaseName = "", const std::string& origDatabaseURL = "", const std::string& databaseReleaseDate = "", const std::string& databaseReleaseIdentifier = "", int sizeOfResidues = -1 );
};

class PepXMLEnzymaticSearchConstraint : public XMLOutputAttrItem {
public:
	PepXMLEnzymaticSearchConstraint ( const std::string& enzyme, int maxNumInternalCleavages, int minNumberTermini );
};

class PepXMLAminoAcidModification : public XMLOutputAttrItem {
public:
	PepXMLAminoAcidModification ( const std::string& aminoacid, double massdiff, const std::string& variable, const std::string& peptideTerminus, const std::string& description, const std::string& binary = "", const std::string& symbol = "" );
};

class PepXMLTerminalModification : public XMLOutputAttrItem {
public:
	PepXMLTerminalModification ( const std::string& terminus, double massdiff, const std::string& variable, const std::string& proteinTerminus, const std::string& description, const std::string& symbol = "" );
};

class PepXMLParameter : public XMLOutputAttrItem {
public:
	PepXMLParameter ( const std::string& name, const std::string& value );
};

class PepXMLAnalysisTimestamp : public XMLOutputContainerItem {
public:
	PepXMLAnalysisTimestamp ( const VectorXMLOutputItemPtr& subItems, const std::string& analysis, const std::string& time, const std::string& id );
};

class PepXMLSpectrumQuery : public XMLOutputContainerItem {
public:
	PepXMLSpectrumQuery ( const VectorXMLOutputItemPtr& subItems, const std::string& spectrum,
		int startScan, int endScan, double precursorNeutralMass, int assumedCharge,
		const std::string& searchSpecification, int index, double retentionTimeSec,
		const std::string& activationMethod );
};

class PepXMLSearchResult : public XMLOutputContainerItem {
public:
	PepXMLSearchResult ( const VectorXMLOutputItemPtr& subItems, int searchID = -1 );
};

class PepXMLSearchHit : public XMLOutputContainerItem {
public:
	PepXMLSearchHit ( const VectorXMLOutputItemPtr& subItems, int hitRank, const std::string& peptide,
		const std::string& peptidePrevAA, const std::string& peptideNextAA, const std::string& protein,
		int numTotProteins, int numMatchedIons, int totNumIons, double calcNeutralPepMass, double massdiff,
		int numTolTerm,	int numMissedCleavages, int isRejected, const std::string& proteinDescr,
		double calcPI, double proteinMW );
};

class PepXMLAlternativeProtein : public XMLOutputAttrItem {
public:
	PepXMLAlternativeProtein ( const std::string& protein, const std::string& proteinDescr, int numTolTerm,
		double proteinMW, const std::string& peptidePrevAA, const std::string& peptideNextAA );
};

class PepXMLModificationInfo : public XMLOutputContainerItem {
public:
	PepXMLModificationInfo ( const VectorXMLOutputItemPtr& subItems, double modNTermMass, double modCTermMass, const std::string& modifiedPeptide = "" );
};

class PepXMLModAminoacidMass : public XMLOutputAttrItem {
public:
	PepXMLModAminoacidMass ( int position, double mass );
};

class PepXMLSearchScore : public XMLOutputAttrItem {
public:
	PepXMLSearchScore ( const std::string& name, double value );
};

#endif /* ! __lu_pep_xml_h */
