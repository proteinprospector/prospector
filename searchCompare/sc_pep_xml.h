/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_pep_xml.h                                                  *
*                                                                             *
*  Created    : December 21st 2010                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2010-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __sc_pep_xml_h
#define __sc_pep_xml_h

#include <lu_pros_form.h>
#include <lu_pep_xml.h>

class ParameterList;

class SCPepXMLReport {
	VectorConstParameterListPtr pList;
	PepXMLMSMSPipelineAnalysis pipeAnal;
	PepXMLAnalysisSummary analSumm;
	PepXMLMSMSRunSummary runSumm;
	PepXMLSampleEnzyme* pxse;
	std::vector <PepXMLSearchSummary*> sSumm;
	std::vector <PepXMLSearchDatabase*> pxsd;
	std::vector <PepXMLEnzymaticSearchConstraint*> pxesc;
	VectorXMLOutputItemPtr enzymeItems;
	std::vector <std::vector <PepXMLAminoAcidModification*> > vpxaam;
	std::vector <std::vector <PepXMLTerminalModification*> > vpxtm;
	std::string getPrecursorMassType ( const std::string& pmc );
	std::string getFragmentMassType ( const std::string& pmc );
	void setSearchSummary ();
	void setSampleEnzyme ();
	void setSearchDatabase ();
	void setEnzymaticSearchConstraint ();
	void setConstantMods ();
	void setVariableMods ();
public:
	SCPepXMLReport ();
	~SCPepXMLReport ();
	void updateFractionName ( const std::string& fractionName );
	void printHeader ( std::ostream& os ) const;
	void printFooter ( std::ostream& os ) const;
};

class SearchResultsPeptideLine;

class SCPepXMLReport2 {
	int num;
	int cur;
	std::vector <const SearchResultsPeptideLine*> vsrpl;
	std::vector <VectorXMLOutputItemPtr> subItems2;
	VectorXMLOutputItemPtr subItems2a;
	std::vector <VectorXMLOutputItemPtr> modInfoSubItems;
	VectorXMLOutputItemPtr pxsh;
	PepXMLSearchResult* pxsr;
	VectorXMLOutputItemPtr subItems4;
	PepXMLSpectrumQuery* pxsq;
	void createSearchHit ( const SearchResultsPeptideLine* s );
	void createSearchResult ( const SearchResultsPeptideLine* s );
	void createSpectrumQuery ( const SearchResultsPeptideLine* s, int num );
	void init ();
	void addAdditionalAccessionNumbers ( const SearchResultsPeptideLine* s );
public:
	SCPepXMLReport2 ( const SearchResultsPeptideLine* s, int num );
	~SCPepXMLReport2 ();
	void add ( const SearchResultsPeptideLine* s );
	void print ( std::ostream& os );
};

#endif /* ! __sc_pep_xml_h */
