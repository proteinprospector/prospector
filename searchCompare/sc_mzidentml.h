/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_mzidentml.h                                                *
*                                                                             *
*  Created    : January 29th 2013                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2013-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __sc_mzidentml_h
#define __sc_mzidentml_h

#include <lu_pros_form.h>
#include <lu_mzidentml.h>

/*
Overview
<MzIdentML>
	<cvList>
	<AnalysisSoftwareList>
	<Provider>
	<AuditCollection>
	<AnalysisSampleCollection>
	<SequenceCollection>
	<AnalysisCollection>
	<AnalysisProtocolCollection>
	<DataCollection>
	<BibliographicReference>
</MzIdentML>
*/

class SCMZIdentMLReport {
	MZIdentML_MZIdentML mzIdentML;
	MZIdentML_CVList* cvList;
	MZIdentML_AnalysisSoftwareList* analysisSoftwareList;
	MZIdentML_Provider* provider;
	MZIdentML_AuditCollection* auditCollection;
	MZIdentML_AnalysisSampleCollection* analysisSampleCollection;

	MZIdentML_BibliographicReference* bibliographicReference;

	void setCVList ();
	void setAnalysisSoftwareList ();
	void setProvider ();
	void setAuditCollection ();
	void setAnalysisSampleCollection ();

	void setBibliographicReference ();
public:
	SCMZIdentMLReport ();
	~SCMZIdentMLReport ();
	void printHeader ( std::ostream& os ) const;
	void printFooter ( std::ostream& os ) const;
};

class SearchResultsProteinLine;
class SearchResultsPeptideLine;

class SCMZIdentML_DBSequence {
	MZIdentML_DBSequence* dbSequence;
public:
	SCMZIdentML_DBSequence ( const SearchResultsProteinLine* s );
	~SCMZIdentML_DBSequence ();
	void print ( std::ostream& os, int ntab ) const;
};

class SCMZIdentML_Peptide {
	MZIdentML_Peptide* peptide;
public:
	SCMZIdentML_Peptide ( const SearchResultsPeptideLine* s, const std::string& id );
	~SCMZIdentML_Peptide ();
	void print ( std::ostream& os, int ntab ) const;
};

#endif /* ! __sc_mzidentml_h */
