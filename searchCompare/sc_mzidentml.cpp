/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_mzidentml.cpp                                              *
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
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#include <lu_aa_info.h>
#include <lu_getfil.h>
#include <lu_param_list.h>
#include <sc_sres_rep.h>
#include <sc_mzidentml.h>

using std::ostream;
using std::string;

SCMZIdentMLReport::SCMZIdentMLReport ()
{
	setCVList ();
	setAnalysisSoftwareList ();
	setProvider ();
	setAuditCollection ();
	setAnalysisSampleCollection ();
	setBibliographicReference ();
}
SCMZIdentMLReport::~SCMZIdentMLReport ()
{
	delete cvList;
	delete analysisSoftwareList;
	delete provider;
	delete auditCollection;
	delete analysisSampleCollection;
	delete bibliographicReference;
}
void SCMZIdentMLReport::setCVList ()
{
	VectorXMLOutputItemPtr subItems;
	subItems.push_back ( new MZIdentML_CV ( "PSI-MS","Proteomics Standards Initiative Mass Spectrometry Vocabularies", "http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo", "2.25.0" ) );
	subItems.push_back ( new MZIdentML_CV ( "UNIMOD","UNIMOD",		"http://www.unimod.org/obo/unimod.obo" ) );
	subItems.push_back ( new MZIdentML_CV ( "UO",	"UNIT-ONTOLOGY","http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo" ) );
	cvList = new MZIdentML_CVList ( subItems );
}
void SCMZIdentMLReport::setAnalysisSoftwareList ()
{
	VectorXMLOutputItemPtr subItems;
	subItems.push_back ( new MZIdentML_AnalysisSoftware ( VectorXMLOutputItemPtr (), "AS_ProteinProspector_Batch-Tag", "Protein Prospector Batch-Tag", PeptidePosition::getParams0 ()->getStringValue ( "version" ) ) );
	subItems.push_back ( new MZIdentML_AnalysisSoftware ( VectorXMLOutputItemPtr (), "AS_ProteinProspector_Search_Compare", "Protein Prospector Search Compare", Version::instance ().getVersion () ) );
	analysisSoftwareList = new MZIdentML_AnalysisSoftwareList ( subItems );
}
void SCMZIdentMLReport::setProvider ()
{
	VectorXMLOutputItemPtr subItems1;
	subItems1.push_back ( new MZIdentML_CVParam ( "MS:1001271", "researcher", "PSI-MS" ) );

	VectorXMLOutputItemPtr subItems2;
	subItems2.push_back ( new MZIdentML_Role ( subItems1 ) );

	VectorXMLOutputItemPtr subItems3;
	subItems3.push_back ( new MZIdentML_ContactRole ( subItems2, "PERSON_DOC_OWNER" ) );

	provider = new MZIdentML_Provider ( subItems3 );
}
void SCMZIdentMLReport::setAuditCollection ()
{
	VectorXMLOutputItemPtr subItems;
	VectorXMLOutputItemPtr subItems2;
	subItems.push_back ( new MZIdentML_Person ( subItems2, "per_pp" ) );
	VectorXMLOutputItemPtr subItems3;
	subItems3.push_back ( new MZIdentML_CVParam ( "MS:1000589", "contact email", "PSI-MS", "ppadmin@cgl.ucsf.edu" ) );
	subItems3.push_back ( new MZIdentML_CVParam ( "MS:1000587", "contact address", "PSI-MS", "Mass Spectrometry Facility, UCSF, 600 16th St, San Francisco, CA 94158-2517" ) );
	subItems3.push_back ( new MZIdentML_userParam ( "contact phone", "415 476 5189" ) );
	subItems3.push_back ( new MZIdentML_userParam ( "contact fax", "415 502 1655" ) );
	subItems.push_back ( new MZIdentML_Organization ( subItems3, "org_pp", "UCSF Mass Spectrometry Facility" ) );
	auditCollection = new MZIdentML_AuditCollection ( subItems );
}
void SCMZIdentMLReport::setAnalysisSampleCollection ()
{
	VectorXMLOutputItemPtr subItems;
	VectorXMLOutputItemPtr subItems2;
	VectorXMLOutputItemPtr subItems2A;
	VectorXMLOutputItemPtr subItems2A1;
	subItems2A1.push_back ( new MZIdentML_CVParam ( "MS:1001267", "software vendor", "PSI-MS" ) );
	subItems2A.push_back ( new MZIdentML_Role ( subItems2A1 ) );
	subItems2.push_back ( new MZIdentML_ContactRole ( subItems2A, "org_pp" ) );
	subItems.push_back ( new MZIdentML_Sample ( subItems2, "sample1" ) );
	analysisSampleCollection = new MZIdentML_AnalysisSampleCollection ( subItems );
}
void SCMZIdentMLReport::setBibliographicReference ()
{
	bibliographicReference = new MZIdentML_BibliographicReference (
		"bib1",
		"Chalkley R. J., Baker P. R., Huang L., Hansen K. C., Allen N. P., Rexach M. and Burlingame A. L.",
		"Comprehensive Analysis of a Multidimensional Liquid Chromatography Mass Spectrometry Dataset Acquired on a Quadrupole Selecting Quadrupole Collision Cell, Time-of-flight Mass Spectrometer. II. New Developments in Protein Prospector Allow for Reliable and Comprehensive Automatic Analysis of Large Datasets",
		"Molecular and Cellular Proteomics",
		"4",
		"8",
		"1194-1204",
		"2005",
		"10.1074/mcp.D500002-MCP200",
		"",
		"American Society for Biochemistry and Molecular Biology",
		"" );
}
void SCMZIdentMLReport::printHeader ( ostream& os ) const
{
	printXMLHeader ( os );
	mzIdentML.printOpenTag ( os, 0 );
	cvList->print ( os, 1 );
	analysisSoftwareList->print ( os, 1 );
	provider->print ( os, 1 );
	auditCollection->print ( os, 1 );
	analysisSampleCollection->print ( os, 1 );
}
void SCMZIdentMLReport::printFooter ( ostream& os ) const
{
	bibliographicReference->print ( os, 1 );
	mzIdentML.printCloseTag ( os, 0 );
}

SCMZIdentML_DBSequence::SCMZIdentML_DBSequence ( const SearchResultsProteinLine* s )
{
	VectorXMLOutputItemPtr subItems;
	//subItems.push_back ( new MZIdentML_Seq ( s->getProteinSequence () ) );
	subItems.push_back ( new MZIdentML_CVParam_ProteinDescription ( s->getName () ) );
	subItems.push_back ( new MZIdentML_CVParam_ProteinTaxonomy ( s->getSpecies () ) );
	subItems.push_back ( new MZIdentML_CVParam_ProteinTaxonomyID ( "" ) );
	dbSequence = new MZIdentML_DBSequence ( subItems, s->getDatabaseMZIdentMLRef (), s->getAcc (), s->getLength () );
}
SCMZIdentML_DBSequence::~SCMZIdentML_DBSequence ()
{
	delete dbSequence;
}
void SCMZIdentML_DBSequence::print ( ostream& os, int ntab ) const
{
	dbSequence->print ( os, ntab );
}

SCMZIdentML_Peptide::SCMZIdentML_Peptide ( const SearchResultsPeptideLine* s, const string& id )
{
	VectorXMLOutputItemPtr subItems;

	const string& pep = s->getDBPeptide ();
	subItems.push_back ( new MZIdentML_PeptideSequence ( pep ) );

	if ( s->getModNTermMass () != 0.0 ) {	// N-term
		VectorXMLOutputItemPtr subItems2;
		subItems2.push_back ( new MZIdentML_CVParam_Unimod ( 0, s->getModNTerm () ) );
		subItems.push_back ( new MZIdentML_Modification ( subItems2, 0, 0, s->getModNTermMass () ) );
	}

	VectorPairIntPairStringDouble vpid;
	s->getModMassesIndiciesAndString ( vpid );
	for ( int i = 0 ; i < vpid.size () ; i++ ) {
		VectorXMLOutputItemPtr subItems3;
		subItems3.push_back ( new MZIdentML_CVParam_Unimod ( 0, vpid [i].second.first ) );
		subItems.push_back ( new MZIdentML_Modification ( subItems3, vpid [i].first, pep[vpid [i].first-1], vpid [i].second.second ) );
	}

	if ( s->getModCTermMass () != 0.0 ) {	// C-term
		VectorXMLOutputItemPtr subItems4;
		subItems4.push_back ( new MZIdentML_CVParam_Unimod ( 0, s->getModCTerm () ) );
		subItems.push_back ( new MZIdentML_Modification ( subItems4, pep.length () + 1, 0, s->getModCTermMass () ) );
	}

	peptide = new MZIdentML_Peptide ( subItems, id );
}
SCMZIdentML_Peptide::~SCMZIdentML_Peptide ()
{
	delete peptide;
}
void SCMZIdentML_Peptide::print ( ostream& os, int ntab ) const
{
	peptide->print ( os, ntab );
}
