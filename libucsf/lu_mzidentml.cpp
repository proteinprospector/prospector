/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mzidentml.cpp                                              *
*                                                                             *
*  Created    : January 4th 2010                                              *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2010-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <algorithm>
#include <lg_string.h>
#include <lg_time.h>
#include <lu_getfil.h>
#include <lu_mzidentml.h>
#include <lu_xml.h>

using std::ostream;
using std::string;
using std::remove;

MZIdentML_AdditionalSearchParams::MZIdentML_AdditionalSearchParams ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "AdditionalSearchParams", subItems )
{
}

MZIdentML_Affiliation::MZIdentML_Affiliation ( const string& organizationRef ) :
	XMLOutputAttrItem ( "Affiliation" )
{
	attr.push_back ( makePairStringString ( "organization_ref", organizationRef ) );
}

MZIdentML_AmbiguousResidue::MZIdentML_AmbiguousResidue ( const VectorXMLOutputItemPtr& subItems, char code ) :
	XMLOutputContainerItem ( "AmbiguousResidue", subItems )
{
	attr.push_back ( makePairStringString ( "code", string ( 1, code ) ) );
}

MZIdentML_AnalysisCollection::MZIdentML_AnalysisCollection ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "AnalysisCollection", subItems )
{
}

MZIdentML_AnalysisData::MZIdentML_AnalysisData ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "AnalysisData", subItems )
{
}

MZIdentML_AnalysisParams::MZIdentML_AnalysisParams ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "AnalysisParams", subItems )
{
}

MZIdentML_AnalysisProtocolCollection::MZIdentML_AnalysisProtocolCollection ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "AnalysisProtocolCollection", subItems )
{
}

MZIdentML_AnalysisSampleCollection::MZIdentML_AnalysisSampleCollection ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "AnalysisSampleCollection", subItems )
{
}

MZIdentML_AnalysisSoftware::MZIdentML_AnalysisSoftware ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& name, const string& version ) :
	XMLOutputContainerItem ( "AnalysisSoftware", subItems, id )
{
	if ( !name.empty () )	attr.push_back ( makePairStringString ( "name", name ) );
	if ( !version.empty () )attr.push_back ( makePairStringString ( "version", version ) );
}

MZIdentML_AnalysisSoftwareList::MZIdentML_AnalysisSoftwareList ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "AnalysisSoftwareList", subItems )
{
}

MZIdentML_AuditCollection::MZIdentML_AuditCollection ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "AuditCollection", subItems )
{
}

MZIdentML_BibliographicReference::MZIdentML_BibliographicReference (
		const string& id,
		const string& authors,
		const string& title,
		const string& publication,
		const string& volume,
		const string& issue,
		const string& pages,
		const string& year,
		const string& doi,
		const string& editor,
		const string& publisher,
		const string& name ) :

	XMLOutputAttrItem ( "BibliographicReference", id )
{
	if ( !authors.empty () )	attr.push_back ( makePairStringString ( "authors", authors ) );
	if ( !title.empty () )		attr.push_back ( makePairStringString ( "title", title ) );
	if ( !publication.empty () )attr.push_back ( makePairStringString ( "publication", publication ) );
	if ( !volume.empty () )		attr.push_back ( makePairStringString ( "volume", volume ) );
	if ( !issue.empty () )		attr.push_back ( makePairStringString ( "issue", issue ) );
	if ( !pages.empty () )		attr.push_back ( makePairStringString ( "pages", pages ) );
	if ( !year.empty () )		attr.push_back ( makePairStringString ( "year", year ) );
	if ( !doi.empty () )		attr.push_back ( makePairStringString ( "doi", doi ) );
	if ( !editor.empty () )		attr.push_back ( makePairStringString ( "editor", editor ) );
	if ( !publisher.empty () )	attr.push_back ( makePairStringString ( "publisher", publisher ) );
	if ( !name.empty () )		attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_ContactRole::MZIdentML_ContactRole ( const VectorXMLOutputItemPtr& subItems, const string& contactRef ) :
	XMLOutputContainerItem ( "ContactRole", subItems )
{
	attr.push_back ( makePairStringString ( "Contact_ref", contactRef ) );
}

MZIdentML_Customizations::MZIdentML_Customizations ( const string& customizations ) :
	XMLOutputValueItem ( "Customizations", customizations )
{
}

MZIdentML_CV::MZIdentML_CV ( const string& id, const string& fullName, const string& URI, const string& version ) :
	XMLOutputAttrItem ( "cv", id )
{
	if ( !fullName.empty () )	attr.push_back ( makePairStringString ( "fullName", fullName ) );
	if ( !URI.empty () )		attr.push_back ( makePairStringString ( "URI", URI ) );
	if ( !version.empty () )	attr.push_back ( makePairStringString ( "version", version ) );
}

MZIdentML_CVList::MZIdentML_CVList ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "cvList", subItems )
{
}

MZIdentML_CVParam::MZIdentML_CVParam ( const string& accession, const string& name, const string& cvRef, const string& value, const string& unitAccession, const string& unitName, const string& unitCvRef ) :
	XMLOutputAttrItem ( "cvParam" )
{
	attr.push_back ( makePairStringString ( "accession", accession ) );
	attr.push_back ( makePairStringString ( "name", name ) );
	attr.push_back ( makePairStringString ( "cvRef", cvRef ) );
	if ( !value.empty () )			attr.push_back ( makePairStringString ( "value", value ) );
	if ( !unitAccession.empty () )	attr.push_back ( makePairStringString ( "unitAccession", unitAccession ) );
	if ( !unitName.empty () )		attr.push_back ( makePairStringString ( "unitName", unitName ) );
	if ( !unitCvRef.empty () )		attr.push_back ( makePairStringString ( "unitCvRef", unitCvRef ) );
}

MZIdentML_CVParam_AlternativeSingleLetterCodes::MZIdentML_CVParam_AlternativeSingleLetterCodes ( const string& value ) :
	MZIdentML_CVParam ( "MS:1001360", "alternate single letter codes", "PSI-MS", value )
{
}

MZIdentML_CVParam_ProteinDescription::MZIdentML_CVParam_ProteinDescription ( const string& value ) :
	MZIdentML_CVParam ( "MS:1001088", "protein description", "PSI-MS", value )
{
}

MZIdentML_CVParam_ProteinTaxonomy::MZIdentML_CVParam_ProteinTaxonomy ( const string& value ) :
	MZIdentML_CVParam ( "MS:1001469", "taxonomy: scientific name", "PSI-MS", value )
{
}

MZIdentML_CVParam_ProteinTaxonomyID::MZIdentML_CVParam_ProteinTaxonomyID ( const string& value ) :
	MZIdentML_CVParam ( "MS:1001467", "taxonomy: NCBI TaxID", "PSI-MS", value )
{
}

MZIdentML_CVParam_Tolerance::MZIdentML_CVParam_Tolerance ( double tol, const string& units, double sysError, const string& sign ) :
	MZIdentML_CVParam ( getAccession ( sign ), getName ( sign ), "PSI-MS", getValue ( tol, sysError, sign ), getUnitAccession ( units ), getUnitName ( units ), "UO" )
{
}
string MZIdentML_CVParam_Tolerance::getAccession ( const string& sign ) const
{
	if ( sign == "+" ) return "MS:1001412";
	if ( sign == "-" ) return "MS:1001413";
	return "";
}
string MZIdentML_CVParam_Tolerance::getName ( const string& sign ) const
{
	if ( sign == "+" ) return "search tolerance plus value";
	if ( sign == "-" ) return "search tolerance minus value";
	return "";
}
string MZIdentML_CVParam_Tolerance::getValue ( double tol, double sysError, const string& sign ) const
{
	double val;
	if ( sign == "+" ) val = tol + sysError;
	if ( sign == "-" ) val = tol - sysError;
	return gen_ftoa ( val, "%.4f" );
}
string MZIdentML_CVParam_Tolerance::getUnitAccession ( const string& units ) const
{
	if ( units == "Da" )	return "UO:0000221";
	if ( units == "%" )		return "UO:0000187";
	if ( units == "ppm" )	return "UO:0000169";
	return "";
}
string MZIdentML_CVParam_Tolerance::getUnitName ( const string& units ) const
{
	if ( units == "Da" )	return "dalton";
	if ( units == "%" )		return "percent";
	if ( units == "ppm" )	return "parts per million";
	return "";
}

MZIdentML_CVParam_Unimod::MZIdentML_CVParam_Unimod ( int accNum, const string& name ) :
	MZIdentML_CVParam ( "UNIMOD:" + gen_itoa ( accNum ), name, "UNIMOD" )
{
}

MZIdentML_DatabaseFilters::MZIdentML_DatabaseFilters ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "DatabaseFilters", subItems )
{
}

MZIdentML_DatabaseName::MZIdentML_DatabaseName ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "DatabaseName", subItems )
{
}

MZIdentML_DataCollection::MZIdentML_DataCollection ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "DataCollection", subItems )
{
}

MZIdentML_DBSequence::MZIdentML_DBSequence ( const VectorXMLOutputItemPtr& subItems, const string& searchDatabaseRef, const string& accession, int length, const string& name ) :
	XMLOutputContainerItem ( "DBSequence", subItems, "DBSeq_" + accession )
{
	if ( !searchDatabaseRef.empty () )	attr.push_back ( makePairStringString ( "SearchDatabase_ref", searchDatabaseRef ) );
	if ( !accession.empty () )			attr.push_back ( makePairStringString ( "accession", accession ) );
	if ( length != 0 )					attr.push_back ( makePairStringString ( "length", gen_itoa ( length ) ) );
	if ( !name.empty () )				attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_Enzyme::MZIdentML_Enzyme ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& cTermGain, const string& minDistance, const string& missedCleavages, const string& nTermGain, const string& name, const string& semiSpecific ) :
	XMLOutputContainerItem ( "Enzyme", subItems, id )
{
	if ( !cTermGain.empty () )			attr.push_back ( makePairStringString ( "cTermGain", cTermGain ) );
	if ( !minDistance.empty () )		attr.push_back ( makePairStringString ( "minDistance", minDistance ) );
	if ( !missedCleavages.empty () )	attr.push_back ( makePairStringString ( "missedCleavages", missedCleavages ) );
	if ( !nTermGain.empty () )			attr.push_back ( makePairStringString ( "nTermGain", nTermGain ) );
	if ( !name.empty () )				attr.push_back ( makePairStringString ( "name", name ) );
	if ( !semiSpecific.empty () )		attr.push_back ( makePairStringString ( "semiSpecific", semiSpecific ) );
}

MZIdentML_EnzymeName::MZIdentML_EnzymeName ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "EnzymeName", subItems )
{
}

MZIdentML_Enzymes::MZIdentML_Enzymes ( const VectorXMLOutputItemPtr& subItems, const string& independent ) :
	XMLOutputContainerItem ( "Enzymes", subItems )
{
	if ( !independent.empty () )	attr.push_back ( makePairStringString ( "independent", independent ) );
}

MZIdentML_FileFormat::MZIdentML_FileFormat ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "FileFormat", subItems )
{
}

MZIdentML_Filter::MZIdentML_Filter ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "Filter", subItems )
{
}

MZIdentML_FilterType::MZIdentML_FilterType ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "FilterType", subItems )
{
}

MZIdentML_FragmentArray::MZIdentML_FragmentArray ( const string& values, const string& measureRef ) :
	XMLOutputAttrItem ( "FragmentArray" )
{
	attr.push_back ( makePairStringString ( "values", values ) );
	attr.push_back ( makePairStringString ( "measure_ref", measureRef ) );
}

MZIdentML_Fragmentation::MZIdentML_Fragmentation ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "Fragmentation", subItems )
{
}

MZIdentML_FragmentationTable::MZIdentML_FragmentationTable ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "FragmentationTable", subItems )
{
}

MZIdentML_FragmentTolerance::MZIdentML_FragmentTolerance ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "FragmentTolerance", subItems )
{
}

MZIdentML_Include::MZIdentML_Include ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "Include", subItems )
{
}

MZIdentML_Inputs::MZIdentML_Inputs ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "Inputs", subItems )
{
}

MZIdentML_InputSpectra::MZIdentML_InputSpectra ( const string& spectraDataRef ) :
	XMLOutputAttrItem ( "InputSpectra" )
{
	if ( !spectraDataRef.empty () )	attr.push_back ( makePairStringString ( "spectraData_ref", spectraDataRef ) );
}

MZIdentML_InputSpectrumIdentifications::MZIdentML_InputSpectrumIdentifications ( const string& spectrumIdentificationListRef ) :
	XMLOutputAttrItem ( "InputSpectrumIdentifications" )
{
	if ( !spectrumIdentificationListRef.empty () )	attr.push_back ( makePairStringString ( "spectrumIdentificationList_ref", spectrumIdentificationListRef ) );
}

MZIdentML_IonType::MZIdentML_IonType ( const VectorXMLOutputItemPtr& subItems, const string& charge, const string& index ) :
	XMLOutputContainerItem ( "IonType", subItems )
{
	attr.push_back ( makePairStringString ( "charge", charge ) );
	if ( !index.empty () )	attr.push_back ( makePairStringString ( "index", index ) );
}

MZIdentML_MassTable::MZIdentML_MassTable ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& msLevel, const string& name ) :
	XMLOutputContainerItem ( "MassTable", subItems, id )
{
	attr.push_back ( makePairStringString ( "msLevel", msLevel ) );
	if ( !name.empty () )	attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_Measure::MZIdentML_Measure ( const VectorXMLOutputItemPtr& subItems, const string& id ) :
	XMLOutputContainerItem ( "Measure", subItems, id )
{
}

MZIdentML_Modification::MZIdentML_Modification ( const VectorXMLOutputItemPtr& subItems, int location, char residues, double monoisotopicMassDelta, double avgMassDelta ) :
	XMLOutputContainerItem ( "Modification", subItems )
{
	if ( location != -1 )				attr.push_back ( makePairStringString ( "location", gen_itoa ( location ) ) );
	if ( residues != 0 )				attr.push_back ( makePairStringString ( "residues", string ( 1, residues ) ) );
	if ( monoisotopicMassDelta != 0.0 )	attr.push_back ( makePairStringString ( "monoisotopicMassDelta", gen_ftoa ( monoisotopicMassDelta, "%.4f" ) ) );
	if ( avgMassDelta != 0.0 )			attr.push_back ( makePairStringString ( "avgMassDelta", gen_ftoa ( avgMassDelta, "%.4f" ) ) );
}

MZIdentML_ModificationParams::MZIdentML_ModificationParams ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "ModificationParams", subItems )
{
}

MZIdentML_MZIdentML::MZIdentML_MZIdentML () :
	XMLOutputContainerItem ( "mzIdentML" )
{
	attr.push_back ( makePairStringString ( "xmlns:xsi",			"http://www.w3.org/2001/XMLSchema-instance" ) );
	attr.push_back ( makePairStringString ( "xsi:schemaLocation",	"http://psidev.info/psi/pi/mzIdentML/1.1 ../schema/mzIdentML1.1.0.xsd" ) );
	attr.push_back ( makePairStringString ( "xmlns",				"http://psidev.info/psi/pi/mzIdentML/1.1" ) );
	attr.push_back ( makePairStringString ( "id",					"" ) );
	attr.push_back ( makePairStringString ( "creationDate",			genCurrentXSDDateAndTimeString () ) );
	attr.push_back ( makePairStringString ( "version",				"1.1.0" ) );
}

MZIdentML_Organization::MZIdentML_Organization ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& name ) :
	XMLOutputContainerItem ( "Organization", subItems, id )
{
	if ( !name.empty () ) attr.push_back ( makePairStringString ( "name", name ) );
}	

MZIdentML_Parent::MZIdentML_Parent ( const string& organizationRef ) :
	XMLOutputAttrItem ( "Parent" )
{
	attr.push_back ( makePairStringString ( "organization_ref", organizationRef ) );
}

MZIdentML_ParentTolerance::MZIdentML_ParentTolerance ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "ParentTolerance", subItems )
{
}

MZIdentML_Peptide::MZIdentML_Peptide ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& name ) :
	XMLOutputContainerItem ( "Peptide", subItems, id )
{
	if ( !name.empty () ) attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_PeptideEvidence::MZIdentML_PeptideEvidence ( const string& id, const string& dBSequenceRef, const string& peptideRef, int start, int end, const string& pre, const string& post, bool isDecoy, const string& frame, const string& name ) :
	XMLOutputAttrItem ( "PeptideEvidence", id )
{
	attr.push_back ( makePairStringString ( "dBSequence_ref", dBSequenceRef ) );
	attr.push_back ( makePairStringString ( "peptide_ref", peptideRef ) );
	if ( start != 0 )		attr.push_back ( makePairStringString ( "start", gen_itoa ( start ) ) );
	if ( end != 0 )			attr.push_back ( makePairStringString ( "end", gen_itoa ( end ) ) );
	if ( !pre.empty () )	attr.push_back ( makePairStringString ( "pre", pre ) );
	if ( !post.empty () )	attr.push_back ( makePairStringString ( "post", post ) );
	if ( isDecoy )			attr.push_back ( makePairStringString ( "isDecoy", "true" ) );
	if ( !frame.empty () )	attr.push_back ( makePairStringString ( "frame", frame ) );
	if ( !name.empty () )	attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_PeptideEvidenceRef::MZIdentML_PeptideEvidenceRef ( const string& peptideEvidenceRef ) :
	XMLOutputAttrItem ( "PeptideEvidenceRef" )
{
	attr.push_back ( makePairStringString ( "peptideEvidence_ref", peptideEvidenceRef ) );
}

MZIdentML_PeptideHypothesis::MZIdentML_PeptideHypothesis ( const VectorXMLOutputItemPtr& subItems, const string& peptideEvidenceRef ) :
	XMLOutputContainerItem ( "PeptideHypothesis", subItems )
{
	attr.push_back ( makePairStringString ( "peptideEvidence_ref", peptideEvidenceRef ) );
}

MZIdentML_PeptideSequence::MZIdentML_PeptideSequence ( const string& peptideSequence ) :
	XMLOutputValueItem ( "PeptideSequence", peptideSequence )
{
}

MZIdentML_Person::MZIdentML_Person ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& firstName, const string& midInitials, const string& lastName, const string& name ) :
	XMLOutputContainerItem ( "Person", subItems, id )
{
	if ( !firstName.empty () )	attr.push_back ( makePairStringString ( "firstName", firstName ) );
	if ( !midInitials.empty () )attr.push_back ( makePairStringString ( "midInitials", midInitials ) );
	if ( !lastName.empty () )	attr.push_back ( makePairStringString ( "lastName", lastName ) );
	if ( !name.empty () )		attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_ProteinAmbiguityGroup::MZIdentML_ProteinAmbiguityGroup ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& name ) :
	XMLOutputContainerItem ( "ProteinAmbiguityGroup", subItems, id )
{
	if ( !name.empty () )	attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_ProteinDetection::MZIdentML_ProteinDetection ( const VectorXMLOutputItemPtr& subItems,
		const string& id,
		const string& proteinDetectionListRef,
		const string& proteinDetectionProtocolRef,
		const string& name ) :

	XMLOutputContainerItem ( "ProteinDetection", subItems, id )
{
	attr.push_back ( makePairStringString ( "proteinDetectionList_ref", proteinDetectionListRef ) );
	attr.push_back ( makePairStringString ( "proteinDetectionProtocol_ref", proteinDetectionProtocolRef ) );
	if ( !name.empty () )	attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_ProteinDetectionHypothesis::MZIdentML_ProteinDetectionHypothesis ( const VectorXMLOutputItemPtr& subItems,
		const string& id,
		const string& passThreshold,
		const string& dBSequenceRef,
		const string& name ) :

	XMLOutputContainerItem ( "ProteinDetectionHypothesis", subItems, id )
{
	attr.push_back ( makePairStringString ( "passThreshold", passThreshold ) );
	if ( !dBSequenceRef.empty () )	attr.push_back ( makePairStringString ( "dBSequence_ref", dBSequenceRef ) );
	if ( !name.empty () )			attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_ProteinDetectionList::MZIdentML_ProteinDetectionList ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& name ) :
	XMLOutputContainerItem ( "ProteinDetectionList", subItems, id )
{
	if ( !name.empty () ) attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_ProteinDetectionProtocol::MZIdentML_ProteinDetectionProtocol ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& analysisSoftwareRef, const string& name ) :
	XMLOutputContainerItem ( "ProteinDetectionProtocol", subItems, id )
{
	attr.push_back ( makePairStringString ( "analysisSoftware_ref", analysisSoftwareRef ) );
	if ( !name.empty () ) attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_Provider::MZIdentML_Provider ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "Provider", subItems )
{
	attr.push_back ( makePairStringString ( "id", "PROVIDER" ) );
}

MZIdentML_SearchDatabase::MZIdentML_SearchDatabase ( const VectorXMLOutputItemPtr& subItems,
		const string& id,
		const string& location,
		const string& name,
		int numDatabaseSequences,
		const string& numResidues,
		const string& releaseDate,
		const string& version ) :

	XMLOutputContainerItem ( "SearchDatabase", subItems, id )
{
	attr.push_back ( makePairStringString ( "location", location ) );
	if ( !name.empty () )				attr.push_back ( makePairStringString ( "name", name ) );
	if ( numDatabaseSequences != -1 )	attr.push_back ( makePairStringString ( "numDatabaseSequences", gen_itoa ( numDatabaseSequences ) ) );
	if ( !numResidues.empty () )		attr.push_back ( makePairStringString ( "numResidues", numResidues ) );
	if ( !releaseDate.empty () )		attr.push_back ( makePairStringString ( "releaseDate", releaseDate ) );
	if ( !version.empty () )			attr.push_back ( makePairStringString ( "version", version ) );
}

MZIdentML_SiteRegexp::MZIdentML_SiteRegexp ( const string& regexp ) :
	XMLOutputValueItem ( "SiteRegexp", regexp )
{
}

MZIdentML_Residue::MZIdentML_Residue ( char code, double mass ) :
	XMLOutputAttrItem ( "Residue" )
{
	attr.push_back ( makePairStringString ( "code", string ( 1, code ) ) );
	attr.push_back ( makePairStringString ( "mass", gen_ftoa ( mass, "%.6f" ) ) );
}

MZIdentML_Role::MZIdentML_Role ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "Role", subItems )
{
}

MZIdentML_Sample::MZIdentML_Sample ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& name ) :
	XMLOutputContainerItem ( "Sample", subItems, id )
{
	if ( !name.empty () )	attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_SearchDatabaseRef::MZIdentML_SearchDatabaseRef ( const string& searchDatabaseRef ) :
	XMLOutputAttrItem ( "SearchDatabaseRef" )
{
	if ( !searchDatabaseRef.empty () )	attr.push_back ( makePairStringString ( "searchDatabase_ref", searchDatabaseRef ) );
}

MZIdentML_SearchModification::MZIdentML_SearchModification ( const VectorXMLOutputItemPtr& subItems,
		const string& fixedMod,
		const string& massDelta,
		const string& residues ) :

	XMLOutputContainerItem ( "SearchModification", subItems )
{
	attr.push_back ( makePairStringString ( "fixedMod", fixedMod ) );
	attr.push_back ( makePairStringString ( "massDelta", massDelta ) );
	attr.push_back ( makePairStringString ( "residues", residues ) );
}

MZIdentML_SearchType::MZIdentML_SearchType ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "SearchType", subItems )
{
}

MZIdentML_Seq::MZIdentML_Seq ( const string& seq ) :
	XMLOutputValueItem ( "Seq", seq )
{
}

MZIdentML_SequenceCollection::MZIdentML_SequenceCollection ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "SequenceCollection", subItems )
{
}

MZIdentML_SoftwareName::MZIdentML_SoftwareName ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "SoftwareName", subItems )
{
}

MZIdentML_SourceFile::MZIdentML_SourceFile ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& location, const string& name ) :
	XMLOutputContainerItem ( "SourceFile", subItems, id )
{
	attr.push_back ( makePairStringString ( "location", location ) );
	if ( !name.empty () ) attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_SpecificityRules::MZIdentML_SpecificityRules ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "SpecificityRules", subItems )
{
}

MZIdentML_SpectraData::MZIdentML_SpectraData ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& location, const string& name ) :
	XMLOutputContainerItem ( "SpectraData", subItems, id )
{
	attr.push_back ( makePairStringString ( "location", location ) );
	if ( !name.empty () ) attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_SpectrumIdentification::MZIdentML_SpectrumIdentification ( const VectorXMLOutputItemPtr& subItems,
		const string& id,
		const string& spectrumIdentificationListRef,
		const string& spectrumIdentificationProtocolRef,
		const string& activityDate,
		const string& name ) :

	XMLOutputContainerItem ( "SpectrumIdentification", subItems, id )
{
	attr.push_back ( makePairStringString ( "spectrumIdentificationList_ref", spectrumIdentificationListRef ) );
	attr.push_back ( makePairStringString ( "spectrumIdentificationProtocol_ref", spectrumIdentificationProtocolRef ) );
	if ( !activityDate.empty () )	attr.push_back ( makePairStringString ( "activityDate", activityDate ) );
	if ( !name.empty () )			attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_SpectrumIdentificationItem::MZIdentML_SpectrumIdentificationItem ( const VectorXMLOutputItemPtr& subItems,
		const string& id, double calculatedMassToCharge, int chargeState, double experimentalMassToCharge,
		bool passThreshold, int rank, const string& peptideRef, const string& massTableRef,
		const string& sampleRef, const string& calculatedPI, const string& name ) :

	XMLOutputContainerItem ( "SpectrumIdentificationItem", subItems, id )
{
	attr.push_back ( makePairStringString ( "calculatedMassToCharge", gen_ftoa ( calculatedMassToCharge, "%.4f" ) ) );
	attr.push_back ( makePairStringString ( "chargeState", gen_itoa ( chargeState ) ) );
	attr.push_back ( makePairStringString ( "experimentalMassToCharge", gen_ftoa ( experimentalMassToCharge, "%.4f" ) ) );
	attr.push_back ( makePairStringString ( "passThreshold", passThreshold ? "true" : "false" ) );
	attr.push_back ( makePairStringString ( "rank", gen_itoa ( rank ) ) );
	if ( !peptideRef.empty () )		attr.push_back ( makePairStringString ( "peptide_ref", peptideRef ) );
	if ( !massTableRef.empty () )	attr.push_back ( makePairStringString ( "massTable_ref", massTableRef ) );
	if ( !sampleRef.empty () )		attr.push_back ( makePairStringString ( "sample_ref", sampleRef ) );
	if ( !calculatedPI.empty () )	attr.push_back ( makePairStringString ( "calculatedPI", calculatedPI ) );
	if ( !name.empty () )			attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_SpectrumIdentificationItemRef::MZIdentML_SpectrumIdentificationItemRef (  const string& spectrumIdentificationItemRef ) :
	XMLOutputAttrItem ( "SpectrumIdentificationItemRef" )
{
	attr.push_back ( makePairStringString ( "spectrumIdentificationItem_ref", spectrumIdentificationItemRef ) );
}

MZIdentML_SpectrumIdentificationList::MZIdentML_SpectrumIdentificationList ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& numSequencesSearched, const string& name ) :
	XMLOutputContainerItem ( "SpectrumIdentificationList", subItems, id )
{
	if ( !numSequencesSearched.empty () )	attr.push_back ( makePairStringString ( "numSequencesSearched", numSequencesSearched ) );
	if ( !name.empty () )					attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_SpectrumIdentificationProtocol::MZIdentML_SpectrumIdentificationProtocol ( const VectorXMLOutputItemPtr& subItems, const string& id, const string& analysisSoftwareRef, const string& name ) :
	XMLOutputContainerItem ( "SpectrumIdentificationProtocol", subItems, id )
{
	attr.push_back ( makePairStringString ( "analysisSoftware_ref", analysisSoftwareRef ) );
	if ( !name.empty () )	attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_SpectrumIdentificationResult::MZIdentML_SpectrumIdentificationResult ( const VectorXMLOutputItemPtr& subItems,
		const string& id,
		const string& spectrumID,
		const string& spectraDataRef,
		const string& name ) :

	XMLOutputContainerItem ( "SpectrumIdentificationResult", subItems, id )
{
	attr.push_back ( makePairStringString ( "spectrumID", spectrumID ) );
	attr.push_back ( makePairStringString ( "spectraData_ref", spectraDataRef ) );
	if ( !name.empty () )	attr.push_back ( makePairStringString ( "name", name ) );
}

MZIdentML_SpectrumIDFormat::MZIdentML_SpectrumIDFormat ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "SpectrumIDFormat", subItems )
{
}

MZIdentML_SubSample::MZIdentML_SubSample ( const string& sampleRef ) :
	XMLOutputAttrItem ( "SubSample" )
{
	attr.push_back ( makePairStringString ( "sample_ref", sampleRef ) );
}

MZIdentML_Threshold::MZIdentML_Threshold ( const VectorXMLOutputItemPtr& subItems ) :
	XMLOutputContainerItem ( "Threshold", subItems )
{
}

MZIdentML_userParam::MZIdentML_userParam ( const string& name, const string& value ) :
	XMLOutputAttrItem ( "userParam" )
{
	attr.push_back ( makePairStringString ( "name", name ) );
	if ( !value.empty () ) attr.push_back ( makePairStringString ( "value", value ) );
}
