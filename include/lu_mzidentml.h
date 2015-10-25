/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mzidentml.h                                                *
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
*  Copyright (2010-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_mzidentml_h
#define __lu_mzidentml_h

#include <ostream>
#include <lu_xml.h>

// 1.  AdditionalSearchParams             XMLOutputContainerItem
/*
<AdditionalSearchParams>
*/
class MZIdentML_AdditionalSearchParams : public XMLOutputContainerItem {
public:
	MZIdentML_AdditionalSearchParams ( const VectorXMLOutputItemPtr& subItems );
};

// 2.  Affiliation                        XMLOutputAttrItem
/*
<Affiliation organization_ref="ORG_MSL"/>
*/
class MZIdentML_Affiliation : public XMLOutputAttrItem {
public:
	MZIdentML_Affiliation ( const std::string& organizationRef );
};

// 3.  AmbiguousResidue                   XMLOutputContainerItem
/*
<AmbiguousResidue code="B">
*/
class MZIdentML_AmbiguousResidue : public XMLOutputContainerItem {
public:
	MZIdentML_AmbiguousResidue ( const VectorXMLOutputItemPtr& subItems, char code );
};

// 4.  AnalysisCollection                 XMLOutputContainerItem
/*
<AnalysisCollection>
*/
class MZIdentML_AnalysisCollection : public XMLOutputContainerItem {
public:
	MZIdentML_AnalysisCollection ( const VectorXMLOutputItemPtr& subItems = VectorXMLOutputItemPtr () );
};

// 5.  AnalysisData                       XMLOutputContainerItem
/*
<AnalysisData>
*/
class MZIdentML_AnalysisData : public XMLOutputContainerItem {
public:
	MZIdentML_AnalysisData ( const VectorXMLOutputItemPtr& subItems );
};

// 6.  AnalysisParams                     XMLOutputContainerItem
/*
<AnalysisParams>
*/
class MZIdentML_AnalysisParams : public XMLOutputContainerItem {
public:
	MZIdentML_AnalysisParams ( const VectorXMLOutputItemPtr& subItems );
};

// 7.  AnalysisProtocolCollection         XMLOutputContainerItem
/*
<AnalysisProtocolCollection>
*/
class MZIdentML_AnalysisProtocolCollection : public XMLOutputContainerItem {
public:
	MZIdentML_AnalysisProtocolCollection ( const VectorXMLOutputItemPtr& subItems = VectorXMLOutputItemPtr () );
};

// 8.  AnalysisSampleCollection           XMLOutputContainerItem
/*
<AnalysisSampleCollection>
*/
class MZIdentML_AnalysisSampleCollection : public XMLOutputContainerItem {
public:
	MZIdentML_AnalysisSampleCollection ( const VectorXMLOutputItemPtr& subItems );
};

// 9.  AnalysisSoftware                   XMLOutputContainerItem
/*
<AnalysisSoftware id="AS_mascot_server"
                  name="Mascot Server" version="2.2.03"
                  uri="http://www.matrixscience.com/search_form_select.html">
*/
class MZIdentML_AnalysisSoftware : public XMLOutputContainerItem {
public:
	MZIdentML_AnalysisSoftware ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& name, const std::string& version );
};

// 10. AnalysisSoftwareList               XMLOutputContainerItem
/*
<AnalysisSoftwareList>
*/
class MZIdentML_AnalysisSoftwareList : public XMLOutputContainerItem {
public:
	MZIdentML_AnalysisSoftwareList ( const VectorXMLOutputItemPtr& subItems );
};

// 11. AuditCollection                    XMLOutputContainerItem
/*
<AuditCollection>
*/
class MZIdentML_AuditCollection : public XMLOutputContainerItem {
public:
	MZIdentML_AuditCollection ( const VectorXMLOutputItemPtr& subItems );
};

// 12. BibliographicReference             XMLOutputAttrItem
/*
<BibliographicReference
	authors="David N. Perkins, Darryl J. C. Pappin, David M. Creasy, John S. Cottrell"
	editor=""
	id="10.1002/(SICI)1522-2683(19991201)20:18&lt;3551::AID-ELPS3551&gt;3.0.CO;2-2"
	name="Probability-based protein identification by searching sequence databases using mass spectrometry data"
	issue="18" pages="3551-3567" publication="Electrophoresis" volume="20" year="1999"
	publisher="Wiley VCH"
	title="Probability-based protein identification by searching sequence databases using mass spectrometry data"/>
*/
class MZIdentML_BibliographicReference : public XMLOutputAttrItem {
public:
	MZIdentML_BibliographicReference (
		const std::string& id,
		const std::string& authors = "",
		const std::string& title = "",
		const std::string& publication = "",
		const std::string& volume = "",
		const std::string& issue = "",
		const std::string& pages = "",
		const std::string& year = "",
		const std::string& doi = "",
		const std::string& editor = "",
		const std::string& publisher = "",
		const std::string& name = "" );
};

// 13. ContactRole                        XMLOutputContainerItem
/*
<ContactRole contact_ref="PERSON_DOC_OWNER">
*/
class MZIdentML_ContactRole : public XMLOutputContainerItem {
public:
	MZIdentML_ContactRole ( const VectorXMLOutputItemPtr& subItems, const std::string& contactRef );
};

// 14. Customizations                     XMLOutputValueItem
/*
<Customizations>No customisations</Customizations>
*/
class MZIdentML_Customizations : public XMLOutputValueItem {
public:
	MZIdentML_Customizations ( const std::string& customizations );
};

// 15. cv                                 XMLOutputAttrItem
//
// <cv id="PSI-MS" fullName="Proteomics Standards Initiative Mass Spectrometry Vocabularies"
//    uri="http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"
//    version="2.25.0"/>
// <cv id="UNIMOD"
//     fullName="UNIMOD"
//     uri="http://www.unimod.org/obo/unimod.obo"/>
// <cv id="UO"
//     fullName="UNIT-ONTOLOGY"
//     uri="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo"/>
//
class MZIdentML_CV : public XMLOutputAttrItem {
public:
	MZIdentML_CV ( const std::string& id, const std::string& fullName = "", const std::string& URI = "", const std::string& version = "" );
};

// 16. cvList                             XMLOutputContainerItem
/*
<cvList>
*/
class MZIdentML_CVList : public XMLOutputContainerItem {
public:
	MZIdentML_CVList ( const VectorXMLOutputItemPtr& subItems );
};

// 17. cvParam                            XMLOutputAttrItem
class MZIdentML_CVParam : public XMLOutputAttrItem {
public:
	MZIdentML_CVParam ( const std::string& accession, const std::string& name, const std::string& cvRef, const std::string& value = "", const std::string& unitAccession = "", const std::string& unitName = "", const std::string& unitCvRef = "" );
};

class MZIdentML_CVParam_AlternativeSingleLetterCodes : public MZIdentML_CVParam {
public:
	MZIdentML_CVParam_AlternativeSingleLetterCodes ( const std::string& value );
};

class MZIdentML_CVParam_ProteinDescription : public MZIdentML_CVParam {
public:
	MZIdentML_CVParam_ProteinDescription ( const std::string& value );
};

class MZIdentML_CVParam_ProteinTaxonomy : public MZIdentML_CVParam {
public:
	MZIdentML_CVParam_ProteinTaxonomy ( const std::string& value );
};

class MZIdentML_CVParam_ProteinTaxonomyID : public MZIdentML_CVParam {
public:
	MZIdentML_CVParam_ProteinTaxonomyID ( const std::string& value );
};

class MZIdentML_CVParam_Tolerance : public MZIdentML_CVParam {
	std::string getAccession ( const std::string& sign ) const;
	std::string getName ( const std::string& sign ) const;
	std::string getValue ( double tol, double sysError, const std::string& sign ) const;
	std::string getUnitAccession ( const std::string& units ) const;
	std::string getUnitName ( const std::string& units ) const;
public:
	MZIdentML_CVParam_Tolerance ( double tol, const std::string& units, double sysError, const std::string& sign );
};

class MZIdentML_CVParam_Unimod : public MZIdentML_CVParam {
public:
	MZIdentML_CVParam_Unimod ( int accNum, const std::string& name );
};

// 18. DatabaseFilters                    XMLOutputContainerItem
/*
<DatabaseFilters>
*/
class MZIdentML_DatabaseFilters : public XMLOutputContainerItem {
public:
	MZIdentML_DatabaseFilters ( const VectorXMLOutputItemPtr& subItems );
};

// 19. DatabaseName                       XMLOutputContainerItem
/*
<DatabaseName>
*/
class MZIdentML_DatabaseName : public XMLOutputContainerItem {
public:
	MZIdentML_DatabaseName ( const VectorXMLOutputItemPtr& subItems );
};

// 20. DataCollection                     XMLOutputContainerItem
/*
<DataCollection>
*/
class MZIdentML_DataCollection : public XMLOutputContainerItem {
public:
	MZIdentML_DataCollection ( const VectorXMLOutputItemPtr& subItems = VectorXMLOutputItemPtr () );
};

// 21. DBSequence                         XMLOutputContainerItem
/*
<DBSequence id="DBSeq_HSP7A_CAEEL" length="640" searchDatabase_ref="SDB_SwissProt" accession="HSP7A_CAEEL">
*/
class MZIdentML_DBSequence : public XMLOutputContainerItem {
public:
	MZIdentML_DBSequence ( const VectorXMLOutputItemPtr& subItems, const std::string& searchDatabaseRef, const std::string& accession, int length = 0, const std::string& name = "" );
};

// 22. Enzyme                             XMLOutputContainerItem
/*
<Enzyme id="ENZ_0" cTermGain="OH" nTermGain="H" missedCleavages="1" semiSpecific="0">
*/
class MZIdentML_Enzyme : public XMLOutputContainerItem {
public:
	MZIdentML_Enzyme ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& cTermGain = "", const std::string& minDistance = "", const std::string& missedCleavages = "", const std::string& nTermGain = "", const std::string& name = "", const std::string& semiSpecific = "" );
};

// 23. EnzymeName                         XMLOutputContainerItem
/*
<EnzymeName>
*/
class MZIdentML_EnzymeName : public XMLOutputContainerItem {
public:
	MZIdentML_EnzymeName ( const VectorXMLOutputItemPtr& subItems );
};

// 24. Enzymes                            XMLOutputContainerItem
/*
<Enzymes independent="0">
*/
class MZIdentML_Enzymes : public XMLOutputContainerItem {
public:
	MZIdentML_Enzymes ( const VectorXMLOutputItemPtr& subItems, const std::string& independent = "" );
};

// 25. FileFormat                         XMLOutputContainerItem
/*
<FileFormat>
*/
class MZIdentML_FileFormat : public XMLOutputContainerItem {
public:
	MZIdentML_FileFormat ( const VectorXMLOutputItemPtr& subItems );
};

// 26. Filter                             XMLOutputContainerItem
/*
<Filter>
*/
class MZIdentML_Filter : public XMLOutputContainerItem {
public:
	MZIdentML_Filter ( const VectorXMLOutputItemPtr& subItems );
};

// 27. FilterType                         XMLOutputContainerItem
/*
<FilterType>
*/
class MZIdentML_FilterType : public XMLOutputContainerItem {
public:
	MZIdentML_FilterType ( const VectorXMLOutputItemPtr& subItems );
};

// 28. FragmentArray                      XMLOutputAttrItem
/*
<FragmentArray values="214.8 286.1 342.8 444.1 814.1 " measure_ref="m_mz"/>
*/
class MZIdentML_FragmentArray : public XMLOutputAttrItem {
public:
	MZIdentML_FragmentArray ( const std::string& values, const std::string& measureRef );
};

// 29. Fragmentation                      XMLOutputContainerItem
/*
<Fragmentation>
*/
class MZIdentML_Fragmentation : public XMLOutputContainerItem {
public:
	MZIdentML_Fragmentation ( const VectorXMLOutputItemPtr& subItems );
};

// 30. FragmentationTable                 XMLOutputContainerItem
/*
<FragmentationTable>
*/
class MZIdentML_FragmentationTable : public XMLOutputContainerItem {
public:
	MZIdentML_FragmentationTable ( const VectorXMLOutputItemPtr& subItems );
};

// 31. FragmentTolerance                  XMLOutputContainerItem
/*
<FragmentTolerance>
*/
class MZIdentML_FragmentTolerance : public XMLOutputContainerItem {
public:
	MZIdentML_FragmentTolerance ( const VectorXMLOutputItemPtr& subItems );
};

// 32. Include                            XMLOutputContainerItem
/*
<Include>
*/
class MZIdentML_Include : public XMLOutputContainerItem {
public:
	MZIdentML_Include ( const VectorXMLOutputItemPtr& subItems );
};

// 33. Inputs                             XMLOutputContainerItem
/*
<Inputs>
*/
class MZIdentML_Inputs : public XMLOutputContainerItem {
public:
	MZIdentML_Inputs ( const VectorXMLOutputItemPtr& subItems );
};

// 34. InputSpectra                       XMLOutputAttrItem
/*
<InputSpectra spectraData_ref="SD_1"/>
*/
class MZIdentML_InputSpectra : public XMLOutputAttrItem {
public:
	MZIdentML_InputSpectra ( const std::string& spectraDataRef = "" );
};

// 35. InputSpectrumIdentifications       XMLOutputAttrItem
/*
<InputSpectrumIdentifications spectrumIdentificationList_ref="SIL_1"/>
*/
class MZIdentML_InputSpectrumIdentifications : public XMLOutputAttrItem {
public:
	MZIdentML_InputSpectrumIdentifications ( const std::string& spectrumIdentificationListRef );
};

// 36. IonType                            XMLOutputContainerItem
/*
<IonType index="1 2 3 4 8" charge="1">
*/
class MZIdentML_IonType : public XMLOutputContainerItem {
public:
	MZIdentML_IonType ( const VectorXMLOutputItemPtr& subItems, const std::string& charge, const std::string& index = "" );
};

// 37. MassTable                          XMLOutputContainerItem
/*
<MassTable id="MT" msLevel="1 2">
*/
class MZIdentML_MassTable : public XMLOutputContainerItem {
public:
	MZIdentML_MassTable ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& msLevel, const std::string& name = "" );
};

// 38. Measure                            XMLOutputContainerItem
/*
<Measure id="m_mz">
*/
class MZIdentML_Measure : public XMLOutputContainerItem {
public:
	MZIdentML_Measure ( const VectorXMLOutputItemPtr& subItems, const std::string& id );
};

// 39. Modification                       XMLOutputContainerItem
/*
<Modification location="11" residues="K" monoisotopicMassDelta="127.063324">
*/
class MZIdentML_Modification : public XMLOutputContainerItem {
public:
	MZIdentML_Modification ( const VectorXMLOutputItemPtr& subItems, int location, char residues, double monoisotopicMassDelta, double avgMassDelta = 0.0 );
};

// 40. ModificationParams                 XMLOutputContainerItem
/*
<ModificationParams>
*/
class MZIdentML_ModificationParams : public XMLOutputContainerItem {
public:
	MZIdentML_ModificationParams ( const VectorXMLOutputItemPtr& subItems );
};

// 41. MzIdentML                          XMLOutputContainerItem
/*
<MzIdentML xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
           xsi:schemaLocation="http://psidev.info/psi/pi/mzIdentML/1.1 ../../schema/mzIdentML1.1.0.xsd"
           xmlns="http://psidev.info/psi/pi/mzIdentML/1.1"
           id="example_mzidentml_1"
           version="1.1.0"
           creationDate="2009-08-18T17:59:55">
*/
class MZIdentML_MZIdentML : public XMLOutputContainerItem {
public:
	MZIdentML_MZIdentML ();
};

// 42. Organization                       XMLOutputContainerItem
/*
<Organization id="ORG_MSL" name="Matrix Science Limited">
*/
class MZIdentML_Organization : public XMLOutputContainerItem {
public:
	MZIdentML_Organization ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& name = "" );
};

// 43. Parent                             XMLOutputAttrItem
/*
<Parent organization_ref="ORG_MSL" />
*/
class MZIdentML_Parent : public XMLOutputAttrItem {
public:
	MZIdentML_Parent ( const std::string& organizationRef );
};

// 44. ParentTolerance                    XMLOutputContainerItem
/*
<ParentTolerance>
*/
class MZIdentML_ParentTolerance : public XMLOutputContainerItem {
public:
	MZIdentML_ParentTolerance ( const VectorXMLOutputItemPtr& subItems );
};

// 45. Peptide                            XMLOutputContainerItem
/*
<Peptide id="peptide_4_3">
*/
class MZIdentML_Peptide : public XMLOutputContainerItem {
public:
	MZIdentML_Peptide ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& name = "" );
};

// 46. PeptideEvidence                    XMLOutputAttrItem
/*
<PeptideEvidence
	peptide_ref="peptide_1_1"
	id="PE_1_1_HSP70_ONCMY_0"
	start="160"
	end="171"
	pre="K"
	post="L"
	isDecoy="false"
	dBSequence_ref="DBSeq_HSP70_ONCMY"/>
*/
class MZIdentML_PeptideEvidence : public XMLOutputAttrItem {
public:
	MZIdentML_PeptideEvidence ( const std::string& id, const std::string& dBSequenceRef, const std::string& peptideRef, int start = 0, int end = 0, const std::string& pre = "", const std::string& post = "", bool isDecoy = false, const std::string& frame = "", const std::string& name = "" );
};

// 47. PeptideEvidenceRef                 XMLOutputAttrItem
/*
<PeptideEvidenceRef peptideEvidence_ref="PE_1_1_HSP70_ONCMY_0"/>
*/
class MZIdentML_PeptideEvidenceRef : public XMLOutputAttrItem {
public:
	MZIdentML_PeptideEvidenceRef ( const std::string& peptideEvidenceRef );
};

// 48. PeptideHypothesis                  XMLOutputContainerItem
/*
<PeptideHypothesis peptideEvidence_ref="PE_3_1_HSP7D_MANSE_0">
*/
class MZIdentML_PeptideHypothesis : public XMLOutputContainerItem {
public:
	MZIdentML_PeptideHypothesis ( const VectorXMLOutputItemPtr& subItems, const std::string& peptideEvidenceRef );
};

// 49. PeptideSequence                    XMLOutputValueItem
/*
<PeptideSequence>MLMMMALLVPLVYTIK</PeptideSequence>
*/
class MZIdentML_PeptideSequence : public XMLOutputValueItem {
public:
	MZIdentML_PeptideSequence ( const std::string& peptideSequence );
};

// 50. Person                             XMLOutputContainerItem
/*
<Person id="PERSON_DOC_OWNER" firstName="" lastName="Some Person">
*/
class MZIdentML_Person : public XMLOutputContainerItem {
public:
	MZIdentML_Person ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& firstName = "", const std::string& midInitials = "", const std::string& lastName = "", const std::string& name = "" );
};

// 51. ProteinAmbiguityGroup              XMLOutputContainerItem
/*
<ProteinAmbiguityGroup id="PAG_hit_1">
*/
class MZIdentML_ProteinAmbiguityGroup : public XMLOutputContainerItem {
public:
	MZIdentML_ProteinAmbiguityGroup ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& name = "" );
};

// 52. ProteinDetection                   XMLOutputContainerItem
/*
<ProteinDetection
	id="PD_1"
	proteinDetectionProtocol_ref="PDP_MascotParser_1"
	proteinDetectionList_ref="PDL_1"
	activityDate="2009-08-18T18:01:21">
*/
class MZIdentML_ProteinDetection : public XMLOutputContainerItem {
public:
	MZIdentML_ProteinDetection ( const VectorXMLOutputItemPtr& subItems,
		const std::string& id,
		const std::string& proteinDetectionListRef,
		const std::string& proteinDetectionProtocolRef,
		const std::string& name = "" );
};

// 53. ProteinDetectionHypothesis         XMLOutputContainerItem
/*
<ProteinDetectionHypothesis
	id="PDH_HSP7D_MANSE_0"
	dBSequence_ref="DBSeq_HSP7D_MANSE"
	passThreshold="true">
*/
class MZIdentML_ProteinDetectionHypothesis : public XMLOutputContainerItem {
public:
	MZIdentML_ProteinDetectionHypothesis ( const VectorXMLOutputItemPtr& subItems,
		const std::string& id,
		const std::string& passThreshold,
		const std::string& dBSequenceRef = "",
		const std::string& name = "" );
};

// 54. ProteinDetectionList               XMLOutputContainerItem
/*
<ProteinDetectionList id="PDL_1">
*/
class MZIdentML_ProteinDetectionList : public XMLOutputContainerItem {
public:
	MZIdentML_ProteinDetectionList ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& name = "" );
};

// 55. ProteinDetectionProtocol           XMLOutputContainerItem
/*
<ProteinDetectionProtocol
	id="PDP_MascotParser_1"
	analysisSoftware_ref="AS_mascot_parser">
*/
class MZIdentML_ProteinDetectionProtocol : public XMLOutputContainerItem {
public:
	MZIdentML_ProteinDetectionProtocol ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& analysisSoftwareRef, const std::string& name = "" );
};

// 56. Provider                           XMLOutputContainerItem
/*
<Provider id="PROVIDER">
*/
class MZIdentML_Provider : public XMLOutputContainerItem {
public:
	MZIdentML_Provider ( const VectorXMLOutputItemPtr& subItems );
};

// 57. Residue                            XMLOutputAttrItem
/*
<Residue code="Y" mass="163.063329"/>
*/
class MZIdentML_Residue : public XMLOutputAttrItem {
public:
	MZIdentML_Residue ( char code, double mass );
};

// 58. Role                               XMLOutputContainerItem
/*
<Role>
*/
class MZIdentML_Role : public XMLOutputContainerItem {
public:
	MZIdentML_Role ( const VectorXMLOutputItemPtr& subItems );
};

// 59. Sample                             XMLOutputContainerItem
/*
<Sample id="sample2" name="name24">
*/
class MZIdentML_Sample : public XMLOutputContainerItem {
public:
	MZIdentML_Sample ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& name = "" );
};

// 60. SearchDatabase                     XMLOutputContainerItem
/*
<SearchDatabase
	location="file:///C:/inetpub/mascot/sequence/SwissProt/current/SwissProt_51.6.fasta"
	id="SDB_SwissProt"
	name="SwissProt"
	numDatabaseSequences="257964"
	numResidues="93947433"
	releaseDate="2011-03-01T21:32:52"
	version="SwissProt_51.6.fasta">
*/
class MZIdentML_SearchDatabase : public XMLOutputContainerItem {
public:
	MZIdentML_SearchDatabase ( const VectorXMLOutputItemPtr& subItems,
		const std::string& id,
		const std::string& location,
		const std::string& name = "",
		int numDatabaseSequences = -1,
		const std::string& numResidues = "",
		const std::string& releaseDate = "",
		const std::string& version = "" );
};

// 61. SearchDatabaseRef                  XMLOutputAttrItem
/*
<SearchDatabaseRef searchDatabase_ref="SDB_SwissProt"/>
*/
class MZIdentML_SearchDatabaseRef : public XMLOutputAttrItem {
public:
	MZIdentML_SearchDatabaseRef ( const std::string& searchDatabaseRef = "" );
};

// 62. SearchModification                 XMLOutputContainerItem
/*
<SearchModification fixedMod="false" massDelta="15.994915" residues="M">
*/
class MZIdentML_SearchModification : public XMLOutputContainerItem {
public:
	MZIdentML_SearchModification ( const VectorXMLOutputItemPtr& subItems,
		const std::string& fixedMod,
		const std::string& massDelta,
		const std::string& residues );
};

// 63. SearchType                         XMLOutputContainerItem
/*
<SearchType>
*/
class MZIdentML_SearchType : public XMLOutputContainerItem {
public:
	MZIdentML_SearchType ( const VectorXMLOutputItemPtr& subItems );
};

// 64. Seq                                XMLOutputValueItem
/*
<Seq>MSKNAIGIDLGTTYSCVGVFMHGKVEIIANDQGNRTTPSYVAFTDTERLIGDAAKNQVAMNPHNTVFDANRLIGRKFDDGSVQSDMKHWPFKVVNAGGGKPKVQVEYKGETKTFTPEEISSMVLVKMKETAEAFLGHAVKDAVITVPAYFNDSQRQATKDSGAIAGLNVLRIINEPTAAAIAYGLDKKGHGERNVLIFDLGGGTFDVSILTIEDGIFEVKSTAGDTHLGEDFDNRMVNHFVAEFKRNDKKDLASNPRALRRLRTACERAKRTLSSSSQASIEIDSLFEGIDFYTNITRARFEELCADLFRSTMDPVEKALRDAKMDKAQVHDIVLVGGSTRIPKVQKLLSDFFSGKELNKSINPDEAVAYGAAVQAAILSGDKSEAVQDLLLLDVAPLSLGIETAGGVMTALIKRNTTIPTKTSETFTTYSDNQPGVLIQVYEGERALTKDNNLLGKFELSGIPPAPRGVPQIEVTFDIDANGILNVSAQDKSTGKQNKITITNDKGRLSKDEIERMVQEAEKYKADDEAQKDRIAAKNALESYAFNMKQTIEDEKLKDKISEEDKKKIQEKCDETVRWLDGNQTAEKDEFEHRQKELESVCNPIITKLYQSAGGMPGGMPGGMPGGAPGAGSTGGGPTIEEVD</Seq>
*/
class MZIdentML_Seq : public XMLOutputValueItem {
public:
	MZIdentML_Seq ( const std::string& seq );
};

// 65. SequenceCollection                 XMLOutputContainerItem
/*
<SequenceCollection>
*/
class MZIdentML_SequenceCollection : public XMLOutputContainerItem {
public:
	MZIdentML_SequenceCollection ( const VectorXMLOutputItemPtr& subItems = VectorXMLOutputItemPtr () );
};

// 66. SiteRegexp                         XMLOutputValueItem
/*
<SiteRegexp><![CDATA[(?<=[KR])(?!P)]]></SiteRegexp>
*/
class MZIdentML_SiteRegexp : public XMLOutputValueItem {
public:
	MZIdentML_SiteRegexp ( const std::string& regexp );
};

// 67. SoftwareName                       XMLOutputContainerItem
/*
<SoftwareName>
*/
class MZIdentML_SoftwareName : public XMLOutputContainerItem {
public:
	MZIdentML_SoftwareName ( const VectorXMLOutputItemPtr& subItems );
};

// 68. SourceFile                         XMLOutputContainerItem
/*
<SourceFile location="file:///../data/F001350.dat" id="SF_1">
*/
class MZIdentML_SourceFile : public XMLOutputContainerItem {
public:
	MZIdentML_SourceFile ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& location, const std::string& name = "" );
};

// 69. SpecificityRules                   XMLOutputContainerItem
/*
<SpecificityRules>
*/
class MZIdentML_SpecificityRules : public XMLOutputContainerItem {
public:
	MZIdentML_SpecificityRules ( const VectorXMLOutputItemPtr& subItems );
};

// 70. SpectraData                        XMLOutputContainerItem
/*
<SpectraData location="file:///dyckall.asc" id="SD_1">
*/
class MZIdentML_SpectraData : public XMLOutputContainerItem {
public:
	MZIdentML_SpectraData ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& location, const std::string& name = "" );
};

// 71. SpectrumIdentification             XMLOutputContainerItem
/*
<SpectrumIdentification
	id="SI"
	spectrumIdentificationProtocol_ref="SIP"
	spectrumIdentificationList_ref="SIL_1"
	activityDate="2008-07-09T18:33:47">
*/
class MZIdentML_SpectrumIdentification : public XMLOutputContainerItem {
public:
	MZIdentML_SpectrumIdentification ( const VectorXMLOutputItemPtr& subItems,
		const std::string& id,
		const std::string& spectrumIdentificationListRef,
		const std::string& spectrumIdentificationProtocolRef,
		const std::string& activityDate = "",
		const std::string& name = "" );
};

// 72. SpectrumIdentificationItem         XMLOutputContainerItem
/*
<SpectrumIdentificationItem
	id="SII_1_1"
	calculatedMassToCharge="671.869886"
	chargeState="2"
	experimentalMassToCharge="671.9"
	peptide_ref="peptide_1_1"
	rank="1"
	passThreshold="true"
	massTable_ref="MT"
	sample_ref="sample1">
*/
class MZIdentML_SpectrumIdentificationItem : public XMLOutputContainerItem {
public:
	MZIdentML_SpectrumIdentificationItem ( const VectorXMLOutputItemPtr& subItems, const std::string& id,
		double calculatedMassToCharge, int chargeState, double experimentalMassToCharge,
		bool passThreshold, int rank, const std::string& peptideRef = "",
		const std::string& massTableRef = "", const std::string& sampleRef = "", const std::string& calculatedPI = "",
		const std::string& name = "" );
};

// 73. SpectrumIdentificationItemRef      XMLOutputAttrItem
/*
<SpectrumIdentificationItemRef spectrumIdentificationItem_ref="SII_1_1"/>
*/
class MZIdentML_SpectrumIdentificationItemRef : public XMLOutputAttrItem {
public:
	MZIdentML_SpectrumIdentificationItemRef ( const std::string& spectrumIdentificationItemRef );
};

// 74. SpectrumIdentificationList         XMLOutputContainerItem
/*
<SpectrumIdentificationList id="SIL_1" numSequencesSearched="71412">
*/
class MZIdentML_SpectrumIdentificationList : public XMLOutputContainerItem {
public:
	MZIdentML_SpectrumIdentificationList ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& numSequencesSearched = "", const std::string& name = "" );
};

// 75. SpectrumIdentificationProtocol     XMLOutputContainerItem
/*
<SpectrumIdentificationProtocol id="SIP" analysisSoftware_ref="AS_mascot_server">
*/
class MZIdentML_SpectrumIdentificationProtocol : public XMLOutputContainerItem {
public:
	MZIdentML_SpectrumIdentificationProtocol ( const VectorXMLOutputItemPtr& subItems, const std::string& id, const std::string& analysisSoftwareRef, const std::string& name = "" );
};

// 76. SpectrumIdentificationResult       XMLOutputContainerItem
/*
<SpectrumIdentificationResult id="SIR_1" spectrumID="query=1" spectraData_ref="SD_1">
*/
class MZIdentML_SpectrumIdentificationResult : public XMLOutputContainerItem {
public:
	MZIdentML_SpectrumIdentificationResult ( const VectorXMLOutputItemPtr& subItems,
		const std::string& id,
		const std::string& spectrumID,
		const std::string& spectraDataRef,
		const std::string& name = "" );
};

// 77. SpectrumIDFormat                   XMLOutputContainerItem
/*
<SpectrumIDFormat>
*/
class MZIdentML_SpectrumIDFormat : public XMLOutputContainerItem {
public:
	MZIdentML_SpectrumIDFormat ( const VectorXMLOutputItemPtr& subItems );
};

// 78. SubSample                          XMLOutputAttrItem
/*
<SubSample sample_ref="sample1"/>
*/
class MZIdentML_SubSample : public XMLOutputAttrItem {
public:
	MZIdentML_SubSample ( const std::string& sampleRef );
};

// 79. Threshold                          XMLOutputContainerItem
/*
<Threshold>
*/
class MZIdentML_Threshold : public XMLOutputContainerItem {
public:
	MZIdentML_Threshold ( const VectorXMLOutputItemPtr& subItems );
};

// 80. userParam                          XMLOutputAttrItem
/*
<userParam name="contact fax" value="+44 (0)20 7224 1344"/>
*/
class MZIdentML_userParam : public XMLOutputAttrItem {
public:
	MZIdentML_userParam ( const std::string& name, const std::string& value = "" );
};

#endif /* ! __lu_mzidentml_h */
