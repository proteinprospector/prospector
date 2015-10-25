/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pep_xml.cpp                                                *
*                                                                             *
*  Created    : December 8th 2010                                             *
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
#include <lg_io.h>
#include <lg_string.h>
#include <lg_time.h>
#include <lu_aa_info.h>
#include <lu_mass_elem.h>
#include <lu_pep_xml.h>

using std::string;
using std::ostream;
using std::ostringstream;

/*
http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v114.xsd
http://www.matrixscience.com/xmlns/schema/pepXML_v18/index.html
*/

/*
tandem K <?xml-stylesheet type="text/xsl" href="pepXML_std.xsl"?>

interact-prob <?xml-stylesheet type="text/xsl" href="/regis/sbeams4/nobackup/edeutsch/iPRG2011/submission/XTK1/interact-prob.pep.xsl"?>

SpectraST <?xml-stylesheet type="text/xsl" href="/proteomics/sw//tpp-trunk-5223/schema/pepXML_std.xsl"?>

SEQUEST <?xml-stylesheet type="text/xsl" href="http://regis-web.systemsbiology.net/pepXML_std.xsl"?>

OMSSA <!DOCTYPE msms_pipeline_analysis PUBLIC "-//NCBI//pepXML/EN" "pepXML.dtd">

MASCOT <?xml-stylesheet type="text/xsl" href="/proteomics/sw/tpp-4.0.0/schema/pepXML_std.xsl"?>

ascore <?xml-stylesheet type="text/xsl" href="http://regis-web.systemsbiology.net/pepXML_std.xsl"?>
*/

void printPepXMLHeader ( ostream& os )
{
	printXMLHeader ( os );
	printXMLStylesheet ( os, "text/xsl", "/prospector/html/pepXML_std.xsl" );
}

/*

Example:

tandem K
--------

<msms_pipeline_analysis
	date="2010:12:01:01:12:03"
	summary_xml="/regis/sbeams4/nobackup/edeutsch/iPRG2011/submission/XTK1/D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed.pep.xml"
	xmlns="http://regis-web.systemsbiology.net/pepXML"
	xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
	xsi:schemaLocation="http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v116.xsd">

interact-prob
-------------

<msms_pipeline_analysis
	date="2010-12-02T23:00:48"
	summary_xml="/regis/sbeams4/nobackup/edeutsch/iPRG2011/submission/XTK1/interact-prob.pep.xml"
	xmlns="http://regis-web.systemsbiology.net/pepXML"
	xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
	xsi:schemaLocation="http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v116.xsd">

SpectraST
---------

<msms_pipeline_analysis
	date="2010-11-29T23:22:34"
	xmlns="http://regis-web.systemsbiology.net/pepXML"
	xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
	xsi:schemaLocation="http://regis-web.systemsbiology.net/pepXML /proteomics/sw//tpp-trunk-5223/schema/pepXML_v116.xsd"
	summary_xml="D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed.pep.xml">

SEQUEST
-------

<msms_pipeline_analysis
	date="2010-12-02T18:33:35"
	xmlns="http://regis-web.systemsbiology.net/pepXML"
	xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
	xsi:schemaLocation="http://regis-web.systemsbiology.net/pepXML /net/pr/vol1/ProteomicsResource/bin/TPP/bin/20101102-TPP-v4.4.1/schema/pepXML_v115.xsd"
	summary_xml="./D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed.pep.xml">

OMSSA
-----

<msms_pipeline_analysis
	date="2010-11-29T23:28:55"
	summary_xml="D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed.pep.xml">

MASCOT
------

<msms_pipeline_analysis
	date="2010-12-02T11:57:20"
	xmlns="http://regis-web.systemsbiology.net/pepXML"
	xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
	xsi:schemaLocation="http://regis-web.systemsbiology.net/pepXML /proteomics/sw/tpp-4.0.0/schema/pepXML_v115.xsd"
	summary_xml="D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed.pep.xml">

ascore
------

<msms_pipeline_analysis
	date="2010-02-16T12:44:08"
	xmlns="http://regis-web.systemsbiology.net/pepXML"
	xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
	xsi:schemaLocation="http://regis-web.systemsbiology.net/pepXML c:\inetpub\wwwrootpepXML_v114.xsd"
	summary_xml="test.pep.xml">

Example:

<msms_pipeline_analysis
	name=" xs:string [0..1] ?"			Summary name (currently not used)
	date=" xs:dateTime [1] ?"			Date pepXML file was written
	summary_xml=" xs:string [1] ?">		Full path self reference

*/

PepXMLMSMSPipelineAnalysis::PepXMLMSMSPipelineAnalysis ( const string& name, const string& summaryXML ) :
	XMLOutputContainerItem ( "msms_pipeline_analysis" )
{
	init ( name, summaryXML );
}
PepXMLMSMSPipelineAnalysis::PepXMLMSMSPipelineAnalysis ( const VectorXMLOutputItemPtr& subItems, const string& name, const string& summaryXML ) :
	XMLOutputContainerItem ( "msms_pipeline_analysis", subItems )
{
	init ( name, summaryXML );
}
void PepXMLMSMSPipelineAnalysis::init ( const string& name, const string& summaryXML )
{
	attr.push_back ( makePairStringString ( "date", genCurrentXSDDateAndTimeString () ) );
	if ( !name.empty () )	attr.push_back ( makePairStringString ( "name", name ) );
	attr.push_back ( makePairStringString ( "summary_xml", summaryXML ) );
	attr.push_back ( makePairStringString ( "xmlns",				"http://regis-web.systemsbiology.net/pepXML" ) );
	attr.push_back ( makePairStringString ( "xmlns:xsi",			"http://www.w3.org/2001/XMLSchema-instance" ) );
	attr.push_back ( makePairStringString ( "xsi:schemaLocation",	"http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v116.xsd" ) );
}

/*
tandem K - not present

interact-prob
-------------

<analysis_summary analysis="database_refresh" time="2010-12-02T23:00:59"/>
<analysis_summary analysis="peptideprophet" time="2010-12-02T23:00:49">
	...peptideprophet_summary
</analysis_summary>

SpectraST - not present

SEQUEST - <analysis_summary analysis="database_refresh" time="2010-12-02T22:47:34"/>

OMSSA - <analysis_summary analysis="database_refresh" time="2010-11-29T23:40:51"/>

MASCOT - <analysis_summary analysis="database_refresh" time="2010-12-02T23:33:55"/>

ascore - not present

*/

PepXMLAnalysisSummary::PepXMLAnalysisSummary ( const string& analysis, const string& t ) :
	XMLOutputContainerItem ( "analysis_summary" )
{
	init ( analysis, t );
}
PepXMLAnalysisSummary::PepXMLAnalysisSummary ( const VectorXMLOutputItemPtr& subItems, const string& analysis, const string& t ) :
	XMLOutputContainerItem ( "analysis_summary", subItems )
{
	init ( analysis, t );
}
void PepXMLAnalysisSummary::init ( const string& analysis, const string& t )
{
	attr.push_back ( makePairStringString ( "analysis", analysis ) );
	attr.push_back ( makePairStringString ( "time", t ) );
}

// msms_run_summary
/*
tandem K
--------

<msms_run_summary
	base_name="/regis/sbeams4/nobackup/edeutsch/iPRG2011/submission/XTK1/D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed"
	search_engine="X! Tandem (k-score)"
	msManufacturer="Thermo Scientific"
	msModel="LTQ Orbitrap XL"
	msIonization="nanoelectrospray"
	msMassAnalyzer="orbitrap"
	msDetector="inductive detector"
	raw_data_type="raw"
	raw_data=".mzXML">

interact-prob
-------------

<msms_run_summary
	base_name="/regis/sbeams4/nobackup/edeutsch/iPRG2011/submission/XTK1/D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed"
	search_engine="X! Tandem (k-score)"
	msManufacturer="Thermo Scientific"
	msModel="LTQ Orbitrap XL"
	msIonization="nanoelectrospray"
	msMassAnalyzer="orbitrap"
	msDetector="inductive detector"
	raw_data_type="raw"
	raw_data=".mzXML">

SpectraST
---------

<msms_run_summary
	base_name="/regis/sbeams4/nobackup/edeutsch/iPRG2011/submission/SST1/D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed"
	msManufacturer="Thermo Scientific"
	msModel="LTQ Orbitrap XL"
	msIonization="nanoelectrospray"
	msMassAnalyzer="orbitrap"
	msDetector="inductive detector"
	raw_data_type=".mzXML"
	raw_data=".mzXML">

SEQUEST
-------

<msms_run_summary
	base_name="./D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed"
	msManufacturer="Thermo Scientific" msModel="LTQ Orbitrap XL"
	msIonization="nanoelectrospray"
	msMassAnalyzer="orbitrap"
	msDetector="inductive detector"
	raw_data_type="raw"
	raw_data=".mzXML">

OMSSA
-----

<msms_run_summary
	base_name="D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed.pep.xml"
	raw_data_type="raw"
	raw_data=".mzXML">

MASCOT
------

<msms_run_summary
	base_name="/regis/sbeams4/nobackup/edeutsch/iPRG2011/submission/MAS1/D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed"
	msManufacturer="Thermo Scientific"
	msModel="LTQ Orbitrap XL"
	msIonization="nanoelectrospray"
	msMassAnalyzer="orbitrap"
	msDetector="inductive detector"
	raw_data_type="raw"
	raw_data=".mzXML">

ascore
------

<msms_run_summary
	base_name="test"
	raw_data_type="raw"
	raw_data=".mzXML">

<msms_run_summary
	base_name=" xs:string [1] ?"					full path file name of mzXML
	raw_data_type=" xs:string [1] ?"				raw data type extension
	raw_data=" xs:string [1] ?"
	msManufacturer=" xs:string [0..1] ?"			ABI / SCIEX,Bruker Daltonics,IonSpec,Shimadzu,ThermoFinnigan,
													Waters
	msModel=" xs:string [0..1] ?"					Bruker Daltonics:
														microFlex,ultraFlex TOF/TOF,ultraFlex,autoFlex TOF/TOF,
														autoFlex II,OmniFlex,HCTPlus,HCT,esquire6000,esquire4000,
														APEX IV,APEX-Q,BioTOF II,BioTOF Q,microTOFLC
													ThermoFinnigan:
														TEMPUS TOF,ELEMENT2,Neptune,TRACE DSQ,Surveyor MSQ,
														TSQ QUANTUM,LCQ Classic,LCQ Deca XP Plus,LCQ Advantage,
														PolarisQ,DELTAplusAdvantage,DELTAplusXP,MAT253,MAT900XP,
														MAT900XP Trap,MAT95XP,MAT95XP Trap,TRITON
													Waters:
														Q-Tof micro,Quattro micro,Quattro UItima,GCT,
														AutoSpec Ultima NT,M@LDI L,M@LDI LR,Q-Tof Ultima,
														IsoPrime,IsoProbe,IsoProbe T,Platform ICP,NG-5400
													Applied Biosystems:
														API 150EX,API 150EX Prep,API 2000,API 3000,API 4000,
														Proteomics Solution 1,Q TRAP,4000 Q TRAP,QSTAR,
														4700 Proteomic Analyzer,SymBiot I,SymBiot XVI,
														Voyager-DE PRO,Voyager-DE STR
													IonSpec:
														HiResMALDI,HiResESI,Ultima,Explorer,OMEGA,OMEGA-2001
	msIonization=" xs:string [0..1] ?"				ESI,MALDI,EI,CI,FAB,APCI									
	msMassAnalyzer=" xs:string [0..1] ?"			TOF,Quadrupole,Quadrupole Ion Trap,FT-ICR,Magnetic Sector
	msDetector=" xs:string [0..1] ?"> [1..*] ?		Channeltron,Daly,EMT,Faraday Cup,Microchannel plate
*/

PepXMLMSMSRunSummary::PepXMLMSMSRunSummary ( const string& baseName, const string& rawDataType, const string& rawData, const string& msManufacturer, const string& msModel, const string& msIonization, const string& msMassAnalyzer, const string& msDetector ) :
	XMLOutputContainerItem ( "msms_run_summary" )
{
	init ( baseName, rawDataType, rawData, msManufacturer, msModel, msIonization, msMassAnalyzer, msDetector );
}
PepXMLMSMSRunSummary::PepXMLMSMSRunSummary ( const VectorXMLOutputItemPtr& subItems, const string& baseName, const string& rawDataType, const string& rawData, const string& msManufacturer, const string& msModel, const string& msIonization, const string& msMassAnalyzer, const string& msDetector ) :
	XMLOutputContainerItem ( "msms_run_summary", subItems )
{
	init ( baseName, rawDataType, rawData, msManufacturer, msModel, msIonization, msMassAnalyzer, msDetector );
}
void PepXMLMSMSRunSummary::init ( const string& baseName, const string& rawDataType, const string& rawData, const string& msManufacturer, const string& msModel, const string& msIonization, const string& msMassAnalyzer, const string& msDetector )
{
	attr.push_back ( makePairStringString ( "base_name", baseName ) );
	attr.push_back ( makePairStringString ( "raw_data_type", rawDataType ) );
	attr.push_back ( makePairStringString ( "raw_data", rawData ) );
	if ( !msManufacturer.empty () )	attr.push_back ( makePairStringString ( "msManufacturer", msManufacturer ) );
	if ( !msModel.empty () )		attr.push_back ( makePairStringString ( "msModel", msModel ) );
	if ( !msIonization.empty () )	attr.push_back ( makePairStringString ( "msIonization", msIonization ) );
	if ( !msMassAnalyzer.empty () )	attr.push_back ( makePairStringString ( "msMassAnalyzer", msMassAnalyzer ) );
	if ( !msDetector.empty () )		attr.push_back ( makePairStringString ( "msDetector", msDetector ) );
}
void PepXMLMSMSRunSummary::setBaseName ( const string& baseName )
{
	attr [0] = makePairStringString ( "base_name", baseName );
}

/*
<sample_enzyme									Defines the net cleavage specificity of an enzyme, chemical reagent, or a mixture of these, for mass spectrometry purposes
	name=" xs:string (length >= 1) [1] ?"		Controlled code name for the enzyme that can be referred to by applications
	description=" xs:string [0..1] ?"			Free text to describe alternative names, special conditions, etc.
	fidelity=" xs:string (value comes from list: {'specific'|'semispecific'|'nonspecific'}) [0..1] ?"		Semispecific means that at least one end of a pepide must conform to the cleavage specificity, (unless the peptide was at the terminus of the parent sequence). Nonspecific means that neither end of a peptide must conform to the cleavage specificity.
	independent=" xs:boolean [0..1] ?"> [1] ?		If there are multiple specificities and independent is true, then a single peptide cannot exhibit one specificity at one terminus and a different specificity at the other. If independent is false, then a single peptide can exhibit mixed specificities.
	<specificity									Component cleavage specificity. Must be at least one specificity unless enzymeType:fidelity is nonspecific 	
		sense=" xs:string (value comes from list: {'C'|'N'}) [1] ?"				Defines whether cleavage occurs on the C-terminal or N-terminal side of the residue(s) listed in cut
		min_spacing=" xs:nonNegativeInteger [0..1] ?"							Minimum separation between adjacent cleavages
		cut=" xs:string (pattern = [A,C-I,K-N,P-T,VWY]+) (length >= 1) [1] ?"	One or more 1-letter residue codes. Enzyme cleaves on the sense side of the residue(s) listed in cut unless one of the residues listed in no_cut is adjacent to the potential cleavage site
		no_cut=" xs:string (length >= 0) [0..1] ?"/> [0..*] ?					Zero or more 1-letter residue codes. Enzyme cleaves on the sense side of the residue(s) listed in cut unless one of the residues listed in no_cut is adjacent to the potential cleavage site
</sample_enzyme>

Examples:

<sample_enzyme name="lysc" fidelity="semispecific">
	<specificity cut="K" no_cut="P" sense="C"/>
</sample_enzyme>

<sample_enzyme name="trypsin">
	<specificity cut="KR" no_cut="P" sense="C"/>
</sample_enzyme>
*/

PepXMLSampleEnzyme::PepXMLSampleEnzyme ( const VectorXMLOutputItemPtr& subItems, const string& name, const string& description, const string& fidelity, int independent ) :
	XMLOutputContainerItem ( "sample_enzyme", subItems )
{
	attr.push_back ( makePairStringString ( "name", name ) );
	if ( !description.empty () )	attr.push_back ( makePairStringString ( "description", description ) );
	if ( !fidelity.empty () )		attr.push_back ( makePairStringString ( "fidelity", fidelity ) );
	if ( independent != -1 )		attr.push_back ( makePairStringString ( "independent", independent == 0 ? "0" : "1" ) );
}
PepXMLSpecificity::PepXMLSpecificity ( const string& sense, int minSpacing, const string& cut, const string& noCut ) :
	XMLOutputAttrItem ( "specificity" )
{
	attr.push_back ( makePairStringString ( "sense",	sense ) );
	if ( minSpacing != -1 ) attr.push_back ( makePairStringString ( "min_spacing", gen_itoa ( minSpacing ) ) );
	attr.push_back ( makePairStringString ( "cut",		cut ) );
	if ( !noCut.empty () )	attr.push_back ( makePairStringString ( "no_cut", noCut ) );
}

/*

<search_summary										Database search settings
	base_name=" xs:string [1] ?"					Full path location of mzXML file for this search run (without the .mzXML extension)
	search_engine=" engineType [1] ?"				SEQUEST, Mascot, COMET, etc
	precursor_mass_type=" massType [1] ?"			average or monoisotopic
	fragment_mass_type=" massType [1] ?"			average or monoisotopic
	out_data_type=" xs:string [1] ?"				Format of file storing the runner up peptides (if not present in pepXML)
	out_data=" xs:string [1] ?"						runner up search hit data type extension (e.g. .tgz)
	search_id=" positiveInt [1] ?"> [1..*] ?		matches id in search hit

Examples:

tandem K
--------

<search_summary
	base_name="/regis/sbeams4/nobackup/edeutsch/iPRG2011/submission/XTK1/D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed"
	search_engine="X! Tandem (k-score)"
	precursor_mass_type="monoisotopic"
	fragment_mass_type="monoisotopic"
	search_id="1">

interact-prob
-------------

<search_summary
	base_name="/regis/sbeams4/nobackup/edeutsch/iPRG2011/submission/XTK1/D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed"
	search_engine="X! Tandem (k-score)"
	precursor_mass_type="monoisotopic"
	fragment_mass_type="monoisotopic"
	search_id="1">

SpectraST
---------

<search_summary
	base_name="/regis/sbeams4/nobackup/edeutsch/iPRG2011/submission/SST1/D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed"
	search_engine="SpectraST"
	precursor_mass_type="monoisotopic"
	fragment_mass_type="monoisotopic"
	out_data_type="out"
	out_data=".tgz"
	search_id="1">

SEQUEST
-------

<search_summary
	base_name="./D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed"
	search_engine="SEQUEST"
	precursor_mass_type="monoisotopic"
	fragment_mass_type="monoisotopic"
	out_data_type="out"
	out_data=".tgz"
	search_id="1">

OMSSA
-----

<search_summary
	base_name="D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed.pep.xml"
	search_engine="OMSSA"
	precursor_mass_type="monoisotopic"
	fragment_mass_type="monoisotopic"
	out_data_type="n/a"
	out_data="n/a"
	search_id="1">

MASCOT
------

<search_summary base_name="/regis/sbeams4/nobackup/edeutsch/iPRG2011/submission/MAS1/D100930_yeast_SCX10S_rak_ft8E_pc_01-fixed"
	search_engine="MASCOT"
	precursor_mass_type="monoisotopic"
	fragment_mass_type="monoisotopic"
	out_data_type="out"
	out_data=".tgz" 
	search_id="1">

ascore
------

<search_summary
	base_name="test"
	search_engine="SEQUEST"
	precursor_mass_type="average"
	fragment_mass_type="average"
	out_data_type="out"
	out_data=".tgz"
	search_id="1">
*/

PepXMLSearchSummary::PepXMLSearchSummary ( const string& baseName, const string& searchEngine, const string& precursorMassType, const string& fragmentMassType, int searchID, const string& outDataType, const string& outData ) :
	XMLOutputContainerItem ( "search_summary" )
{
	init ( baseName, searchEngine, precursorMassType, fragmentMassType, searchID, outDataType, outData );
}
PepXMLSearchSummary::PepXMLSearchSummary ( const VectorXMLOutputItemPtr& subItems, const string& baseName, const string& searchEngine, const string& precursorMassType, const string& fragmentMassType, int searchID, const string& outDataType, const string& outData ) :
	XMLOutputContainerItem ( "search_summary", subItems )
{
	init ( baseName, searchEngine, precursorMassType, fragmentMassType, searchID, outDataType, outData );
}
void PepXMLSearchSummary::init ( const string& baseName, const string& searchEngine, const string& precursorMassType, const string& fragmentMassType, int searchID, const string& outDataType, const string& outData )
{
	attr.push_back ( makePairStringString ( "base_name", baseName ) );
	attr.push_back ( makePairStringString ( "search_engine", searchEngine ) );
	attr.push_back ( makePairStringString ( "precursor_mass_type", precursorMassType ) );
	attr.push_back ( makePairStringString ( "fragment_mass_type", fragmentMassType ) );
	attr.push_back ( makePairStringString ( "search_id", gen_itoa ( searchID ) ) );
	if ( !outDataType.empty () )	attr.push_back ( makePairStringString ( "out_data_type", outDataType ) );
	if ( !outData.empty () )		attr.push_back ( makePairStringString ( "out_data", outData ) );
}
void PepXMLSearchSummary::setBaseName ( const string& baseName )
{
	attr [0] = makePairStringString ( "base_name", baseName );
}

/*
<search_database
	local_path=" xs:string [1] ?"										Full path address of database on local computer
	URL=" xs:string [0..1]"
	database_name=" xs:string [0..1]"
	orig_database_url=" xs:string [0..1]"
	database_release_date=" xs:dateTime [0..1]"
	database_release_identifier=" xs:string [0..1]"
	size_in_db_entries=" xs:nonNegativeInteger [0..1]"
	size_of_residues=" xs:nonNegativeInteger [0..1]"
	type=" xs:string (value comes from list: {'AA'|'NA'}) [1] ?"/> [0..1]	Database type (AA=amino acid, NA=nucleic acid)

Example:

<search_database
	local_path="/regis/sbeams4/nobackup/edeutsch/iPRG2011/dbase/iPRG2011_TargetDecoy.fasta"
	type="AA"/>
<search_database
	local_path="/regis/sbeams4/nobackup/edeutsch/iPRG2011/dbase/iPRG2011_TargetDecoy.fasta"
	size_in_db_entries="13332"
	type="AA"/>
*/

PepXMLSearchDatabase::PepXMLSearchDatabase ( const string& localPath, const string& type, int sizeInDBEntries, const string& url, const string& databaseName, const string& origDatabaseURL, const string& databaseReleaseDate, const string& databaseReleaseIdentifier, int sizeOfResidues ) :
	XMLOutputAttrItem ( "search_database" )
{
	attr.push_back ( makePairStringString ( "local_path",			localPath ) );
	attr.push_back ( makePairStringString ( "type",					type ) );
	attr.push_back ( makePairStringString ( "size_in_db_entries",	gen_itoa ( sizeInDBEntries ) ) );

	if ( !url.empty () )						attr.push_back ( makePairStringString ( "URL",							url ) );
	if ( !databaseName.empty () )				attr.push_back ( makePairStringString ( "database_name",				databaseName ) );
	if ( !origDatabaseURL.empty () )			attr.push_back ( makePairStringString ( "orig_database_url",			origDatabaseURL ) );
	if ( !databaseReleaseDate.empty () )		attr.push_back ( makePairStringString ( "database_release_date",		databaseReleaseDate ) );
	if ( !databaseReleaseIdentifier.empty () )	attr.push_back ( makePairStringString ( "database_release_identifier",	databaseReleaseIdentifier ) );

	if ( sizeOfResidues != -1 ) attr.push_back ( makePairStringString ( "size_of_residues", gen_itoa ( sizeOfResidues ) ) );
}

/*
<enzymatic_search_constraint
	enzyme=" xs:string [1]"
	max_num_internal_cleavages=" xs:nonNegativeInteger [1] ?"
	min_number_termini=" xs:nonNegativeInteger [1] ?"/> [0..1] ?

tandem K - <enzymatic_search_constraint enzyme="trypsin_k" max_num_internal_cleavages="2" min_number_termini="1" />

interact-prob - <enzymatic_search_constraint enzyme="trypsin_k" max_num_internal_cleavages="2" min_number_termini="1"/>

SpectraST - not present

SEQUEST - <enzymatic_search_constraint enzyme="Trypsin_K" max_num_internal_cleavages="2" min_number_termini="1"/>

OMSSA - <enzymatic_search_constraint enzyme="Lys-C" max_num_internal_cleavages="3" min_number_termini="1"/>

MASCOT - <enzymatic_search_constraint enzyme="lysc" max_num_internal_cleavages="0" min_number_termini="1"/>

ascore - not present
*/

PepXMLEnzymaticSearchConstraint::PepXMLEnzymaticSearchConstraint ( const string& enzyme, int maxNumInternalCleavages, int minNumberTermini ) :
	XMLOutputAttrItem ( "enzymatic_search_constraint" )
{
	attr.push_back ( makePairStringString ( "enzyme",						enzyme ) );
	attr.push_back ( makePairStringString ( "max_num_internal_cleavages",	gen_itoa ( maxNumInternalCleavages ) ) );
	attr.push_back ( makePairStringString ( "min_number_termini",			gen_itoa ( minNumberTermini ) ) );
}

/*
<aminoacid_modification								Modified aminoacid, static or variable
	aminoacid=" xs:string [1]"
	massdiff=" xs:string [1] ?"						Mass difference with respect to unmodified aminoacid, must begin with either + (nonnegative) or - [e.g. +1.05446 or -2.3342]
	mass=" xs:float [1] ?"							Mass of modified aminoacid
	variable=" xs:string [1] ?"						Y if both modified and unmodified aminoacid could be present in the dataset, N if only modified aminoacid can be present
	peptide_terminus=" xs:string [0..1] ?"			whether modification can reside only at protein terminus (specified 'n', 'c', or 'nc')
	symbol=" aa_symbolType [0..1] ?"				Special symbol used by search engine to designate this modification
	binary=" xs:string [0..1] ?"					Y if each peptide must have only modified or unmodified aminoacid, N if a peptide may contain both modified and unmodified aminoacid
	description=" xs:string [0..1]"/> [0..*] 

Example:

<aminoacid_modification
	aminoacid="B"
	mass="114.534940"
	massdiff="114.534940"
	variable="N"/>
*/

PepXMLAminoAcidModification::PepXMLAminoAcidModification ( const string& aminoacid, double massdiff, const string& variable, const string& peptideTerminus, const string& description, const string& binary, const string& symbol ) :
	XMLOutputAttrItem ( "aminoacid_modification" )
{
	double mass = AAInfo::getInfo ().getMonoisotopicMass ( aminoacid [0] ) + massdiff;
	attr.push_back ( makePairStringString ( "aminoacid",aminoacid ) );
	attr.push_back ( makePairStringString ( "massdiff",	gen_ftoa ( massdiff, "%.5f" ) ) );
	attr.push_back ( makePairStringString ( "mass",		gen_ftoa ( mass, "%.5f" ) ) );
	attr.push_back ( makePairStringString ( "variable",	variable ) );
	if ( !peptideTerminus.empty () )attr.push_back ( makePairStringString ( "peptide_terminus",	peptideTerminus ) );
	if ( !description.empty () )	attr.push_back ( makePairStringString ( "description",		description ) );
	if ( !binary.empty () )			attr.push_back ( makePairStringString ( "binary",			binary ) );
	if ( !symbol.empty () )			attr.push_back ( makePairStringString ( "symbol",			symbol ) );
}

/*

<terminal_modification							Modification to the N or C terminus, static or variable
	terminus=" xs:string [1] ?"					n for N-terminus, c for C-terminus
	massdiff=" xs:string [1] ?"					Mass difference with respect to unmodified terminus
	mass=" xs:float [1] ?"						Mass of modified terminus
	variable=" xs:string [1] ?"					Y if both modified and unmodified terminus could be present in the dataset, N if only modified terminus can be present
	symbol=" term_symbolType [0..1] ?"			Special symbol used by search engine to designate this modification
	protein_terminus=" xs:string [1] ?"			whether modification can reside only at protein terminus (specified n or c)
	description=" xs:string [0..1]"/> [0..*] ?

Example:

*/

PepXMLTerminalModification::PepXMLTerminalModification ( const string& terminus, double massdiff, const string& variable, const string& proteinTerminus, const string& description, const string& symbol ) :
	XMLOutputAttrItem ( "terminal_modification" )
{
	double mass = massdiff + formula_to_monoisotopic_mass ( terminus == "n" ? "H" : "O H" );
	attr.push_back ( makePairStringString ( "terminus",terminus ) );
	attr.push_back ( makePairStringString ( "massdiff",	gen_ftoa ( massdiff, "%.5f" ) ) );
	attr.push_back ( makePairStringString ( "mass",		gen_ftoa ( mass, "%.5f" ) ) );
	attr.push_back ( makePairStringString ( "variable",	variable ) );
	if ( !proteinTerminus.empty () )attr.push_back ( makePairStringString ( "protein_terminus",	proteinTerminus ) );
	if ( !description.empty () )	attr.push_back ( makePairStringString ( "description",		description ) );
	if ( !symbol.empty () )			attr.push_back ( makePairStringString ( "symbol",			symbol ) );
}

/*

Example:

<parameter
	name="ITOLU"
	value="Da"/>
*/

PepXMLParameter::PepXMLParameter ( const string& name, const string& value ) :
	XMLOutputAttrItem ( "parameter" )
{
	attr.push_back ( makePairStringString ( "name", name ) );
	attr.push_back ( makePairStringString ( "value", value ) );
}

/*
<analysis_timestamp									Reference for analysis applied to current run (time corresponds with analysis_summary/@time, id corresponds with analysis_result/@id)
	time=" xs:dateTime [1] ?"						Date of analysis
	analysis=" xs:string [1] ?"						Analysis name
	id=" positiveInt [1] ?"> [0..*] ?				Unique identifier for each type of analysis
	Allow any elements from any namespace (lax validation). [0..1]
</analysis_timestamp>

Example:	SEQUEST, OMSSA, interact-prob, Mascot

<analysis_timestamp
	analysis="database_refresh"
	time="2010-12-02T23:33:55"
	id="1">
*/

PepXMLAnalysisTimestamp::PepXMLAnalysisTimestamp ( const VectorXMLOutputItemPtr& subItems, const string& analysis, const string& time, const string& id ) :
	XMLOutputContainerItem ( "analysis_timestamp", subItems )
{
	attr.push_back ( makePairStringString ( "analysis", analysis ) );
	attr.push_back ( makePairStringString ( "time", genCurrentXSDDateAndTimeString () ) );
	attr.push_back ( makePairStringString ( "id", id ) );
}

/*
<spectrum_query
	spectrum=" xs:string [1]"
	start_scan=" xs:unsignedInt [1] ?"				first scan number integrated into MS/MS spectrum
	end_scan=" xs:unsignedInt [1] ?"				last scan number integrated into MS/MS spectrum
	precursor_neutral_mass=" xs:float [1]"
	assumed_charge=" xs:nonNegativeInteger [1] ?"	Precursor ion charge used for search
	search_specification=" xs:string [0..1] ?"		Search constraint applied specifically to this query
	index=" positiveInt [1] ?">						Unique identifier

	retention_time_sec								retention time associated with start_scan
	activation_method								Activation or fragmentation method: ETD, ECD, CID, etc

Example:

<spectrum_query
	spectrum="Shelly_AlphaYeast_9.0639.0639.4"
	start_scan="639"
	end_scan="639"
	precursor_neutral_mass="2875.1264"
	assumed_charge="4"
	index="1">
*/

PepXMLSpectrumQuery::PepXMLSpectrumQuery ( const VectorXMLOutputItemPtr& subItems, const string& spectrum,
	int startScan, int endScan, double precursorNeutralMass, int assumedCharge,
	const string& searchSpecification, int index, double retentionTimeSec,
	const string& activationMethod ) :

	XMLOutputContainerItem ( "spectrum_query", subItems )
{
	attr.push_back ( makePairStringString ( "spectrum", spectrum ) );
	attr.push_back ( makePairStringString ( "start_scan", gen_itoa ( startScan ) ) );
	attr.push_back ( makePairStringString ( "end_scan", gen_itoa ( endScan ) ) );
	attr.push_back ( makePairStringString ( "precursor_neutral_mass", gen_ftoa ( precursorNeutralMass, "%.4f" ) ) );
	attr.push_back ( makePairStringString ( "assumed_charge", gen_itoa ( assumedCharge ) ) );
	if ( !searchSpecification.empty () )	attr.push_back ( makePairStringString ( "search_specification", searchSpecification ) );
	attr.push_back ( makePairStringString ( "index", gen_itoa ( index ) ) );
	if ( retentionTimeSec != -1.0 )	attr.push_back ( makePairStringString ( "retention_time_sec", gen_ftoa ( retentionTimeSec, "%.4f" ) ) );
	if ( !activationMethod.empty () )		attr.push_back ( makePairStringString ( "activation_method", activationMethod ) );
}

/*
<search_result
	search_id=" positiveInt [0..1] ?"> [0..*]

	<search_hit ...
*/

PepXMLSearchResult::PepXMLSearchResult ( const VectorXMLOutputItemPtr& subItems, int searchID ) :

	XMLOutputContainerItem ( "search_result", subItems )
{
	if ( searchID != -1 )	attr.push_back ( makePairStringString ( "search_id", gen_itoa ( searchID ) ) );
}

/*
<search_hit
	hit_rank=" positiveInt [1]"
	peptide=" xs:string [1] ?"							Peptide aminoacid sequence (with no indicated modifications)
	peptide_prev_aa=" xs:string [0..1] ?"				Aminoacid preceding peptide (- if none)
	peptide_next_aa=" xs:string [0..1] ?"				Aminoacid following peptide (- if none)
	protein=" xs:string [1]"
	num_tot_proteins=" xs:unsignedInt [1] ?"			Number of unique proteins in search database containing peptide
	num_matched_ions=" xs:nonNegativeInteger [0..1] ?"	Number of peptide fragment ions found in spectrum
	tot_num_ions=" xs:nonNegativeInteger [0..1] ?"		Number of peptide fragment ions predicted for peptide
	calc_neutral_pep_mass=" xs:float [1]"
	massdiff=" xs:string [1] ?"							Mass(precursor ion) - Mass(peptide)
	num_tol_term=" xs:nonNegativeInteger [0..1] ?"		Number of peptide termini consistent with cleavage by sample enzyme
	num_missed_cleavages=" xs:integer [0..1] ?"			Number of sample enzyme cleavage sites internal to peptide
	is_rejected=" xs:nonNegativeInteger (value comes from list: {'0'|'1'}) [0..1] ?"	Potential use in future for user manual validation (0 or 1)
	protein_descr=" xs:string [0..1] ?"					Extracted from search database
	calc_pI=" xs:string [0..1]"
	protein_mw=" xs:double [0..1]"> [0..*] ?

	<alternative_protein ...
	<modification_info ...
	<search_score ...
	<analysis_result ...
	<parameter ...
</search_hit>

Example:

<search_hit
	hit_rank="1"
	peptide="RFSSEHPDPVETSIPEQAAEIAEELSK"
	peptide_prev_aa="R"
	peptide_next_aa="Q"
	protein="ORFP:YGL120C"
	num_tot_proteins="1"
	num_matched_ions="36"
	tot_num_ions="104"
	calc_neutral_pep_mass="3075.4074"
	massdiff="+0.021570"
	num_tol_term="2"
	num_missed_cleavages="0"
	is_rejected="0">
*/

PepXMLSearchHit::PepXMLSearchHit ( const VectorXMLOutputItemPtr& subItems, int hitRank, const string& peptide,
	const string& peptidePrevAA, const string& peptideNextAA, const string& protein, int numTotProteins,
	int numMatchedIons, int totNumIons, double calcNeutralPepMass, double massdiff, int numTolTerm,
	int numMissedCleavages, int isRejected, const string& proteinDescr, double calcPI, double proteinMW ) :

	XMLOutputContainerItem ( "search_hit", subItems )
{
	attr.push_back ( makePairStringString ( "hit_rank", gen_itoa ( hitRank ) ) );
	attr.push_back ( makePairStringString ( "peptide", peptide ) );
	string prevAA = peptidePrevAA.substr ( peptidePrevAA.length () - 1 );
	string nextAA = peptideNextAA.substr ( 0, 1 );
	if ( prevAA == "*" ) prevAA = "";
	if ( nextAA == "*" ) nextAA = "";
	if ( !prevAA.empty () ) attr.push_back ( makePairStringString ( "peptide_prev_aa", prevAA ) );
	if ( !nextAA.empty () ) attr.push_back ( makePairStringString ( "peptide_next_aa", nextAA ) );
	attr.push_back ( makePairStringString ( "protein", protein ) );
	attr.push_back ( makePairStringString ( "num_tot_proteins", gen_itoa ( numTotProteins ) ) );
	if ( numMatchedIons != -1 )	attr.push_back ( makePairStringString ( "num_matched_ions", gen_itoa ( numMatchedIons ) ) );
	if ( totNumIons != -1 )		attr.push_back ( makePairStringString ( "tot_num_ions", gen_itoa ( totNumIons ) ) );
	attr.push_back ( makePairStringString ( "calc_neutral_pep_mass", gen_ftoa ( calcNeutralPepMass, "%.4f" ) ) );
	attr.push_back ( makePairStringString ( "massdiff", gen_ftoa ( massdiff, "%.5f" ) ) );
	if ( numTolTerm != -1 )			attr.push_back ( makePairStringString ( "num_tol_term", gen_itoa ( numTolTerm ) ) );
	if ( numMissedCleavages != -1 )	attr.push_back ( makePairStringString ( "num_missed_cleavages", gen_itoa ( numMissedCleavages ) ) );
	if ( isRejected != -1 )			attr.push_back ( makePairStringString ( "is_rejected", gen_itoa ( isRejected ) ) );
	if ( !proteinDescr.empty () )	attr.push_back ( makePairStringString ( "protein_descr", proteinDescr ) );
	if ( calcPI != 0.0 )		attr.push_back ( makePairStringString ( "calc_pI", gen_ftoa ( calcPI, "%.1f" ) ) );
	if ( proteinMW != 0.0 )		attr.push_back ( makePairStringString ( "protein_mw", gen_ftoa ( proteinMW, "%.1f" ) ) );
}

/*
<alternative_protein
	protein=" xs:string [1]"
	protein_descr=" xs:string [0..1]"
	num_tol_term=" xs:nonNegativeInteger [0..1]"
	protein_mw=" xs:double [0..1]"/> [0..*] ?
	peptide_prev_aa=" xs:string [0..1] ?"				Aminoacid preceding peptide (- if none)
	peptide_next_aa=" xs:string [0..1] ?"				Aminoacid following peptide (- if none)

Example:

<alternative_protein
	protein="DECOY_sp|P38873|KOG1_YEAST"
	protein_descr="Decoy sequence"
	num_tol_term="1"
	peptide_prev_aa="P"
	peptide_next_aa="Q"/>
*/

PepXMLAlternativeProtein::PepXMLAlternativeProtein ( const string& protein, const string& proteinDescr,
	int numTolTerm, double proteinMW, const string& peptidePrevAA, const string& peptideNextAA ) :

	XMLOutputAttrItem ( "alternative_protein" )
{
	attr.push_back ( makePairStringString ( "protein", protein ) );
	if ( !proteinDescr.empty () )	attr.push_back ( makePairStringString ( "protein_descr", proteinDescr ) );
	if ( numTolTerm != -1 )			attr.push_back ( makePairStringString ( "num_tol_term", gen_itoa ( numTolTerm ) ) );
	if ( proteinMW != 0.0 )			attr.push_back ( makePairStringString ( "protein_mw", gen_ftoa ( proteinMW, "%.1f" ) ) );
	attr.push_back ( makePairStringString ( "peptide_prev_aa", peptidePrevAA ) );
	attr.push_back ( makePairStringString ( "peptide_next_aa", peptideNextAA ) );
}

/*
<modification_info
	mod_nterm_mass=" xs:double [0..1] ?"				Mass of modified N terminus
	mod_cterm_mass=" xs:double [0..1] ?"				Mass of modified C terminus
	modified_peptide=" xs:string [0..1] ?"> [0..1] ?	Peptide sequence (with indicated modifications)

	mod_aminoacid_mass...

</modification_info>

<mod_aminoacid_mass
	position=" xs:nonNegativeInteger [1] ?"				modified aminoacid position in peptide [ranging from 1 to peptide length]
	mass=" xs:double [1] ?"/> [0..*]					modified mass of aminoacid

Example:

<modification_info>
	<mod_aminoacid_mass position="1" mass="147.035404"/>
	<mod_aminoacid_mass position="5" mass="115.026936"/>
</modification_info>

<modification_info
	mod_nterm_mass="111.03258"
	modified_peptide="Q[111]SPHKKRN[115]SNNK">
	<mod_aminoacid_mass position="8" mass="115.02693"/>
</modification_info>

*/

PepXMLModificationInfo::PepXMLModificationInfo ( const VectorXMLOutputItemPtr& subItems, double modNTermMass, double modCTermMass, const string& modifiedPeptide ) :

	XMLOutputContainerItem ( "modification_info", subItems )
{
	static double hAdjust = formula_to_monoisotopic_mass ( "H" );
	static double ohAdjust = formula_to_monoisotopic_mass ( "O H" );

	if ( modNTermMass != 0.0 )		attr.push_back ( makePairStringString ( "mod_nterm_mass", gen_ftoa ( modNTermMass + hAdjust, "%.5f" ) ) );
	if ( modCTermMass != 0.0 )		attr.push_back ( makePairStringString ( "mod_cterm_mass", gen_ftoa ( modCTermMass + ohAdjust, "%.5f" ) ) );
	if ( !modifiedPeptide.empty () )attr.push_back ( makePairStringString ( "modified_peptide", modifiedPeptide ) );
}

PepXMLModAminoacidMass::PepXMLModAminoacidMass ( int position, double mass ) :

	XMLOutputAttrItem ( "mod_aminoacid_mass" )
{
	attr.push_back ( makePairStringString ( "position", gen_itoa ( position ) ) );
	attr.push_back ( makePairStringString ( "mass", gen_ftoa ( mass, "%.5f" ) ) );
}

/*

<search_score> nameValueType </search_score> [0..*]

Example:

<search_score
	name="ionscore"
	value="10.27"/>
*/

PepXMLSearchScore::PepXMLSearchScore ( const string& name, double value ) :

	XMLOutputAttrItem ( "search_score" )
{
	attr.push_back ( makePairStringString ( "name", name ) );
	ostringstream ost;
	genPrintSigFig ( ost, value, 3 );
	attr.push_back ( makePairStringString ( "value", ost.str () ) );
}
