/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_srch_form.cpp                                              *
*                                                                             *
*  Created    : December 9th 2004                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lg_time.h>
#include <lgen_math.h>
#include <lgen_file.h>
#ifdef MYSQL_DATABASE
#include <ld_init.h>
#endif
#include <lu_table.h>
#include <lu_srch_form.h>
#include <lu_prog.h>
#include <lu_html.h>
#include <lu_cgi_val.h>
#include <lu_fas_sp.h>
#include <lu_check_db.h>
#include <lu_mut_mtrx.h>
#include <lu_fas_enz.h>
#include <lu_getfil.h>
#include <lu_usermod.h>
#include <lu_pk_filter.h>
#include <lu_param_list.h>
#include <lu_const_mod.h>
#include <lu_form_valid.h>
#include <lu_get_link.h>
using std::ostream;
using std::copy;
using std::back_inserter;
using std::string;
using std::vector;
using std::endl;

class FormItemProtLowMass : public FormItemText {
public:
	FormItemProtLowMass ( FormValidatingJavascript* fvj, const string& ms = "", const string& value = "1000" ) :
		FormItemText ( "Protein MW (Da)", "allman.htm#prot_MW", getName ( ms ), 6, 10, value, fvj->addPositiveFloatingPointValidator (getName ( ms ), "Start Protein MW Range") ) {}
	static string getName ( const string& ms = "" ) { return ms + string ( "prot_low_mass" ); }
	static string getSuffix ( const string& ms )
	{
		if ( ms == "ms_" ) return " (MS)";
		else if ( ms == "msms_" ) return " (MS/MS)";
		else return "";
	}
};

class FormItemProtHighMass : public FormItemText {
public:
	FormItemProtHighMass ( FormValidatingJavascript* fvj, const string& ms = "", const string& value = "125000" ) :
		FormItemText ( "to", "", getName ( ms ), 6, 10, value, fvj->addPositiveFloatingPointValidator (getName ( ms ), "End Protein MW Range") ) {}
	static string getName ( const string& ms = "" ) { return ms + string ( "prot_high_mass" ); }
};

class FormItemFullMWRange : public FormItemCheckbox {
public:
	FormItemFullMWRange ( const string& ms = "", bool value = true ) :
		FormItemCheckbox ( "All", "", getName ( ms ), value ) {}
	static string getName ( const string& ms = "" ) { return ms + string ( "full_mw_range" ); }
};

class FormItemLowPI : public FormItemText {
public:
	FormItemLowPI ( FormValidatingJavascript* fvj ) :
		FormItemText ( "Protein pI", "allman.htm#prot_pi", getName (), 6, 6, "3.0", fvj->addPositiveFloatingPointValidator (getName (), "Start pI Range") ) {}
	static string getName () { return "low_pi"; }
};

class FormItemHighPI : public FormItemText {
public:
	FormItemHighPI ( FormValidatingJavascript* fvj ) :
		FormItemText ( "to", "", getName (), 6, 6, "10.0", fvj->addPositiveFloatingPointValidator (getName (), "End pI Range") ) {}
	static string getName () { return "high_pi"; }
};

class FormItemFullPIRange : public FormItemCheckbox {
public:
	FormItemFullPIRange () :
		FormItemCheckbox ( "All", "", getName (), true ) {}
	static string getName () { return "full_pi_range"; }
};

void ResourceLinkTable::getNumRows ()
{
	numRows = gen_round_up_int_divide ( resourceList.size (), numColumns );
}
void ResourceLinkTable::printHTML ( ostream& os ) const
{
	os << "<div id=\"header\">" << endl;
		for ( int i = 0 ; i < numRows ; i++ ) {
			for ( int j = 0 ; j < numColumns ; j++ ) {
				int index = ( i * numColumns ) + j;
					if ( index < resourceList.size () ) {
						const string& f = resourceList [index].first;
						const string& s = resourceList [index].second;
						os << "<a ";
						//if ( f [0] != '.' ) os << "target=\"_blank\" ";	// ../mshome shouldn't launch a new window
						os << "href=\"";
						os << f;
						os << "\">";
						os << s;
						os << "</a>";
						os << endl;
						if ( index != resourceList.size () - 1 && index != numColumns - 1 ) os << " | ";
					}
			}
			os << "<br />" << endl;
		}
	os << "</div>" << endl;
}
OldResourceLinkTable::OldResourceLinkTable ()
{
	getResourceList ();
}
void OldResourceLinkTable::getResourceList ( const string& s )
{
	resourceList.push_back ( PairStringString ( "../mshome.htm", "Home" ) );
	resourceList.push_back ( getPair ( "msfitstandard", "MS-Fit" ) );
	resourceList.push_back ( getPair ( "mstagstandard", "MS-Tag" ) );
	resourceList.push_back ( getPair ( "msseq", "MS-Seq" ) );
	resourceList.push_back ( getPair ( "mspattern", "MS-Pattern" ) );
	resourceList.push_back ( getPair ( "msbridgestandard", "MS-Bridge" ) );

	resourceList.push_back ( getPair ( "msdigest", "MS-Digest" ) );
	resourceList.push_back ( getPair ( "msproduct", "MS-Product" ) );
	resourceList.push_back ( getPair ( "mscomp", "MS-Comp" ) );
	resourceList.push_back ( getPair ( "dbstat", "DB-Stat" ) );
	resourceList.push_back ( getPair ( "msisotope", "MS-Isotope" ) );
	resourceList.push_back ( getPair ( "mshomology", "MS-Homology" ) );
	resourceList.push_back ( getPair ( "msviewer", "MS-Viewer" ) );

#ifdef MYSQL_DATABASE
	vector <PairStringString> pss;
	pss.push_back ( PairStringString ( "select_project", "1" ) );
	resourceList.push_back ( getPair ( "batchtag", "Batch-Tag", pss ) );
	resourceList.push_back ( getPair ( "batchtagweb", "Batch-Tag Web" ) );
	resourceList.push_back ( getPair ( "search_compare", "Search Compare", pss ) );
	resourceList.push_back ( getPair ( "results_management", "Results Management", pss ) );
	resourceList.push_back ( getPair ( "search_table", "Search Table" ) );
	resourceList.push_back ( getPair ( "logout", "Logout" ) );
	resourceList.push_back ( getPair ( "delete_batchtag_cookie", "Delete Batch-Tag Cookie" ) );
	resourceList.push_back ( getPair ( "delete_search_compare_cookie", "Delete Search Compare Cookie" ) );
#endif

	numColumns = 13;
	getNumRows ();
}
PairStringString ResourceLinkTable::getPair ( const string& form, const string& name )
{
	string s = ProgramLink::getURLStart ( "msform" );
	s += "?form=";
	s += form;
	return PairStringString ( s, name );
}
PairStringString ResourceLinkTable::getPair ( const string& form, const string& name, const vector <PairStringString>& p )
{
	string s = ProgramLink::getURLStart ( "msform" );
	s += "?form=";
	s += form;
	for ( int i = 0 ; i < p.size () ; i++ ) {
		s += "&amp;";
		s += p [i].first;
		s += "=";
		s += escapeURL ( p [i].second );
	}
	return PairStringString ( s, name );
}
FormItemOutputType::FormItemOutputType ( bool tabDelimitedTextOption ) :
	FormItemSelect ( "Output", "allman.htm#output", getName (), tabDelimitedTextOption ? options2 : options, "HTML" )
{
}
FormItemOutputType::FormItemOutputType ( bool dummy1, bool dummy2 ) :
	FormItemSelect ( "Output", "allman.htm#output", getName (), options3, "HTML" )
{
}
const char* FormItemDNAReadingFrame::options [] = { "1", "2", "3", "4", "5", "6", 0 };
FormItemDNAReadingFrame::FormItemDNAReadingFrame () :
	FormItemSelect ( "Reading Frame (DNA Only)", "allman.htm#frame", getName (), options, "1" )
{
}

FormItemTaxonomy::FormItemTaxonomy () :
	FormItemSelectMultiple ( "Taxonomy", "allman.htm#taxonomy", getName (), StringVector (), StringVector (), 6 )
{
	value.push_back ( "All" );
	StringVector speciesList = TaxonomyGroupNames::instance ().getList ();
	StringVector speciesList2 = TaxonomyNames::instance ().getList ();
	copy ( speciesList.begin (), speciesList.end (), back_inserter ( value ) );
	copy ( speciesList2.begin (), speciesList2.end (), back_inserter ( value ) );
	select.push_back ( "All" );
}
FormItemTaxonomyRemove::FormItemTaxonomyRemove () :
	FormItemCheckbox ( "Remove Selected Taxonomy", "allman.htm#taxonomy_remove", getName (), false )
{
}
FormItemTaxonomyNames::FormItemTaxonomyNames () :
	FormItemTextArea ( "Taxonomy<br />Names", "allman.htm#taxonomy_names", getName (), 3, 30, StringVector () )
{
}
ProteinMWForm::ProteinMWForm ( FormValidatingJavascript* fvj, const VectorConstParameterListPtr& params, const string& ms ) :
	fvj ( fvj ),
	ms ( ms )
{
	create ( params );
}
void ProteinMWForm::createItems ()
{
	formItemMap [FormItemProtLowMass::getName (ms)]	= new FormItemProtLowMass ( fvj, ms );
	formItemMap [FormItemProtHighMass::getName (ms)]= new FormItemProtHighMass ( fvj, ms );
	formItemMap [FormItemFullMWRange::getName (ms)]	= new FormItemFullMWRange ( ms, ms == "ms_" ? false : true );
}
void ProteinMWForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemProtLowMass::getName (ms)]->setValue ( p );
		formItemMap [FormItemProtHighMass::getName (ms)]->setValue ( p );
		formItemMap [FormItemFullMWRange::getName (ms)]->setValue ( p );
	}
}
void ProteinMWForm::printHTML ( ostream& os )
{
	formItemMap [FormItemProtLowMass::getName (ms)]->printHTML ( os );
	formItemMap [FormItemProtHighMass::getName (ms)]->printHTML ( os );
	formItemMap [FormItemFullMWRange::getName (ms)]->printHTML ( os );
}
ProteinPIForm::ProteinPIForm ( FormValidatingJavascript* fvj, const VectorConstParameterListPtr& params ) :
	fvj ( fvj )
{
	create ( params );
}
void ProteinPIForm::createItems ()
{
	formItemMap [FormItemLowPI::getName ()]			= new FormItemLowPI (fvj);
	formItemMap [FormItemHighPI::getName ()]		= new FormItemHighPI (fvj);
	formItemMap [FormItemFullPIRange::getName ()]	= new FormItemFullPIRange ();
}
void ProteinPIForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemLowPI::getName ()]->setValue ( p );
		formItemMap [FormItemHighPI::getName ()]->setValue ( p );
		formItemMap [FormItemFullPIRange::getName ()]->setValue ( p );
	}
}
void ProteinPIForm::printHTML ( ostream& os )
{
	formItemMap [FormItemLowPI::getName ()]->printHTML ( os );
	formItemMap [FormItemHighPI::getName ()]->printHTML ( os );
	formItemMap [FormItemFullPIRange::getName ()]->printHTML ( os );
}
FormItemNames::FormItemNames () :
	FormItemTextArea ( "Names", "allman.htm#names", getName (), 3, 30, StringVector () )
{
}
FormItemAccessionNums::FormItemAccessionNums () :
	FormItemTextArea ( "Accession<br />Numbers", "allman.htm#accession_numbers", getName (), 3, 10, StringVector () )
{
}
FormItemAddAccessionNumbers::FormItemAddAccessionNumbers () :
	FormItemTextArea ( "Additional<br />Accession<br />Numbers", "allman.htm#add_accession_numbers", getName (), 3, 10, StringVector () )
{
}

FormItemSaveParameters::FormItemSaveParameters () :
	FormItemCheckbox ( "Save Settings", "", getName (), false )
{
}

const char* FormItemDataPlotted::options [] = { "Centroid", "Raw and Centroid", "Raw", 0 };
const char* FormItemAllowNonSpecific::options [] = { "at 0 termini", "at 1 termini", "at 2 termini", "at N termini", "at C termini", "N termini-1=D", 0 };

const char* FormItemParentMassConvert::optionsMS [] = { "monoisotopic", "average", 0 };
const char* FormItemParentMassConvert::optionsMSMS [] = { "monoisotopic", "average", "Par(mi)Frag(av)", "Par(av)Frag(mi)", 0 };
const char* FormItemMSParentMassToleranceUnits::options [] = { "Da", "%", "ppm", "mmu", 0 };
const char* FormItemMSMSParentMassToleranceUnits::options [] = { "Da", "%", "ppm", "mmu", 0 };
const char* FormItemFragmentMassesToleranceUnits::options [] = { "Da", "%", "ppm", "mmu", 0 };
const char* FormItemMSMSPrecursorCharge::options [] = { "8", "7", "6", "5", "4", "3", "2", "1", "Automatic", 0 };
const char* FormItemCompMaskType::options [] = { "OR", "AND", 0 };
const char* FormItemMSMSPrecursorChargeRange::options [] = { "1", "2", "3", "4", "5",
															"1 2 3",
															"2 3",
															"1 2 3 4",
															"2 3 4",
															"1 2 3 4 5",
															"2 3 4 5",
															0 };

class FormItemDNAFrameTranslation : public FormItemSelect {
	static const char* options [];
public:
	FormItemDNAFrameTranslation () :
		FormItemSelect ( "DNA Frame Translation", "allman.htm#frame", getName (), options, "3" ) {}
	static string getName () { return "dna_frame_translation"; }
};
const char* FormItemDNAFrameTranslation::options [] = { "6", "3", "-3", "1", "-1", 0 };

FormItemNTerminusAALimit::FormItemNTerminusAALimit ( FormValidatingJavascript* fvj ) :
	FormItemText ( "N Term AA Limit", "allman.htm#n_term_aa_limit", getName (), 5, 5, "", fvj->addPositiveNonZeroIntegerAllowBlankValidator ( getName (), "N Term AA Limit" ) )
{
}

class FormItemResultsFromFile : public FormItemCheckbox {
public:
	FormItemResultsFromFile ( bool set = false ) :
		FormItemCheckbox ( "Hits", "allman.htm#saving", getName (), set ) {}
	static string getName () { return "results_from_file"; }
};

class FormItemInputProgramName : public FormItemSelect {
	static const char* options [];
public:
	FormItemInputProgramName () :
		FormItemSelect ( "From", "", getName (), options, "msfit" ) {}
	static string getName () { return "input_program_name"; }
};
const char* FormItemInputProgramName::options [] = { "msfit", "mspattern", "mshomology", "mstag", "msseq", 0 };

class FormItemInputFilename : public FormItemText {
public:
	FormItemInputFilename () :
		FormItemText ( "Name", "", getName (), 8, 100, "lastres" ) {}
	static string getName () { return "input_filename"; }
};

const char* FormItemOutputType::options [] = { "HTML", "XML", 0 };
const char* FormItemOutputType::options2 [] = { "HTML", "XML", "Tab delimited text", 0 };
const char* FormItemOutputType::options3 [] = { "HTML", "XML", "Tab delimited text", "mgf", 0 };

const char* FormItemExpectationCalculationMethod::options [] = { "None", "Linear Tail Fit"/*, "Method of Moments", "Closed Form Max Likelihood"*/, 0 };

const char* FormItemMSMSPkFilter::options [] = { "Max MSMS Pks", "Max MSMS Pks / 100 Da", "Unprocessed MSMS", 0 };
const char* FormItemMSMSPkFilter::options2 [] = { "Max MSMS Pks", "Max MSMS Pks / 100 Da", 0 };

const char* FormItemRawType::options [] = { "MS Precursor", "MS Full Scan", "MS/MS", "Quantitation", 0 };

class FormItemResultsToFile : public FormItemCheckbox {
public:
	FormItemResultsToFile ( bool checked ) :
		FormItemCheckbox ( "Hits to file", "allman.htm#saving", getName (), checked ) {}
	static string getName () { return "results_to_file"; }
};

class FormItemOutputFilename : public FormItemText {
public:
	FormItemOutputFilename ( FormValidatingJavascript* fvj, bool btag ) :
		FormItemText ( btag ? "Results Name" : "Name", "", getName (), btag ? 30 : 12, 58, btag ? "" : "lastres", fvj->addFilenameValidator ( getName (), btag ? "Results Name" : "Name" ) ) {}
	static string getName () { return "output_filename"; }
};
SearchResultsForm::SearchResultsForm ( const VectorConstParameterListPtr& params, bool set ) :
	set ( set )
{
	create ( params );
}
void SearchResultsForm::createItems ()
{
	formItemMap [FormItemResultsFromFile::getName ()]	= new FormItemResultsFromFile ( set );
	formItemMap [FormItemInputProgramName::getName ()]	= new FormItemInputProgramName ();
	formItemMap [FormItemInputFilename::getName ()]		= new FormItemInputFilename ();
}
void SearchResultsForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemResultsFromFile::getName ()]->setValue ( p );
		formItemMap [FormItemInputProgramName::getName ()]->setValue ( p );
		formItemMap [FormItemInputFilename::getName ()]->setValue ( p );
	}
}
void SearchResultsForm::printHTML ( ostream& os )
{
	formItemMap [FormItemResultsFromFile::getName ()]->printHTML ( os );
	formItemMap [FormItemInputProgramName::getName ()]->printHTML ( os );
	formItemMap [FormItemInputFilename::getName ()]->printHTML ( os );
}
void SearchResultsForm::printHTMLFAIndex ( ostream& os )
{
	formItemMap [FormItemResultsFromFile::getName()]->printHTMLHidden ( os );
	formItemMap [FormItemInputProgramName::getName ()]->printHTML ( os );
	formItemMap [FormItemInputFilename::getName ()]->printHTML ( os );
}
SaveResultsForm::SaveResultsForm ( const VectorConstParameterListPtr& params, FormValidatingJavascript* fvj, bool defaultSave, bool tabDelimitedTextOption ) :
	fvj ( fvj ),
	defaultSave ( defaultSave ),
	tabDelimitedTextOption ( tabDelimitedTextOption )
{
	create ( params );
}
void SaveResultsForm::createItems ()
{
	if ( !defaultSave ) {
		formItemMap [FormItemOutputType::getName ()]	= new FormItemOutputType ( tabDelimitedTextOption );
		formItemMap [FormItemResultsToFile::getName ()]	= new FormItemResultsToFile ( defaultSave ? true : false );
	}
	formItemMap [FormItemOutputFilename::getName ()]= new FormItemOutputFilename ( fvj, defaultSave );
}
void SaveResultsForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [FormItemOutputFilename::getName ()]->setValue ( p );
	}
}
void SaveResultsForm::printHTML ( ostream& os )
{
	if ( !defaultSave ) {
		formItemMap [FormItemOutputType::getName ()]->printHTML ( os );
		formItemMap [FormItemResultsToFile::getName ()]->printHTML ( os );
	}
	formItemMap [FormItemOutputFilename::getName ()]->printHTML ( os );
}
MSToleranceForm::MSToleranceForm ( FormValidatingJavascript* fvj, const VectorConstParameterListPtr& params ) :
	fvj ( fvj )
{
	create ( params );
}
void MSToleranceForm::createItems ()
{
	formItemMap [FormItemMSParentMassTolerance::getName ()] = new FormItemMSParentMassTolerance ();
	formItemMap [FormItemMSParentMassToleranceUnits::getName ()] = new FormItemMSParentMassToleranceUnits ();
	formItemMap [FormItemMassSystematicError::getName ("ms_")]	= new FormItemMassSystematicError (fvj, "ms_");
}
void MSToleranceForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		formItemMap [FormItemMSParentMassTolerance::getName ()]->setValue ( p );
		formItemMap [FormItemMSParentMassToleranceUnits::getName ()]->setValue ( p );
		formItemMap [FormItemMassSystematicError::getName ("ms_")]->setValue ( p );
	}
}
void MSToleranceForm::printHTML ( ostream& os )
{
	formItemMap [FormItemMSParentMassTolerance::getName ()]->printHTML ( os );
	formItemMap [FormItemMSParentMassToleranceUnits::getName ()]->printHTML ( os );
	formItemMap [FormItemMassSystematicError::getName ("ms_")]->printHTML ( os );
	os << "<br />" << endl;
}
MSMSToleranceForm::MSMSToleranceForm ( FormValidatingJavascript* fvj, const VectorConstParameterListPtr& params, bool pChrgItem, bool pChrgRangeItem, int num ) :
	fvj ( fvj ),
	pChrgItem ( pChrgItem ),
	pChrgRangeItem ( pChrgRangeItem ),
	num ( num )
{
	create ( params );
}
void MSMSToleranceForm::createItems ()
{
	if ( pChrgItem ) {
		formItemMap [FormItemMSMSPrecursorCharge::getName ()]		= new FormItemMSMSPrecursorCharge ();
	}
	if ( pChrgRangeItem ) {
		formItemMap [FormItemMSMSPrecursorChargeRange::getName ()]	= new FormItemMSMSPrecursorChargeRange ();
	}
	formItemMap [FormItemParentMassConvert::getName (num)]			= new FormItemParentMassConvert ( "MSMS", num );
	formItemMap [FormItemMSMSParentMassTolerance::getName ()]		= new FormItemMSMSParentMassTolerance ( fvj );
	formItemMap [FormItemMSMSParentMassToleranceUnits::getName ()]	= new FormItemMSMSParentMassToleranceUnits ();
	formItemMap [FormItemFragmentMassesTolerance::getName (num)]		= new FormItemFragmentMassesTolerance ( fvj, "300", num );
	formItemMap [FormItemFragmentMassesToleranceUnits::getName (num)]	= new FormItemFragmentMassesToleranceUnits ( "ppm", num );
	formItemMap [FormItemMassSystematicError::getName ("msms_")]	= new FormItemMassSystematicError (fvj, "msms_");
}
void MSMSToleranceForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		if ( pChrgRangeItem ) {
			formItemMap [FormItemMSMSPrecursorChargeRange::getName ()]->setValue ( p );
		}
		formItemMap [FormItemParentMassConvert::getName (num)]->setValue ( p );
		formItemMap [FormItemMSMSParentMassTolerance::getName ()]->setValue ( p );
		formItemMap [FormItemMSMSParentMassToleranceUnits::getName ()]->setValue ( p );
		formItemMap [FormItemFragmentMassesTolerance::getName (num)]->setValue ( p );
		formItemMap [FormItemFragmentMassesToleranceUnits::getName (num)]->setValue ( p );
		formItemMap [FormItemMassSystematicError::getName ("msms_")]->setValue ( p );
	}
}
void MSMSToleranceForm::printHTML ( ostream& os )
{
	if ( pChrgItem ) {
		formItemMap [FormItemMSMSPrecursorCharge::getName ()]->printHTML ( os );
		os << "<br />" << endl;
	}
	if ( pChrgRangeItem ) {
		formItemMap [FormItemMSMSPrecursorChargeRange::getName ()]->printHTML ( os );
		os << "<br />" << endl;
		os << "(only applicable when unspecified in peak lists)<br />" << endl;
	}
	formItemMap [FormItemParentMassConvert::getName (num)]->printHTML ( os );
	os << "<br />" << endl;
	formItemMap [FormItemMSMSParentMassTolerance::getName ()]->printHTML ( os );
	formItemMap [FormItemMSMSParentMassToleranceUnits::getName ()]->printHTML ( os );
	formItemMap [FormItemMassSystematicError::getName ("msms_")]->printHTML ( os );
	os << "<br />" << endl;
	formItemMap [FormItemFragmentMassesTolerance::getName (num)]->printHTML ( os );
	formItemMap [FormItemFragmentMassesToleranceUnits::getName (num)]->printHTML ( os );
}
FormItemMassSystematicError::FormItemMassSystematicError ( FormValidatingJavascript* fvj, const string& ms ) :
	FormItemText ( "Sys Err", "allman.htm#systematic_error", getName ( ms ), 5, 6, "0" )
{
}
FormItemMSMSPrecursorCharge::FormItemMSMSPrecursorCharge () :
	FormItemSelect ( "Precursor Charge", "allman.htm#precursor_charge", getName (), options, "Automatic" )
{
}
FormItemMSMSPrecursorChargeRange::FormItemMSMSPrecursorChargeRange () :
	FormItemSelect ( "Precursor Charge Range", "allman.htm#precursor_charge_range", getName (), options, "2 3" )
{
}
SearchForm::SearchForm ( const VectorConstParameterListPtr& params, const string& searchName, FormValidatingJavascript* fvj, int num ) :
	fvj ( fvj ),
	searchName ( searchName ),
	searchResultsForm ( params ),
	saveResultsForm ( params, fvj, searchName == "batchtag", searchName == "mshomology" ),
	proteinPIForm ( fvj, params ),
	proteinMWForm ( 0 ),
	msProteinMWForm ( 0 ),
	msmsProteinMWForm ( 0 ),
	num ( num )
{
	pattTypeSearch = searchName == "mspattern" || searchName == "mshomology";
	protSearch = searchName == "dbstat" || pattTypeSearch;
	fitTypeSearch = searchName == "msfit";
	msSearch = fitTypeSearch;
	msmsSearch = searchName == "mstag" || searchName == "msseq" || searchName == "batchtag";
	if ( protSearch )	proteinMWForm = new ProteinMWForm ( fvj, params );
	if ( msSearch )		msProteinMWForm = new ProteinMWForm ( fvj, params, "ms_" );
	if ( msmsSearch )	msmsProteinMWForm = new ProteinMWForm (  fvj, params, "msms_" );
	create ( params );
}
void SearchForm::createItems ()
{
	formItemMap [FormItemMultipleDatabase::getName ()] = new FormItemMultipleDatabase ();
	formItemMap [FormItemUserProteinSequence::getName ()] = new FormItemUserProteinSequence ();

	if ( searchName == "dbstat" )
		formItemMap [FormItemDNAReadingFrame::getName ()] = new FormItemDNAReadingFrame ();
	else
		formItemMap [FormItemDNAFrameTranslation::getName ()] = new FormItemDNAFrameTranslation ();

	formItemMap [FormItemNTerminusAALimit::getName ()] = new FormItemNTerminusAALimit ( fvj );

	formItemMap [FormItemAccessionNums::getName ()]		= new FormItemAccessionNums ();

	if ( pattTypeSearch ) {
		formItemMap [FormItemEnzyme::getName ()] = new FormItemEnzyme ("No enzyme");
	}
	else {
		formItemMap [FormItemEnzyme::getName ()]			= new FormItemEnzyme ("Trypsin", fitTypeSearch);
		if ( !fitTypeSearch && searchName != "dbstat" ) {
			formItemMap [FormItemAllowNonSpecific::getName ()] = new FormItemAllowNonSpecific ();
		}
		formItemMap [FormItemMissedCleavages::getName ()]	= new FormItemMissedCleavages ( fvj, searchName == "dbstat" ? "0" : "1" );
		formItemMap [FormItemConstMod::getName (num)]	= new FormItemConstMod ( true, false, num );
	}
	formItemMap [FormItemTaxonomy::getName ()]		= new FormItemTaxonomy ();
	formItemMap [FormItemTaxonomyRemove::getName ()]= new FormItemTaxonomyRemove ();
	formItemMap [FormItemTaxonomyNames::getName ()]	= new FormItemTaxonomyNames ();
	formItemMap [FormItemNames::getName ()]					= new FormItemNames ();
	formItemMap [FormItemAddAccessionNumbers::getName ()]	= new FormItemAddAccessionNumbers ();
}
void SearchForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemMultipleDatabase::getName ()]->setValue ( p );
		formItemMap [FormItemUserProteinSequence::getName()]->setValue ( p );
		if ( searchName == "dbstat" )
			formItemMap [FormItemDNAReadingFrame::getName ()]->setValue ( p );
		else
			formItemMap [FormItemDNAFrameTranslation::getName ()]->setValue ( p );

		formItemMap [FormItemNTerminusAALimit::getName ()]->setValue ( p );

		formItemMap [FormItemAccessionNums::getName ()]->setValue ( p );

		formItemMap [FormItemEnzyme::getName ()]->setValue ( p );
		if ( !pattTypeSearch ) {
			if ( !fitTypeSearch && searchName != "dbstat" ) {
				formItemMap [FormItemAllowNonSpecific::getName ()]->setValue ( p );
			}
			formItemMap [FormItemMissedCleavages::getName ()]->setValue ( p );
			formItemMap [FormItemConstMod::getName (num)]->setValue ( p );
		}
		if ( p->getStringValue ( FormItemTaxonomy::getName () ) == "NOT HUMAN" ) {
			StringVector sv;
			sv.push_back ( "HOMO SAPIENS" );
			formItemMap [FormItemTaxonomy::getName ()]->setValue ( sv );
			formItemMap [FormItemTaxonomyRemove::getName ()]->setValue ( true );
		}
		else {
			formItemMap [FormItemTaxonomy::getName ()]->setValue ( p );
			formItemMap [FormItemTaxonomyRemove::getName ()]->setValue ( p );
		}
		formItemMap [FormItemTaxonomyNames::getName ()]->setValue ( p );
		formItemMap [FormItemNames::getName ()]->setValue ( p );
		formItemMap [FormItemAddAccessionNumbers::getName ()]->setValue ( p );
	}
}
void SearchForm::printHTML ( ostream& os )
{
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true );
			formItemMap [FormItemMultipleDatabase::getName ()]->printHTML ( os );
			ExpandableJavascriptBlock ejb ( "<a href=\"../html/instruct/allman.htm#user_protein_sequence\">User Protein Sequence</a>", false );
			ejb.printHeader ( os );
			formItemMap [FormItemUserProteinSequence::getName()]->printHTML ( os );
			ejb.printFooter ( os );
			os << "<br />" << endl;
			if ( searchName == "dbstat" )
				formItemMap [FormItemDNAReadingFrame::getName ()]->printHTML ( os );
			else
				formItemMap [FormItemDNAFrameTranslation::getName ()]->printHTML ( os );
			formItemMap [FormItemNTerminusAALimit::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemTaxonomy::getName ()]->printHTML ( os );
			saveResultsForm.printHTML ( os );
		tableHeaderEnd ( os );
		printHTMLTopRight ( os );
	tableRowEnd ( os );
	printHTMLPreSearch ( os );
	showHiddenItems ( os );
}
void SearchForm::printHTMLTopRight ( ostream& os )
{
	tableHeaderStart ( os, "", "left", true );
		formItemMap [FormItemEnzyme::getName ()]->printHTML ( os );
		if ( !pattTypeSearch ) {
			if ( !fitTypeSearch && searchName != "dbstat" ) {
				formItemMap [FormItemAllowNonSpecific::getName ()]->printHTML ( os );
				os << "<br />" << endl;
			}
			formItemMap [FormItemMissedCleavages::getName ()]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemConstMod::getName (num)]->printHTML ( os );
		}
	tableHeaderEnd ( os );
}
void SearchForm::printHTMLPreSearch ( ostream& os )
{
	tableRowStart ( os );
		tableHeaderStart ( os, "", "left", true, 2 );
			ExpandableJavascriptBlock ejb ( "Pre-Search Parameters" );
			ejb.printHeader ( os );
			tableStart ( os );
				tableRowStart ( os );
					tableHeaderStart ( os );
						if ( protSearch ) {
							proteinMWForm->printHTML ( os );
							os << "<br />" << endl;
						}
						if ( msSearch ) {
							msProteinMWForm->printHTML ( os );
							os << "<br />" << endl;
						}
						if ( msmsSearch ) {
							msmsProteinMWForm->printHTML ( os );
							os << "<br />" << endl;
						}
						proteinPIForm.printHTML ( os );
						os << "<br />" << endl;
						searchResultsForm.printHTML ( os );
						os << "<br />" << endl;
						formItemMap [FormItemTaxonomyRemove::getName ()]->printHTML ( os );
					tableHeaderEnd ( os );
					tableHeaderStart ( os );
						tableStart ( os );
							tableRowStart ( os );
								tableHeaderStart ( os );
									formItemMap [FormItemTaxonomyNames::getName ()]->printHTML ( os );
								tableHeaderEnd ( os );
								tableHeaderStart ( os );
									formItemMap [FormItemAccessionNums::getName ()]->printHTML ( os );
								tableHeaderEnd ( os );
							tableRowEnd ( os );
						tableEnd ( os );
						tableStart ( os );
							tableRowStart ( os );
								tableHeaderStart ( os );
									formItemMap [FormItemNames::getName ()]->printHTML ( os );
								tableHeaderEnd ( os );
								tableHeaderStart ( os );
									formItemMap [FormItemAddAccessionNumbers::getName ()]->printHTML ( os );
								tableHeaderEnd ( os );
							tableRowEnd ( os );
						tableEnd ( os );
					tableHeaderEnd ( os );
				tableRowEnd ( os );
			tableEnd ( os );
			ejb.printFooter ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}
void SearchForm::printCGI ( ostream& os ) const
{
	ProspectorForm::printCGI ( os );
	searchResultsForm.printCGI ( os );
	saveResultsForm.printCGI ( os );
	proteinPIForm.printCGI ( os );
	if ( proteinMWForm )	proteinMWForm->printCGI ( os );
	if ( msProteinMWForm )	msProteinMWForm->printCGI ( os );
	if ( msmsProteinMWForm )msmsProteinMWForm->printCGI ( os );
}
void SearchForm::printHTMLJavascriptHidden ( ostream& os ) const
{
	ProspectorForm::printHTMLJavascriptHidden ( os );
	searchResultsForm.printHTMLJavascriptHidden ( os );
	saveResultsForm.printHTMLJavascriptHidden ( os );
	proteinPIForm.printHTMLJavascriptHidden ( os );
	if ( proteinMWForm )	proteinMWForm->printHTMLJavascriptHidden ( os );
	if ( msProteinMWForm )	msProteinMWForm->printHTMLJavascriptHidden ( os );
	if ( msmsProteinMWForm )msmsProteinMWForm->printHTMLJavascriptHidden ( os );
}
FormItemSearchName::FormItemSearchName ( const string& value ) :
	FormItemText ( "", "", getName (), 20, 20, value )
{
}
FormItemReportTitle::FormItemReportTitle ( const string& value ) :
	FormItemText ( "", "", getName (), 20, 20, value )
{
}
FormItemVersion::FormItemVersion () :
	FormItemText ( "", "", getName (), 20, 20, Version::instance ().getVersion () )
{
}
FormItemDatabase::FormItemDatabase () :
	FormItemSelect ( "Database", "allman.htm#database", getName (), SeqdbDir::instance ().getDatabaseList ( false ), getBestSubstituteDatabase ( "SwissProt" ) )
{
}
FormItemDatabase::FormItemDatabase ( bool userDefault ) :
	FormItemSelect ( "Database", "allman.htm#database", getName (), SeqdbDir::instance ().getUserDatabaseList ( false ), userDefault ? "User Protein" : getBestSubstituteDatabase ( "SwissProt" ) )
{
}
void FormItemDatabase::setValue ( const ParameterList* p, const string& n )
{
	select = p->getStringValue ( ( n == "" ) ? _name : n, select );
	if ( select == "User Protein" ) return;
	if ( genFileExists ( SeqdbDir::instance ().getDatabasePath ( select ) ) == false )
		select = getBestSubstituteDatabase ( select );
}
FormItemMultipleDatabase::FormItemMultipleDatabase () :
	FormItemSelectMultiple ( "Database", "allman.htm#database", getName (), SeqdbDir::instance ().getUserDatabaseList ( true ), StringVector (), 4 )
{
	select.push_back ( getBestSubstituteDatabase ( "SwissProt" ) );
}
void FormItemMultipleDatabase::setValue ( const ParameterList* p, const string& n )
{
	StringVector s = p->getStringVectorValue ( ( n == "" ) ? _name : n );
	SetString ss;
	for ( StringVectorSizeType i = 0 ; i < s.size () ; i++ ) {
		const string& si = s [i];
		if ( si == "User Protein" || genFileExists ( SeqdbDir::instance ().getDatabasePath ( si ) ) )
			ss.insert ( si );
		else if ( isSuffix ( si, ".concat" ) ) {
			string file1 = si.substr ( 0, si.length () - 7 );
			string file2;
			if ( isSuffix ( file1, ".random" ) )	file2 = file1.substr ( 0, file1.length () - 7 );
			if ( isSuffix ( file1, ".reverse" ) )	file2 = file1.substr ( 0, file1.length () - 8 );
			bool f1 = genFileExists ( SeqdbDir::instance ().getDatabasePath ( file1 ) );
			bool f2 = genFileExists ( SeqdbDir::instance ().getDatabasePath ( file2 ) );
			if ( f1 && f2 ) ss.insert ( si );
		}
		else
			ss.insert ( getBestSubstituteDatabase ( si ) );
	}
	select.clear ();
	for ( SetStringConstIterator j = ss.begin () ; j != ss.end () ; j++ ) {
		select.push_back ( *j );
	}
}
FormItemEnzyme::FormItemEnzyme ( const string& v, bool noNoEnzyme ) :
	FormItemSelect ( "Digest", "allman.htm#Enzyme", getName (), StringVector (), v )
{
	if ( !noNoEnzyme ) value.push_back ( "No enzyme" );
	StringVector enzymeList = DigestTable::instance ().getNames ();
	copy ( enzymeList.begin (), enzymeList.end (), back_inserter ( value ) );
}
FormItemAllowNonSpecific::FormItemAllowNonSpecific () :
	FormItemSelect ( "Non-Specific", "allman.htm#non_specific", getName (), options, options [0] )
{
}
FormItemMissedCleavages::FormItemMissedCleavages ( FormValidatingJavascript* fvj, const string& val ) :
	FormItemText ( "Max. Missed Cleavages", "allman.htm#Enzyme", getName (), 2, 2, val, fvj->addPositiveIntegerValidator ( getName (), "Max. Missed Cleavages" ) )
{
}
FormItemS::FormItemS ( int num ) :
	FormItemCheckbox ( "", "", getName ( num ), true )
{
}
FormItemNTerm::FormItemNTerm ( int num, bool showLabel ) :
	FormItemSelect ( showLabel ? "N term" : "", "allman.htm#term_groups", getName ( num ), StringVector (), "" )
{
	value = ConstMod::getNTermNames ();
}
FormItemCTerm::FormItemCTerm ( int num, bool showLabel ) :
	FormItemSelect ( showLabel ? "C term" : "", "allman.htm#term_groups", getName ( num ), StringVector (), "" )
{
	value = ConstMod::getCTermNames ();
}
FormItemConstMod::FormItemConstMod ( bool showTerminalMods, bool nonSelected, int num ) :
	FormItemSelectMultiple ( "Constant<br />Mods", "allman.htm#const_mods", getName ( num ), showTerminalMods ? ConstMod::getNames () : ConstMod::getNonTerminalNames (), StringVector (), 6 )
{
	if ( !nonSelected ) {
		select.push_back ( "Carbamidomethyl (C)" );
	}
}
FormItemConstMod::FormItemConstMod () :		// MS-Viewer
	FormItemSelectMultiple ( "Constant Mods<br />(if unspecified<br />in results file)", "viewerman.htm#modifications", getName (), ConstMod::getNames (), StringVector (), 6 )
{
}
FormItemComment::FormItemComment () :
	FormItemText ( "Sample ID (comment)", "allman.htm#comment", getName (), 30, 256, "" )
{
}
FormItemDisplayGraph::FormItemDisplayGraph ( bool set ) :
	FormItemCheckbox ( "Display Graph", "allman.htm#display_graph", getName (), set )
{
}
FormItemDetailedReport::FormItemDetailedReport () :
	FormItemCheckbox ( "Detailed Report", "", getName (), true )
{
}
const char* FormItemModAA::optionsFit [] = {	"Peptide N-terminal Gln to pyroGlu",
												"Oxidation of M",
												"Protein N-terminus Acetylated",
												"Acrylamide Modified Cys",
												"User Defined 1",
												"User Defined 2",
												"User Defined 3",
												"User Defined 4",
												0 };

const char* FormItemModAA::optionsDigest [] = {	"Peptide N-terminal Gln to pyroGlu",
												"Oxidation of M",
												"Protein N-terminus Acetylated",
												"Acrylamide Modified Cys",
												0 };
FormItemModAA::FormItemModAA ( const string& type ) :
	FormItemSelectMultiple ( "Possible<br />Modifications", "allman.htm#mod_AA", getName (), type == "Fit" ? optionsFit : optionsDigest, StringVector (), 4 )
{
	copy ( value.begin (), value.begin () + 3, back_inserter ( select ) );
}
FormItemUserName::FormItemUserName ( const string& num ) :
	FormItemSelect ( string ( "User Def Mod " ) + num, "allman.htm#mod_AA", getName ( num ), Usermod::getNames (), Usermod::getNames () [0] )
{
}
FormItemMaxReportedHits::FormItemMaxReportedHits ( const string& ms, const string& value, FormValidatingJavascript* fvj, bool allowBlank ) :
	FormItemText ( "Maximum Reported Hits", "allman.htm#max_reported_hits", getName ( ms ), 6, 6, value, allowBlank ? fvj->addPositiveNonZeroIntegerAllowBlankValidator ( getName ( ms ), "Maximum Reported Hits" ) : fvj->addPositiveNonZeroIntegerValidator ( getName ( ms ), "Maximum Reported Hits" ) )
{
}
FormItemMinMatches::FormItemMinMatches ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Min. # peptides required to match", "fitman.htm#min_matches", getName (), 2, 2, "4", fvj->addPositiveNonZeroIntegerValidator ( getName (), "Min. # peptides required to match" ) )
{
}
FormItemMowseOn::FormItemMowseOn () :
	FormItemCheckbox ( "Report MOWSE Scores", "fitman.htm#scoring", getName (), true )
{
}
FormItemMowsePfactor::FormItemMowsePfactor ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Pfactor", "", getName (), 3, 3, "0.4", fvj->addPositiveFloatingPointValidator ( getName (), "Pfactor" ) )
{
}
#ifdef CHEM_SCORE
FormItemChemScore::FormItemChemScore () :
	FormItemCheckbox ( "Chem Score", "allman.htm#chem_score", getName (), false )
{
}
FormItemMetOxFactor::FormItemMetOxFactor ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Met Ox Factor", "", getName (), 3, 3, "1.0", fvj->addPositiveFloatingPointValidator (getName (), "Met Ox Factor") )
{
}
#endif
FormItemMaxHits::FormItemMaxHits ( const string& val ) :
	FormItemText ( "Maximum Hits", "allman.htm#max_hits", getName (), 8, 8, val )
{
}
FormItemMSMSMaxModifications::FormItemMSMSMaxModifications ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Max Mods", "allman.htm#max_mods", getName (), 2, 2, "2", fvj->addPositiveIntegerValidator ( getName (), "Max Mods" ) + "; setMaxOptions (this.form)" )
{
}
FormItemMSMSMaxPeptidePermutations::FormItemMSMSMaxPeptidePermutations ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Max Peptide Permutations", "allman.htm#max_peptide_permutations", getName (), 7, 7, "", fvj->addPositiveNonZeroIntegerAllowBlankValidator ( getName (), "Max Peptide Permutations" ) )
{
}
FormItemModAAList::FormItemModAAList ( int siz ) :
	FormItemSelectMultiple ( "Variable<br />Mods", "allman.htm#variable_mods", getName ( "msms_" ), StringVector (), StringVector (), siz )
{
}
FormItemModAAList::FormItemModAAList ( bool msms, bool digest, int siz ) :
	FormItemSelectMultiple ( "Variable<br />Mods", "allman.htm#variable_mods", getName ( msms ? "msms_" : "" ), Usermod::getNames (), StringVector (), siz )
{
	if ( msms || digest ) {
		select.push_back ( "Acetyl (Protein N-term)" );
		select.push_back ( "Gln->pyro-Glu (N-term Q)" );
		select.push_back ( "Met-loss (Protein N-term M)" );
		select.push_back ( "Met-loss+Acetyl (Protein N-term M)" );
		select.push_back ( "Oxidation (M)" );
		if ( msms ) {
			select.push_back ( "Acetyl+Oxidation (Protein N-term M)" );
		}
	}
}
FormItemParentMassConvert::FormItemParentMassConvert ( const string& type, int num ) :
	FormItemSelect ( "Masses are", "allman.htm#mass_type", getName ( num ), type == "MS" ? optionsMS : optionsMSMS, "monoisotopic" )
{
}
FormItemMSParentMassTolerance::FormItemMSParentMassTolerance () :
	FormItemText ( "Tol", "allman.htm#mass_tolerance", getName (), 5, 6, "20" )
{
}
FormItemMSParentMassToleranceUnits::FormItemMSParentMassToleranceUnits ( const string& v ) :
	FormItemSelect ( "", "allman.htm#mass_tolerance", getName (), options, v )
{
}
FormItemMSMSParentMassTolerance::FormItemMSMSParentMassTolerance ( FormValidatingJavascript* fvj, const string& value ) :
	FormItemText ( "Parent Tol", "allman.htm#mass_tolerance", getName (), 4, 10, value, fvj->addPositiveFloatingPointOrExponentValidator ( getName (), "Parent Tol" ) )
{
}
FormItemMSMSParentMassToleranceUnits::FormItemMSMSParentMassToleranceUnits ( const string& v ) :
	FormItemSelect ( "", "allman.htm#mass_tolerance", getName (), options, v )
{
}
FormItemFragmentMassesTolerance::FormItemFragmentMassesTolerance ( FormValidatingJavascript* fvj, const string& v, int num ) :
	FormItemText ( "Frag Tol", "allman.htm#mass_tolerance", getName ( num ), 5, 10, v, fvj->addPositiveFloatingPointOrExponentValidator ( getName ( num ), "Frag Tol" ) )
{
}
FormItemFragmentMassesToleranceUnits::FormItemFragmentMassesToleranceUnits ( const string& v, int num ) :
	FormItemSelect ( "", "allman.htm#mass_tolerance", getName ( num ), options, v )
{
}
FormItemMaxSavedTagHits::FormItemMaxSavedTagHits ( FormValidatingJavascript* fvj, const string& value ) :
	FormItemText ( "# Saved Tag Hits", "", getName (), 5, 8, value, fvj->addPositiveNonZeroIntegerValidator ( getName (), "# Saved Tag Hits" ) )
{
}
FormItemUseInstrumentIonTypes::FormItemUseInstrumentIonTypes ( bool val, int num ) :
	FormItemCheckbox ( "Use instrument specific defaults to override ion types below", "", getName ( num ), val )
{
}
FormItemCreateParams::FormItemCreateParams () :
	FormItemCheckbox ( "Create Parameter File", "", getName (), false )
{
}
FormItemScriptFilename::FormItemScriptFilename () :
	FormItemText ( "Filename", "", getName (), 20, 50, "script" )
{
}
FormItemParentContaminantMasses::FormItemParentContaminantMasses () :
	FormItemTextArea ( "Contaminant<br />Masses", "allman.htm#contaminant_masses", getName (), 4, 20, StringVector () )
{
}
FormItemCompMaskType::FormItemCompMaskType ( const string& defaultLogic ) :
	FormItemSelect ( "", "allman.htm#present_amino_acids", getName (), options, defaultLogic )
{
}
FormItemUserProteinSequence::FormItemUserProteinSequence () :
	FormItemTextArea ( "", "", getName (), 6, 70, StringVector () )
{
}
FormItemUserProteinSequence::FormItemUserProteinSequence ( bool databaseSequence, const string& sequence ) :
	FormItemTextArea ( "User Protein<br />Sequence", "allman.htm#user_protein_digest", getName (), 12, 80, StringVector () )
{
	static char defaultInstructions [] =
		"Paste or type sequence in this field. Delete these instructions first.\n\n"
		"Tabs, returns, and spaces are ignored.\n\n"
		"USE CAPITAL LETTERS. The following lower case letters can be used:\n"
		"  U       - Selenocysteine\n"
		"  s,t,y   - Phosphorylated S,T,Y\n"
		"  u,v,w,x - user specified amino acids (provide elemental composition --->)\n\n"
		"Do NOT use the letters B, J, O, X, or Z.\n\n"
		"Do NOT use 3-letter code amino acid symbols.\n";
	static char defaultInstructions2 [] =
		"Paste or type sequence in this field. Delete these instructions first.\n"
		"Tabs, returns, and spaces are ignored.\n"
		"USE CAPITAL LETTERS.\n"
		"Do NOT use the letter O.\n"
		"Use B for D or N.\n"
		"Use U for Selenocysteine.\n"
		"Use X for Unknown Amino Acid.\n"
		"Use Z for E or Q.\n"
		"Do NOT use 3-letter code amino acid symbols.\n"
		"DNA sequences can be added to DNA databases.\n";

	if ( sequence.empty () ) {
		if ( databaseSequence )
			value.push_back ( defaultInstructions2 );
		else
			value.push_back ( defaultInstructions );
	}
	else
		value.push_back ( sequence );
}
FormItemIonType::FormItemIonType ( const string& label, const string& value, bool checked ) :
	FormItemCheckbox ( label, "", getName (), checked, value )
{
}
FormItemAlternative::FormItemAlternative () :
	FormItemCheckbox ( "Alt", "", getName (), false )
{
}
FormItemDiscriminating::FormItemDiscriminating () :
	FormItemCheckbox ( "Discriminating", "", getName (), false )
{
}
const char* FormItemMaxLosses::options [] = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "All", 0 };
FormItemMaxLosses::FormItemMaxLosses () :
	FormItemSelect ( "Max Losses", "prodman.htm#Max_losses", getName (), options, "1" )
{
}
FormItemMultiZInternal::FormItemMultiZInternal () :
	FormItemCheckbox ( "Multi Z Internal", "prodman.htm#Multi_z_internal", getName (), false )
{
}
FormItemCalibrate::FormItemCalibrate () :
	FormItemCheckbox ( "Cal", "prodman.htm#Cal", getName (), false )
{
}
FormItemCalTolerance::FormItemCalTolerance ( const string& units, const string& val ) :
	FormItemText ( "Cal Tol (" + units + ")", "prodman.htm#Cal", getName (), 4, 10, val )
{
}
FormItemDataPlotted::FormItemDataPlotted ( const string& method ) :
	FormItemSelect ( "", "", getName (), options, method )
{
}
FormItemHideProteinSequence::FormItemHideProteinSequence ( bool flag ) :
	FormItemCheckbox ( "Hide Protein Sequences", "allman.htm#hide_protein_sequence", getName (), flag )
{
}
FormItemExpectationCalculationMethod::FormItemExpectationCalculationMethod ( const string& method ) :
	FormItemSelect ( "Expectation Calc Method", "", getName (), options, method )
{
}
FormItemMSMSPkFilter::FormItemMSMSPkFilter ( bool allowUnprocessed, int num ) :
	FormItemSelect ( "", "", getName ( num ), allowUnprocessed ? options : options2, options [0] )
{
}
FormItemMaxMSMSPeaks::FormItemMaxMSMSPeaks ( FormValidatingJavascript* fvj, bool allowBlank, int num ) :
	FormItemText ( "", "", getName ( num ), 4, 6, "", allowBlank ? fvj->addPositiveNonZeroIntegerAllowBlankValidator ( getName ( num ), "Max MSMS Peaks" ) : fvj->addPositiveNonZeroIntegerValidator ( getName ( num ), "Max MSMS Peaks" ) )
{
}
FormItemMSMSMinPrecursorMass::FormItemMSMSMinPrecursorMass ( FormValidatingJavascript* fvj, bool allowBlank ) :
	FormItemText ( "", "", getName (), 4, 6, "200", allowBlank ? fvj->addPositiveNonZeroIntegerAllowBlankValidator ( getName (), "Max MSMS Peaks" ) : fvj->addPositiveNonZeroIntegerValidator ( getName (), "Max MSMS Peaks" ) )
{
}
FormItemRawType::FormItemRawType ( const string& rawType ) :
	FormItemSelect ( "Raw Type", "", getName (), options, rawType )
{
}
FormItemIsotopePurity::FormItemIsotopePurity ( const string& element, int isotope, FormValidatingJavascript* fvj, double value ) :
	FormItemText ( "<sup>" + gen_itoa ( isotope ) + "</sup>" + element + "%", "allman.htm#percent_element", getName ( element, isotope ), 4, 4, value == 100.0 ? "100" : gen_ftoa ( value, "%.1f" ), fvj->addPositiveFloatingPointValidator ( getName ( element, isotope ), gen_itoa ( isotope ) + element + "%" ) )
{
}
FormItemReporterIonWindow::FormItemReporterIonWindow ( FormValidatingJavascript* fvj, const string& value ) :
	FormItemText ( "Reporter Ion Window (Da)", "", getName (), 6, 10, value, fvj->addPositiveFloatingPointValidator ( getName (), "Reporter Ion Window" ) )
{
}
FormItemLinkSearchType::FormItemLinkSearchType ( const string& val, const string& chargeFunction, int num ) :
	FormItemSelect ( "Link Search Type", "bridgeman.htm#link_search_type", getName ( num ), LinkInfo::getNameList (), val.empty () ? "Disulfide (C)" : val, 1, chargeFunction )
{
}
FormItemMaxLinkMolecules::FormItemMaxLinkMolecules ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Maximum Link Molecules", "bridgeman.htm#maximum_link_molecules", getName (), 2, 2, "2", fvj->addPositiveIntegerValidator ( getName (), "Maximum Link Molecules" ) )
{
}

class LinkAAOptions {
	StringVector linkOptions;
	StringVector modOptions;
	string linkChoice;
	LinkAAOptions ();
public:
	static LinkAAOptions& instance ();
	StringVector getLinkOptions () { return linkOptions; }
	StringVector getModOptions () { return modOptions; }
	string getLinkChoice () { return linkChoice; }
};

LinkAAOptions& LinkAAOptions::instance ()
{
	static LinkAAOptions d;
	return d;
}
LinkAAOptions::LinkAAOptions ()
{
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( "link_aa.txt" ) );
	string line;
	while ( ifs.getUncommentedLine ( line ) ) {
		linkOptions.push_back ( line );
	}
	linkChoice = linkOptions [0];
	SetString setAA;
	for ( int i = 0 ; i < linkOptions.size () ; i++ ) {
		StringVector sv = genGetSeparatedValues ( linkOptions [i], "->" );
		for ( int j = 0 ; j < sv.size () ; j++ ) {
			StringVector sv2 = genGetSeparatedValues ( sv [j], "," );
			for ( int k = 0 ; k < sv2.size () ; k++ ) {
				setAA.insert ( sv2 [k] );
			}
		}
	}
	copy ( setAA.begin (), setAA.end (), back_inserter ( modOptions ) );
}

FormItemLinkAA::FormItemLinkAA ( const string& chargeFunction, int num ) :
	FormItemSelect ( "Link AAs", "", getName (num), LinkAAOptions::instance ().getLinkOptions (), LinkAAOptions::instance ().getLinkChoice (), 1, chargeFunction )
{
}
FormItemComposition::FormItemComposition ( const string& label, const string& n, int num, const string& v, FormValidatingJavascript* fvj ) :
	FormItemText ( label + string ( " Elem Comp" ), "", getName ( n, num ), 20, 30, v, fvj->addElementalModificationFormulaAllowBlankValidator ( getName ( n, num ), label + string ( " Elem Comp" ) ) )
{
}
FormItemAccurateMass::FormItemAccurateMass ( const string& n, int num, FormValidatingJavascript* fvj ) :
	FormItemText ( "or Accurate Mass", "", getName ( n, num ), 10, 10, "", fvj->addSignedFloatingPointAllowBlankValidator ( getName ( n ), string ( "Accurate Mass" ) ) )
{
}
FormItemLabel::FormItemLabel ( const string& label, const string& n, int num, const string& v ) :
	FormItemText ( label + string ( " Label" ), "", getName ( n, num ), 20, 30, v )
{
}
const char* FormItemSpecificity::modOptions [] = {	"",
													"A","C","D",
													"DE",
													"E","F","G","H","I","K",
													"Uncleaved K","Uncleaved KR",
													"L","M","N","P","Q","R",
													"Uncleaved R",
													"S",
													"ST",
													"STY",
													"T","V","W","Y",
													"N-term",
													"N-term C",
													"N-term E",
													"N-term G",
													"N-term Q",
													"Protein N-term",
													"Protein N-term M",
													"C-term",
													"C-term M",
													"Neutral loss",
													0 };


FormItemSpecificity::FormItemSpecificity ( const string& n, int num ) :
	FormItemSelect ( "Specificity", "", getName ( n, num ), modOptions, modOptions [0] )
{
}
const char* FormItemModFile::modFileOptions [] = {	"All", "Protein", "Peptide", "Site", 0 };

FormItemModFile::FormItemModFile () :
	FormItemSelect ( "", "", getName (), modFileOptions, modFileOptions [0] )
{
}
const char* FormItemMotifOffset::offsetOptions [] = { "Off", "0", "-1", "+1", "-2", "+2", "-3", "+3", 0 };

FormItemMotifOffset::FormItemMotifOffset ( const string& n, int num ) :
	FormItemSelect ( "Motif", "", getName ( n, num ), offsetOptions, offsetOptions [0] )
{
}
FormItemMotif::FormItemMotif ( FormValidatingJavascript* fvj, const string& n, int num ) :
	FormItemText ( "", "", getName ( n, num ), 10, 20, "N[^P][ST]", fvj->addMSPatternValidator ( getName ( n, num ), "Motif" ) )
{
}
const char* FormItemLimit::limitOptions [] = { "Common", "Rare", "Label 1", "Label 2", "Label 3", 0 };

FormItemLimit::FormItemLimit ( const string& n, int num ) :
	FormItemSelect ( "", "", getName ( n, num ), limitOptions, limitOptions [0] )
{
}
void FormItemLimit::setOptions ( const ParameterList* p, const string& n )
{
	StringVector sv;
	sv.push_back ( "Common" );
	sv.push_back ( "Rare" );
	sv.push_back ( "Label 1" );
	sv.push_back ( "Label 2" );
	sv.push_back ( "Label 3" );
	int maxMods = p->getIntValue ( n );
	for ( int i = 1 ; i < maxMods ; i++ ) {
		sv.push_back ( "Max " + gen_itoa (i) );
	}
	FormItemSelect::setOptions ( sv );
}
FormItemAAModified::FormItemAAModified ( const string& n, const string& v ) :
	FormItemSelect ( "", "", getName ( n ), LinkAAOptions::instance ().getModOptions (), v )
{
}
class FormItemUserAAComposition : public FormItemText {
public:
	FormItemUserAAComposition ( const string& aa, const string& num, const string& value, FormValidatingJavascript* fvj ) :
		FormItemText ( string ("User Specified AA Elem Comp  (") + aa + string (")"), "allman.htm#user_AA", getName (num), 30, 30, value, fvj->addElementalFormulaAllowBlankValidator (getName (num), string ("User Specified AA Elem Comp  (") + aa + string (")")) ) {}
	static string getName ( const string& num )
	{
		if ( num != "" )
			return string ( "user_aa_" ) + num + string ( "_composition" );
		else
			return string ( "user_aa_composition" );
	}
};
UserAAForm::UserAAForm ( const VectorConstParameterListPtr& params, FormValidatingJavascript* fvj ) :
	fvj ( fvj )
{
	create ( params );
}
void UserAAForm::createItems ()
{
	formItemMap [FormItemUserAAComposition::getName ("")]	= new FormItemUserAAComposition ( "u", "", "C2 H3 N1 O1", fvj );
	formItemMap [FormItemUserAAComposition::getName ("2")]	= new FormItemUserAAComposition ( "v", "2", "", fvj );
	formItemMap [FormItemUserAAComposition::getName ("3")]	= new FormItemUserAAComposition ( "w", "3", "", fvj );
	formItemMap [FormItemUserAAComposition::getName ("4")]	= new FormItemUserAAComposition ( "x", "4", "", fvj );
}
void UserAAForm::printHTML ( ostream& os )
{
	formItemMap [FormItemUserAAComposition::getName("")]->printHTML ( os );
	os << "<br />" << endl;
	formItemMap [FormItemUserAAComposition::getName("2")]->printHTML ( os );
	os << "<br />" << endl;
	formItemMap [FormItemUserAAComposition::getName("3")]->printHTML ( os );
	os << "<br />" << endl;
	formItemMap [FormItemUserAAComposition::getName("4")]->printHTML ( os );
	os << "<br />" << endl;
}
void UserAAForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemUserAAComposition::getName("")]->setValue ( p );
		formItemMap [FormItemUserAAComposition::getName("2")]->setValue ( p );
		formItemMap [FormItemUserAAComposition::getName("3")]->setValue ( p );
		formItemMap [FormItemUserAAComposition::getName("4")]->setValue ( p );
	}
}

class FormItemCompIon : public FormItemCheckbox {
public:
	FormItemCompIon ( const string& value, const string& name ) :
		FormItemCheckbox ( "", "", name, false, value ) {}
	FormItemCompIon ( const string& label, const string& value, bool checked, const string& name ) :
		FormItemCheckbox ( label, "", name, checked, value ) {}
	static string getCompIonName () { return "comp_ion"; }
	static string getAAExcludeName () { return "aa_exclude"; }
	static string getAAAddName () { return "aa_add"; }
};

CompIonForm::CompIonForm ( const VectorConstParameterListPtr& params, const string& aaList ) :
	aaList ( aaList )
{
}
void CompIonForm::createItems ()
{
	for ( StringSizeType i = 0 ; i < aaList.size () ; i++ ) {
		formItemMap [getName () + gen_itoa ( i )]	= new FormItemCompIon ( aaList.substr ( i, 1 ), getName () );
	}
}
void CompIonForm::setValues ( const VectorConstParameterListPtr& params )
{
}
PresentCompIonForm::PresentCompIonForm ( const VectorConstParameterListPtr& params, const string& aaList, const string& prefix ) :
	CompIonForm ( params, aaList ),
	prefix ( prefix )
{
	create ( params );
}
void PresentCompIonForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		for ( StringSizeType i = 0 ; i < aaList.size () ; i++ ) {
			formItemMap [getName () + gen_itoa ( i )]->setValue ( p );
		}
	}
}
string PresentCompIonForm::getName () const { return prefix + FormItemCompIon::getCompIonName (); }
void PresentCompIonForm::printHTML ( ostream& os )
{
	tableStart ( os, false, "", "0" );
		tableRowStart ( os );
			for ( StringSizeType i = 0 ; i < aaList.size () ; i++ ) {
				tableHeaderStart ( os, "", "center", true );
					string num = gen_itoa ( i );
					formItemMap [getName () + gen_itoa ( i )]->printHTML ( os );
					os << "<br />" << aaList [i] << endl;
				tableHeaderEnd ( os );
			}
		tableRowEnd ( os );
	tableEnd ( os );
}
void PresentCompIonForm::printHTML2Rows ( ostream& os )
{
	int numCols = ( aaList.size () / 2 ) + 1;
	tableStart ( os, false, "", "0" );
		tableRowStart ( os );
			for ( StringSizeType i = 0 ; i < numCols ; i++ ) {
				tableHeaderStart ( os, "", "center", true );
					string num = gen_itoa ( i );
					formItemMap [getName () + gen_itoa ( i )]->printHTML ( os );
					os << "<br />" << aaList [i] << endl;
				tableHeaderEnd ( os );
			}
		tableRowEnd ( os );
		tableRowStart ( os );
			for ( StringSizeType j = numCols ; j < aaList.size () ; j++ ) {
				tableHeaderStart ( os, "", "center", true );
					string num = gen_itoa ( j );
					formItemMap [getName () + gen_itoa ( j )]->printHTML ( os );
					os << "<br />" << aaList [j] << endl;
				tableHeaderEnd ( os );
			}
		tableRowEnd ( os );
	tableEnd ( os );
}
MSCompExtraIonForm::MSCompExtraIonForm ( const VectorConstParameterListPtr& params ) :
	CompIonForm ( params, "ACDEFGHIKLMNPQRSTVUWYmhstyuvwx" )
{
	create ( params );
}
string MSCompExtraIonForm::getName () const { return FormItemCompIon::getCompIonName (); }
void MSCompExtraIonForm::printHTML ( ostream& os )
{
	tableStart ( os, true, "", "0" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, aaList.size () );
				os << "Extra Known Amino Acids (Check box for each additional amino acid known to be present)" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		for ( int i = 0 ; i < 3 ; i++ ) {
			tableRowStart ( os );
				for ( StringSizeType j = 0 ; j < aaList.size () ; j++ ) {
					tableDataStart ( os, "", "center", true, 1 );
						string num = gen_itoa ( j );
						formItemMap [getName () + gen_itoa ( j )]->printHTML ( os );
						os << "<br />" << aaList [j] << endl;
					tableDataEnd ( os );
				}
			tableRowEnd ( os );
		}
	tableEnd ( os );
}
MSCompAbsentIonForm::MSCompAbsentIonForm ( const VectorConstParameterListPtr& params ) :
	CompIonForm ( params, "ACDEFGHIKLMNPQRSTUVWY" )
{
	create ( params );
}
string MSCompAbsentIonForm::getName () const { return FormItemCompIon::getAAExcludeName (); }
void MSCompAbsentIonForm::printHTML ( ostream& os )
{
	tableStart ( os, true, "", "0" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, aaList.size () );
				os << "Absent Amino Acids (Check to Prevent Possible Inclusion)" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			for ( StringSizeType j = 0 ; j < aaList.size () ; j++ ) {
				tableDataStart ( os, "", "center", true, 1 );
					string num = gen_itoa ( j );
					formItemMap [getName () + gen_itoa ( j )]->printHTML ( os );
					os << "<br />" << aaList [j] << endl;
				tableDataEnd ( os );
			}
		tableRowEnd ( os );
	tableEnd ( os );
}
MSCompModifiedIonForm::MSCompModifiedIonForm ( const VectorConstParameterListPtr& params ) :
	CompIonForm ( params, "mhstyuvwx" )
{
	create ( params );
}
string MSCompModifiedIonForm::getName () const { return FormItemCompIon::getAAAddName (); }
void MSCompModifiedIonForm::printHTML ( ostream& os )
{
	static const char* labels [] = { "HTML", "XML", 0 };
	tableStart ( os, true, "", "0" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, aaList.size () );
				os << "Modified Amino Acids Possibly Present (Check to Allow Inclusion)" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			for ( StringSizeType i = 0 ; i < aaList.size () ; i++ ) {
				tableDataStart ( os, "", "center", true, 1 );
					string num = gen_itoa ( i );
					formItemMap [getName () + gen_itoa ( i )]->printHTML ( os );
					os << "<br />" << aaList [i] << endl;
				tableDataEnd ( os );
			}
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				os << "Met-ox" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 1 );
				os << "Homoserine<br />lactone" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 3 );
				os << "Phosphorylated<br />S,T,Y" << endl;
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true, 4 );
				os << "User Specified AAs<br />(Provide Elemental Composition Below)" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
}
ImmoniumIonForm::ImmoniumIonForm ( const VectorConstParameterListPtr& params )
{
	create ( params );
}
string ImmoniumIonForm::getName () const { return FormItemCompIon::getCompIonName (); }
void ImmoniumIonForm::createItems ()
{
	formItemMap [getName()+string("60")]	= new FormItemCompIon ( "60", "S", false, getName () );
	formItemMap [getName()+string("70")]	= new FormItemCompIon ( "70", "[RP]", false, getName () );
	formItemMap [getName()+string("72")]	= new FormItemCompIon ( "72", "V", false, getName () );
	formItemMap [getName()+string("74")]	= new FormItemCompIon ( "74", "T", false, getName () );
	formItemMap [getName()+string("86")]	= new FormItemCompIon ( "86", "[LI]", false, getName () );
	formItemMap [getName()+string("87")]	= new FormItemCompIon ( "87", "[NR]", false, getName () );
	formItemMap [getName()+string("88")]	= new FormItemCompIon ( "88", "D", false, getName () );
	formItemMap [getName()+string("129")]	= new FormItemCompIon ( "<b>84,101</b>,129", "[KQ]", false, getName () );
	formItemMap [getName()+string("102")]	= new FormItemCompIon ( "102", "E", false, getName () );

	formItemMap [getName()+string("104")]	= new FormItemCompIon ( "104", "M", false, getName () );
	formItemMap [getName()+string("110")]	= new FormItemCompIon ( "110", "H", false, getName () );
	formItemMap [getName()+string("120")]	= new FormItemCompIon ( "120", "F", false, getName () );
	formItemMap [getName()+string("126")]	= new FormItemCompIon ( "126", "P", false, getName () );
	formItemMap [getName()+string("136")]	= new FormItemCompIon ( "136", "Y", false, getName () );
	formItemMap [getName()+string("170")]	= new FormItemCompIon ( "117,130,<b>159</b>,170", "W", false, getName () );
	formItemMap [getName()+string("185")]	= new FormItemCompIon ( "<b>70</b>,73,87,100,<b>112</b>,185", "R", false, getName () );

	formItemMap [getName()+string("mR")]	= new FormItemCompIon ( "-SRH,-R", "C", false, getName () );
	formItemMap [getName()+string("m80")]	= new FormItemCompIon ( "-80", "[sty]", false, getName () );
	formItemMap [getName()+string("m98")]	= new FormItemCompIon ( "-98", "[st]", false, getName () );
	formItemMap [getName()+string("m64")]	= new FormItemCompIon ( "-64", "m", false, getName () );
	formItemMap [getName()+string("m17")]	= new FormItemCompIon ( "-17", "R", false, getName () );
}
void ImmoniumIonForm::setMSCompImmoniumDefaults ()
{
	formItemMap [getName()+string("72")]->setValue ( true );
	formItemMap [getName()+string("102")]->setValue ( true );
	formItemMap [getName()+string("110")]->setValue ( true );
	formItemMap [getName()+string("136")]->setValue ( true );
	formItemMap [getName()+string("170")]->setValue ( true );
	formItemMap [getName()+string("185")]->setValue ( true );
}
void ImmoniumIonForm::printHTML ( ostream& os )
{
	tableStart ( os, true, "", "0" );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 10 );
				os << "AA Composition (based on immonium and related ions)" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("60")]->printHTML ( os );
				os << "<br />S" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("70")]->printHTML ( os );
				os << "<br />R,P" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("72")]->printHTML ( os );
				os << "<br />V" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("74")]->printHTML ( os );
				os << "<br />T" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("86")]->printHTML ( os );
				os << "<br />L,I" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("87")]->printHTML ( os );
				os << "<br />N,R" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("88")]->printHTML ( os );
				os << "<br />D" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 2 );
				formItemMap [getName()+string("129")]->printHTML ( os );
				os << "<br />K,Q" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("102")]->printHTML ( os );
				os << "<br />E" << endl;
			tableDataEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("104")]->printHTML ( os );
				os << "<br />M" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("110")]->printHTML ( os );
				os << "<br />H" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("120")]->printHTML ( os );
				os << "<br />F" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("126")]->printHTML ( os );
				os << "<br />P" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 1 );
				formItemMap [getName()+string("136")]->printHTML ( os );
				os << "<br />Y" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 2 );
				formItemMap [getName()+string("170")]->printHTML ( os );
				os << "<br />W" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 3 );
				formItemMap [getName()+string("185")]->printHTML ( os );
				os << "<br />R" << endl;
			tableDataEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 10 );
				os << "AA Composition (based on loss from parent ion)" << endl;
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableDataStart ( os, "", "center", true, 2 );
				formItemMap [getName()+string("mR")]->printHTML ( os );
				os << "<br />C<br />(select&nbsp;modified&nbsp;Cys)" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 2 );
				formItemMap [getName()+string("m80")]->printHTML ( os );
				os << "<br />PO<sub>4</sub>&nbsp;-&nbsp;S,T,Y<br />(homology)" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 2 );
				formItemMap [getName()+string("m98")]->printHTML ( os );
				os << "<br />PO<sub>4</sub>&nbsp;-&nbsp;S,T<br />(homology)" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 2 );
				formItemMap [getName()+string("m64")]->printHTML ( os );
				os << "<br />Mso<br />(homology)" << endl;
			tableDataEnd ( os );
			tableDataStart ( os, "", "center", true, 2 );
				formItemMap [getName()+string("m17")]->printHTML ( os );
				os << "<br />R<br />&nbsp;" << endl;
			tableDataEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableDataStart ( os, "", "left", true, 10 );
				os << "Note: For ions with multiple residues (i.e. 86: means I or L)." << endl;
				os << "<br />" << endl;
				os << "Note: Requires presence of at least one AA for each composition ion selected." << endl;
			tableDataEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
}
void ImmoniumIonForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [getName()+string("60")]->setValue ( p );
		formItemMap [getName()+string("70")]->setValue ( p );
		formItemMap [getName()+string("72")]->setValue ( p );
		formItemMap [getName()+string("74")]->setValue ( p );
		formItemMap [getName()+string("86")]->setValue ( p );
		formItemMap [getName()+string("87")]->setValue ( p );
		formItemMap [getName()+string("88")]->setValue ( p );
		formItemMap [getName()+string("129")]->setValue ( p );
		formItemMap [getName()+string("102")]->setValue ( p );

		formItemMap [getName()+string("104")]->setValue ( p );
		formItemMap [getName()+string("110")]->setValue ( p );
		formItemMap [getName()+string("120")]->setValue ( p );
		formItemMap [getName()+string("126")]->setValue ( p );
		formItemMap [getName()+string("136")]->setValue ( p );
		formItemMap [getName()+string("170")]->setValue ( p );
		formItemMap [getName()+string("185")]->setValue ( p );

		formItemMap [getName()+string("mR")]->setValue ( p );
		formItemMap [getName()+string("m80")]->setValue ( p );
		formItemMap [getName()+string("m98")]->setValue ( p );
		formItemMap [getName()+string("m64")]->setValue ( p );
		//formItemMap [getName()+string("m17")]->setValue ( p );	// R in the default file should set "185" not this one
	}
}
HeaderForm::HeaderForm ( const VectorConstParameterListPtr& params, const string& searchName, const string& reportTitle ) :
	searchName ( searchName ),
	reportTitle ( reportTitle )
{
	create ( params );
}
void HeaderForm::createItems ()
{
	formItemMap [FormItemSearchName::getName ()]	= new FormItemSearchName ( searchName );
	formItemMap [FormItemReportTitle::getName ()]	= new FormItemReportTitle ( reportTitle );
	formItemMap [FormItemVersion::getName ()]	= new FormItemVersion ();
}
void HeaderForm::printHTML ( ostream& os )
{
	showHiddenItems ( os );
}

class FormItemModRangeType : public FormItemSelect {
	static const char* options [];
public:
	FormItemModRangeType () :
		FormItemSelect ( "", "", getName (), options, "Da", 1, "setModRangeTypeDefaults( this.form )" ) {}
	static string getName () { return "mod_range_type"; }
};
const char* FormItemModRangeType::options [] = { "Da", "m/z", 0 };

class FormItemModStartNominal : public FormItemText {
public:
	FormItemModStartNominal () :
		FormItemText ( "", "", getName (), 5, 5, "-100" ) {}
	static string getName () { return "mod_start_nominal"; }
};

class FormItemModEndNominal : public FormItemText {
public:
	FormItemModEndNominal () :
		FormItemText ( "", "", getName (), 5, 5, "100" ) {}
	static string getName () { return "mod_end_nominal"; }
};

class FormItemModDefect : public FormItemText {
public:
	FormItemModDefect () :
		FormItemText ( "Defect", "", getName (), 7, 7, "0.00048" ) {}
	static string getName () { return "mod_defect"; }
};

class FormItemModMaxZ : public FormItemSelect {
	static const char* options [];
public:
	FormItemModMaxZ () :
		FormItemSelect ( "Max Z", "", getName (), options, "4" ) {}
	static string getName () { return "mod_max_z"; }
};
const char* FormItemModMaxZ::options [] = { "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 0 };

class FormItemModUncleaved : public FormItemCheckbox {
public:
	FormItemModUncleaved () :
		FormItemCheckbox ( "Uncleaved", "", getName (), false ) {}
	static string getName () { return "mod_uncleaved"; }
};

class FormItemModNTermType : public FormItemSelect {
	static const char* options [];
public:
	FormItemModNTermType () :
		FormItemSelect ( "", "", getName (), options, "Peptide" ) {}
	static string getName () { return "mod_n_term_type"; }
};
const char* FormItemModNTermType::options [] = { "Peptide", "Protein", 0 };

class FormItemModNTerm : public FormItemCheckbox {
public:
	FormItemModNTerm () :
		FormItemCheckbox ( "N Term", "", getName (), false ) {}
	static string getName () { return "mod_n_term"; }
};

class FormItemModCTermType : public FormItemSelect {
	static const char* options [];
public:
	FormItemModCTermType () :
		FormItemSelect ( "", "", getName (), options, "Peptide" ) {}
	static string getName () { return "mod_c_term_type"; }
};
const char* FormItemModCTermType::options [] = { "Peptide", "Protein", 0 };

class FormItemModCTerm : public FormItemCheckbox {
public:
	FormItemModCTerm () :
		FormItemCheckbox ( "C Term", "", getName (), false ) {}
	static string getName () { return "mod_c_term"; }
};

class FormItemModNeutralLoss : public FormItemCheckbox {
public:
	FormItemModNeutralLoss () :
		FormItemCheckbox ( "Neutral Loss", "", getName (), false ) {}
	static string getName () { return "mod_neutral_loss"; }
};

class FormItemModRare : public FormItemCheckbox {
public:
	FormItemModRare () :
		FormItemCheckbox ( "Rare", "", getName (), false ) {}
	static string getName () { return "mod_rare"; }
};

class FormItemModMotifAA : public FormItemSelect {
	static const char* options [];
public:
	FormItemModMotifAA ( int num, const string& aa ) :
		FormItemSelect ( "Motif " + gen_itoa ( num ), "", getName (num), options, aa ) {}
	static string getName ( int num ) { return "mod_motif_aa_" + gen_itoa ( num ); }
};
const char* FormItemModMotifAA::options [] = { "A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y",0 };

class FormItemModMotifOffset : public FormItemSelect {
	static const char* options [];
public:
	FormItemModMotifOffset ( int num ) :
		FormItemSelect ( "Motif", "", getName (num), options, options [0] ) {}
	static string getName ( int num ) { return "mod_motif_offset_" + gen_itoa ( num ); }
};
const char* FormItemModMotifOffset::options [] = { "Off", "0", "-1", "+1", "-2", "+2", "-3", "+3", 0 };

class FormItemModMotif : public FormItemText {
public:
	FormItemModMotif ( int num, const string& defaultRegExp, FormValidatingJavascript* fvj ) :
		FormItemText ( "", "", getName (num), 10, 10, defaultRegExp, fvj->addMSPatternValidator ( getName (num), "Mass Mod Motif " + gen_itoa ( num ) ) ) {}
	static string getName ( int num ) { return "mod_motif_" + gen_itoa ( num ); }
};

namespace {				// MS-Bridge setting form
void printSetBridgeLinkSearchDefaults ( ostream& os )
{
	os << "function setBridgeLinkSearchDefaults ( form ) {" << endl;
	StringVector allUMods = LinkInfo::getAllUsermods ();
	for ( StringVectorSizeType a = 0 ; a < allUMods.size () ; a++ ) {
		os << "\t" << "removeSelectOption ( \"mod_AA\", '" << allUMods [a] << "' );" << endl;
	}
	os << "\t" << "var val = getSelectValue(form.link_search_type);" << endl;
	os << "\t" << "if ( val != 'No Link' ) {" << endl;

	StringVector linkNames = LinkInfo::getNameList ();
	for ( StringVectorSizeType i = 0 ; i < linkNames.size () ; i++ ) {
		if ( linkNames [i] != "No Link" ) {
			os << "\t\t" << "if ( val == '" << linkNames [i] << "' ) {" << endl;
			if ( linkNames [i] != "User Defined Link" ) {
				StringVector uMods = LinkInfo::getUsermods ( i );
				for ( StringVectorSizeType k = 0 ; k < uMods.size () ; k++ ) {
					os << "\t\t\t" << "addSelectOption ( \"mod_AA\", '" << uMods [k] << "', '', '', '', '' );" << endl;
				}
			}
			os << "\t\t" << "}" << endl;
		}
	}
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
}

namespace {				// MS-Tag setting form
void printSetMassModOptions ( ostream& os, const StringVector& linkAAList, const string& tab )
{
	for ( StringVectorSizeType i = 0 ; i < linkAAList.size () ; i++ ) {
		string linkAA = linkAAList [i];
		if ( linkAA == "Protein N-term" ) {
			os << tab << "form.mod_n_term_type.value = 'Protein';" << endl;
			os << tab << "cbSet(form." << FormItemModNTerm::getName () << ", '', true);" << endl;
		}
		else if ( linkAA == "Protein C-term" ) {
			os << tab << "form.mod_c_term_type.value = 'Protein';" << endl;
			os << tab << "cbSet(form." << FormItemModCTerm::getName () << ", '', true);" << endl;
		}
		else if ( linkAA == "Peptide N-term" ) {
			os << tab << "form.mod_n_term_type.value = 'Peptide';" << endl;
			os << tab << "cbSet(form." << FormItemModNTerm::getName () << ", '', true);" << endl;
		}
		else if ( linkAA == "Peptide C-term" ) {
			os << tab << "form.mod_c_term_type.value = 'Peptide';" << endl;
			os << tab << "cbSet(form." << FormItemModCTerm::getName () << ", '', true);" << endl;
		}
		else
			os << tab << "cbSet(form.mod_comp_ion" << ", '" << linkAA << "', true);" << endl;
	}
}
void printCBSet ( ostream& os )
{
	os << "function cbSet ( field, value, flag ) {" << endl;
	os << "\t" << "if ( field.length == undefined ) {" << endl;
	os << "\t\t" << "field.checked = flag;" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;
	os << "\t\t" << "for ( i=0; i < field.length; i++ ) {" << endl;
	os << "\t\t\t" << "if ( field[i].value == value ) {" << endl;
	os << "\t\t\t\t" << "field[i].checked = flag;" << endl;
	os << "\t\t\t" << "}" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
void printClearModAAs ( ostream& os )
{
	os << "function clearModAAs ( form ) {" << endl;
	os << "\t" << "setCheckboxes(form." << "mod_comp_ion"  << ", false);" << endl;
	os << "\t" << "form.mod_n_term_type.value = 'Peptide';" << endl;
	os << "\t" << "cbSet(form." << FormItemModNTerm::getName ()  << ", '', false);" << endl;
	os << "\t" << "form.mod_c_term_type.value = 'Peptide';" << endl;
	os << "\t" << "cbSet(form." << FormItemModCTerm::getName ()  << ", '', false);" << endl;
	os << "\t" << "cbSet(form." << FormItemModNeutralLoss::getName ()  << ", '', false);" << endl;
	os << "}" << endl;
}
void printSetUserDefinedLinkAA ( ostream& os, int n )
{
	string num = ( n == 1 ) ? "" : gen_itoa ( n ); 
	string linkSearchTypeName = "link_search_type" + num;
	os << "function setUserDefinedLinkAA ( form ) {" << endl;
	os << "\t" << "var val = getSelectValue(form." << linkSearchTypeName << ");" << endl;
	os << "\t" << "if ( val == 'User Defined Link' ) {" << endl;
		os << "\t\t" << "clearModAAs ( form );" << endl;
		os << "\t\t" << "setUserDefinedLinkAADefaults ( form );" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
void printSetUserDefinedLinkAADefaults ( ostream& os )
{
	os << "function setUserDefinedLinkAADefaults ( form ) {" << endl;
	os << "\t" << "var val = getSelectValue(form.link_aa);" << endl;
	StringVector linkOptions = LinkAAOptions::instance ().getLinkOptions ();
	for ( StringVectorSizeType i = 0 ; i < linkOptions.size () ; i++ ) {
		os << "\t" << "if ( val == '" << linkOptions [i] << "' ) {" << endl;
		printSetMassModOptions ( os, LinkInfo::getLinkAAs ( linkOptions [i] ), "\t\t" );
		os << "\t" << "}" << endl;
	}
	os << "}" << endl;

}
void printClearLinkMods ( ostream& os )
{
	os << "function clearLinkMods ( form ) {" << endl;
	StringVector allUMods = LinkInfo::getAllUsermods ();
	for ( StringVectorSizeType a = 0 ; a < allUMods.size () ; a++ ) {
		os << "\t" << "removeSelectOption ( \"msms_mod_AA\", '" << allUMods [a] << "' );" << endl;
	}
	os << "}" << endl;
}
void printSetLinkSearchDefaults ( ostream& os, int n )
{
	string num = ( n == 1 ) ? "" : gen_itoa ( n ); 
	string linkSearchTypeName = "link_search_type" + num;
	os << "function setLinkSearchDefaults ( form ) {" << endl;
	os << "\t" << "var val = getSelectValue(form." << linkSearchTypeName << ");" << endl;
	os << "\t" << "form.mod_start_nominal.value = -100;" << endl;						// First reset everything
	os << "\t" << "form.mod_end_nominal.value = 100;" << endl;
	os << "\t" << "form.msms_max_peaks.value = \"\";" << endl;
	os << "\t" << "clearModAAs ( form );" << endl;
	os << "\t" << "cbSet(form." << FormItemModUncleaved::getName ()  << ", '', false);" << endl;
	os << "\t" << "clearLinkMods(form);" << endl;
	os << "\t" << "if ( val != 'No Link' ) {" << endl;

	StringVector linkNames = LinkInfo::getNameList ();
	DoubleVector linkMasses = LinkInfo::getLinkMasses ();
	os << "\t\t" << "cbSet(form." << FormItemModUncleaved::getName ()  << ", '', true);" << endl;
	for ( StringVectorSizeType i = 0 ; i < linkNames.size () ; i++ ) {
		if ( linkNames [i] != "No Link" ) {
			os << "\t\t" << "if ( val == '" << linkNames [i] << "' ) {" << endl;
			os << "\t\t\t" << "form.mod_start_nominal.value = " << 300 + (int) ( linkMasses [i] + 0.5 ) << ";" << endl;
			os << "\t\t\t" << "form.mod_end_nominal.value = " << 3900 + (int) ( linkMasses [i] + 0.5 ) << ";" << endl;
			os << "\t\t\t" << "form.msms_max_peaks.value = 80;" << endl;
			if ( linkNames [i] == "User Defined Link" ) {
				os << "\t\t\t" << "setUserDefinedLinkAADefaults ( form );" << endl;
			}
			else {
				printSetMassModOptions ( os, LinkInfo::getLinkAAs ( i ), "\t\t\t" );
				StringVector uMods = LinkInfo::getUsermods ( i );
				for ( StringVectorSizeType k = 0 ; k < uMods.size () ; k++ ) {
					os << "\t\t\t" << "addSelectOption ( \"msms_mod_AA\", '" << uMods [k] << "', 'Common', 'All', 'Off', '' );" << endl;
					//os << "\t\t\t" << "multiSelectSet(form.msms_mod_AA" << ", '" << uMods [k] << "', true);" << endl;
				}
			}
			os << "\t\t" << "}" << endl;
		}
	}
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
void printSetModRangeTypeDefaults ( ostream& os )
{
	os << "function setModRangeTypeDefaults ( form ) {" << endl;
	os << "\t" << "var val = getSelectValue(form.mod_range_type);" << endl;
	os << "\t" << "clearModAAs ( form );" << endl;
	os << "\t" << "cbSet(form." << FormItemModUncleaved::getName ()  << ", '', false);" << endl;
	os << "\t" << "clearLinkMods(form);" << endl;
	os << "\t" << "if ( val == 'Da' ) {" << endl;
	os << "\t\t" << "form.mod_start_nominal.value = -100;" << endl;
	os << "\t\t" << "form.mod_end_nominal.value = 100;" << endl;
	os << "\t\t" << "hidediv ( 'div_mm_1' );" << endl;
	os << "\t\t" << "showdiv ( 'div_mm_2' );" << endl;
	os << "\t\t" << "showdiv ( 'div_mm_3' );" << endl;
	os << "\t\t" << "showdiv ( 'div_mm_4' );" << endl;
	os << "\t\t" << "showdiv ( 'div_cl' );" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;
	os << "\t\t" << "form.mod_start_nominal.value = -2;" << endl;
	os << "\t\t" << "form.mod_end_nominal.value = 2;" << endl;
	os << "\t\t" << "showdiv ( 'div_mm_1' );" << endl;
	os << "\t\t" << "hidediv ( 'div_mm_2' );" << endl;
	os << "\t\t" << "hidediv ( 'div_mm_3' );" << endl;
	os << "\t\t" << "hidediv ( 'div_mm_4' );" << endl;
	os << "\t\t" << "hidediv ( 'div_cl' );" << endl;
	os << "\t" << "cbSet(form." << FormItemModNeutralLoss::getName ()  << ", '', true);" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
}
void msBridgeJavascriptFunctions ( ostream& os )
{
	startJavascript ( os );

	printSetBridgeLinkSearchDefaults ( os );

	endJavascript ( os );
}
void massModCrosslinkingJavascriptFunctions ( ostream& os, int n )
{
	startJavascript ( os );

	printCBSet ( os );
	printClearModAAs ( os );
	printSetUserDefinedLinkAADefaults ( os );
	printSetUserDefinedLinkAA ( os, n );
	printClearLinkMods ( os );
	printSetLinkSearchDefaults ( os, n );
	printSetModRangeTypeDefaults ( os );

	endJavascript ( os );
}
const int MassModificationForm::numMotifs = 2;
MassModificationForm::MassModificationForm ( const VectorConstParameterListPtr& params, FormValidatingJavascript* fvj ) :
	compIonForm ( params, "ACDEFGHIKLMNPQRSTVWY", "mod_" ),
	fvj ( fvj )
{
	create ( params );
}
void MassModificationForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemModRangeType::getName ()]->setValue ( p );
		formItemMap [FormItemModStartNominal::getName ()]->setValue ( p );
		formItemMap [FormItemModEndNominal::getName ()]->setValue ( p );
		formItemMap [FormItemModDefect::getName ()]->setValue ( p );
		formItemMap [FormItemModMaxZ::getName ()]->setValue ( p );
		formItemMap [FormItemModUncleaved::getName ()]->setValue ( p );
		formItemMap [FormItemModNTermType::getName ()]->setValue ( p );
		formItemMap [FormItemModNTerm::getName ()]->setValue ( p );
		formItemMap [FormItemModCTermType::getName ()]->setValue ( p );
		formItemMap [FormItemModCTerm::getName ()]->setValue ( p );
		formItemMap [FormItemModNeutralLoss::getName ()]->setValue ( p );
		formItemMap [FormItemModRare::getName ()]->setValue ( p );
		for ( int i = 1 ; i <= numMotifs ; i++ ) {
			formItemMap [FormItemModMotifAA::getName (i)]->setValue ( p );
			formItemMap [FormItemModMotifOffset::getName (i)]->setValue ( p );
			formItemMap [FormItemModMotif::getName (i)]->setValue ( p );
		}
	}
}
void MassModificationForm::createItems ()
{
	formItemMap [FormItemModRangeType::getName ()]		= new FormItemModRangeType ();
	formItemMap [FormItemModStartNominal::getName ()]	= new FormItemModStartNominal ();
	formItemMap [FormItemModEndNominal::getName ()]		= new FormItemModEndNominal ();
	formItemMap [FormItemModDefect::getName ()]			= new FormItemModDefect ();
	formItemMap [FormItemModMaxZ::getName ()]			= new FormItemModMaxZ ();
	formItemMap [FormItemModUncleaved::getName ()]		= new FormItemModUncleaved ();
	formItemMap [FormItemModNTermType::getName ()]		= new FormItemModNTermType ();
	formItemMap [FormItemModNTerm::getName ()]			= new FormItemModNTerm ();
	formItemMap [FormItemModCTermType::getName ()]		= new FormItemModCTermType ();
	formItemMap [FormItemModCTerm::getName ()]			= new FormItemModCTerm ();
	formItemMap [FormItemModNeutralLoss::getName ()]	= new FormItemModNeutralLoss ();
	formItemMap [FormItemModRare::getName ()]			= new FormItemModRare ();
	for ( int i = 1 ; i <= numMotifs ; i++ ) {
		formItemMap [FormItemModMotifAA::getName (i)]		= new FormItemModMotifAA ( i, i == 1 ? "N" : "S" );
		formItemMap [FormItemModMotifOffset::getName (i)]	= new FormItemModMotifOffset ( i );
		formItemMap [FormItemModMotif::getName (i)]			= new FormItemModMotif ( i, i == 1 ? "N[^P][ST]" : ".", fvj );
	}
}
void MassModificationForm::setMMVisualizationFlags ( const string& val, bool& div_mm_1, bool& div_mm_2, bool& div_mm_3, bool& div_mm_4 ) const
{
	if ( val == "Da" ) {
		div_mm_1 = false;
		div_mm_2 = true;
		div_mm_3 = true;
		div_mm_4 = true;
	}
	else {
		div_mm_1 = true;
		div_mm_2 = false;
		div_mm_3 = false;
		div_mm_4 = false;
	}
}
void MassModificationForm::printHTML ( ostream& os )
{
	ExpandableJavascriptBlock ejb ( "Mass Modifications" );
	ejb.printHeader ( os );
	tableStart ( os, false );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				os << "Range ";
				formItemMap [FormItemModRangeType::getName ()]->printHTML ( os );
				os << " ";
				formItemMap [FormItemModStartNominal::getName ()]->printHTML ( os );
				os << "to ";
				formItemMap [FormItemModEndNominal::getName ()]->printHTML ( os );
				formItemMap [FormItemModDefect::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
			bool div_mm_1, div_mm_2, div_mm_3, div_mm_4;
			setMMVisualizationFlags ( formItemMap [FormItemModRangeType::getName ()]->getValue ( 0 ), div_mm_1, div_mm_2, div_mm_3, div_mm_4 );
			tableHeaderStart ( os, "", "center", true );
				divStart ( os, "div_mm_1", div_mm_1 );
					formItemMap [FormItemModMaxZ::getName ()]->printHTML ( os );
				divEnd ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "left", true );
				divStart ( os, "div_mm_2", div_mm_2 );
					StringVector sv;
					sv.push_back ( "mod_comp_ion" );
					sv.push_back ( FormItemModNTerm::getName () );
					sv.push_back ( FormItemModCTerm::getName () );
					sv.push_back ( FormItemModNeutralLoss::getName () );
					CheckboxSettingJavascript csj ( sv );
					csj.print ( os );
				divEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true, 4 );
				divStart ( os, "div_mm_3", div_mm_3 );
					compIonForm.printHTML ( os );
				divEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true, 3 );
				divStart ( os, "div_mm_4", div_mm_4 );
					formItemMap [FormItemModNTermType::getName ()]->printHTML ( os );
					formItemMap [FormItemModNTerm::getName ()]->printHTML ( os );
					formItemMap [FormItemModCTermType::getName ()]->printHTML ( os );
					formItemMap [FormItemModCTerm::getName ()]->printHTML ( os );
					formItemMap [FormItemModUncleaved::getName ()]->printHTML ( os );
				divEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "left", true );
				formItemMap [FormItemModNeutralLoss::getName ()]->printHTML ( os );
				formItemMap [FormItemModRare::getName ()]->printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
		for ( int i = 1 ; i <= numMotifs ; i++ ) {
			tableRowStart ( os );
				tableHeaderStart ( os, "", "left", true );
					formItemMap [FormItemModMotifAA::getName (i)]->printHTML ( os );
					formItemMap [FormItemModMotifOffset::getName (i)]->printHTML ( os );
					formItemMap [FormItemModMotif::getName (i)]->printHTML ( os );
				tableHeaderEnd ( os );
			tableRowEnd ( os );
		}
	tableEnd ( os );
	ejb.printFooter ( os );
	showHiddenItems ( os );
}
void MassModificationForm::printCGI ( ostream& os ) const
{
	ProspectorForm::printCGI ( os );
	compIonForm.printCGI ( os );
}
void MassModificationForm::printHTMLJavascriptHidden ( ostream& os ) const
{
	ProspectorForm::printHTMLJavascriptHidden ( os );
	compIonForm.printHTMLJavascriptHidden ( os );
}
bool MassModificationForm::getDivCL ()
{
	string vv1 = formItemMap [FormItemModRangeType::getName ()]->getValue ( 0 );
	bool div_cl = true;
	if ( vv1 != "Da" ) div_cl = true;
	return div_cl;
}

class VariableModsSettingJavascript {
protected:
	string prefix;
	StringVector typeOptions;
	string phosphoSTYMenuString;
	void printAddSelectOption ( ostream& os ) const;
	void printRemoveSelectOption ( ostream& os ) const;
	void printChangeVariableModType ( ostream& os ) const;
	void addSingleVariableModType ( ostream& os, const string& type, bool first ) const;
	void printRemoveOption ( ostream& os ) const;
private:
	FormItemSelect fis1;
	FormItemSelectMultiple fis2;
	virtual void printFunctions ( ostream& os ) const;
	virtual void printAddOption ( ostream& os ) const;
public:
	VariableModsSettingJavascript ( const string& prefix = "" );
	virtual void setValues ( const ParameterList* p ) {}
	virtual void printHTML ( ostream& os ) const;
};

class MSTagVariableModsSettingJavascript : public VariableModsSettingJavascript {
	FormItemLimit fil;
	//FormItemModFile fimf;
	FormItemMotifOffset mf;
	FormItemMotif motif;
	virtual void printFunctions ( ostream& os ) const;
	virtual void printAddOption ( ostream& os ) const;
	void printRemoveMaxOptions ( ostream& os ) const;
	void printAddMaxOptions ( ostream& os ) const;
	void printSetMaxOptions ( ostream& os ) const;
public:
	MSTagVariableModsSettingJavascript ( FormValidatingJavascript* fvj );
	virtual void setValues ( const ParameterList* p );
	virtual void printHTML ( ostream& os ) const;
};

VariableModsSettingJavascript::VariableModsSettingJavascript ( const string& prefix ) :
	prefix ( prefix ),
	typeOptions ( Usermod::getTypeOptions () ),
	phosphoSTYMenuString ( Usermod::getPhosphoSTYMenuString () ),
	fis1 ( "", "", prefix + "mod_AA_types", typeOptions, phosphoSTYMenuString, 1, "changeVariableModType (this.form)" ),
	fis2 ( "", "", prefix + "mod_AA_list", Usermod::getNames (phosphoSTYMenuString), StringVector (), 4 )
{
}
void VariableModsSettingJavascript::printHTML ( ostream& os ) const
{
	printFunctions ( os );
	os << "Add (+) or remove (-) mods using menus & buttons below<br />" << endl;
	tableStart ( os );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				fis1.printHTML ( os );
			tableHeaderEnd ( os );
			tableHeaderStart ( os, "", "center", true );
				fis2.printHTML ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );
	os << "<input type=\"button\" value=\"+\" onclick=\"" << "addOption (this.form)\" />" << endl;
	os << "<input type=\"button\" value=\"-\" onclick=\"" << "removeOption (this.form)\" />" << endl;
}
void VariableModsSettingJavascript::printFunctions ( ostream& os ) const
{
	startJavascript ( os );
	printAddSelectOption ( os );
	printRemoveSelectOption ( os );
	printChangeVariableModType ( os );
	printAddOption ( os );
	printRemoveOption ( os );
	endJavascript ( os );
}
void VariableModsSettingJavascript::printAddSelectOption ( ostream& os ) const
{
	os << "function addSelectOption" << " ( name, val1, val2, val3, val4, val5 ) {" << endl;
	os << "\t" << "var x = document.getElementsByName(name) [0];" << endl;
	os << "\t" << "var pos = x.length;" << endl;		// Default is to add at the end
	os << "\t" << "for ( var i = 0; i < x.length ; i++ ) {" << endl;
	os << "\t\t" << "var text = x.options [i].text;" << endl;
	os << "\t\t" << "var idx = text.lastIndexOf ( \")\" );" << endl;
	os << "\t\t" << "if ( idx != text.length - 1 ) {" << endl;
	os << "\t\t\t" << "text = text.substring ( 0, idx+1 );" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t\t" << "if ( text == val1 ) {" << endl;
	os << "\t\t\t" << "x.remove(i);" << endl;
	os << "\t\t\t" << "pos = i;" << endl;
	os << "\t\t\t" << "break;" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t\t" << "else if ( text > val1 ) {" << endl;
	os << "\t\t\t" << "pos = i;" << endl;
	os << "\t\t\t" << "break;" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "var op = document.createElement(\"option\");" << endl;
	os << "\t" << "var v = val1;" << endl;
	os << "\t" << "if ( val2 != \"Common\" && val2 != \"\" ) {" << endl;
	os << "\t\t" << "v += \" - \" + val2;" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "if ( val3 != \"All\" && val3 != \"\" ) {" << endl;
	os << "\t\t" << "v += \" - \" + val3;" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "if ( val4 != \"Off\" && val4 != \"\" && val5 != \"\" ) {" << endl;
	os << "\t\t" << "v += \" - Motif \" + val4 + \" \" + val5;" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "op.text = v;" << endl;
	os << "\t" << "x.add(op, pos);" << endl;
	os << "\t" << "x.options [pos].selected = true;" << endl;
	os << "}" << endl;
}
void VariableModsSettingJavascript::printRemoveSelectOption ( ostream& os ) const
{
	os << "function removeSelectOption" << " ( name, val1 ) {" << endl;
	os << "\t" << "var x = document.getElementsByName(name) [0];" << endl;
	os << "\t" << "for ( var i = 0; i < x.length ; i++ ) {" << endl;
	os << "\t\t" << "var text = x.options [i].text;" << endl;
	os << "\t\t" << "var idx = text.lastIndexOf ( \")\" );" << endl;
	os << "\t\t" << "if ( idx != text.length - 1 ) {" << endl;
	os << "\t\t\t" << "text = text.substring ( 0, idx+1 );" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t\t" << "if ( text == val1 ) {" << endl;
	os << "\t\t\t" << "x.remove(i);" << endl;
	os << "\t\t\t" << "return;" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
void VariableModsSettingJavascript::printChangeVariableModType ( ostream& os ) const
{
	os << "function changeVariableModType ( form ) {" << endl;
	os << "\t" << "var val = getSelectValue(form." + prefix + "mod_AA_types);" << endl;
	os << "\t" << "var x = document.getElementsByName(\"" + prefix + "mod_AA_list\") [0];" << endl;
	os << "\t" << "removeOptionsFromSelect ( x );" << endl;
	for ( int i = 0 ; i < typeOptions.size () ; i++ ) {
		addSingleVariableModType ( os, typeOptions [i], i == 0 );
	}
	os << "}" << endl;
}
void VariableModsSettingJavascript::addSingleVariableModType ( ostream& os, const string& type, bool first ) const
{
	StringVector mods = Usermod::getNames ( type );
	os << "\t";
	if ( !first ) os << "else ";
	os << "if ( val == '" << type << "' ) {" << endl;
	for ( int i = 0 ; i < mods.size () ; i++ ) {
		os << "\t\t" << "addOptionToSelect ( x, \"" << mods [i] << "\" );" << endl;
	}
	os << "\t" << "}" << endl;
}
void VariableModsSettingJavascript::printAddOption ( ostream& os ) const
{
	os << "function addOption" << " ( form ) {" << endl;
	os << "\t" << "var x = document.getElementsByName(\"" + prefix + "mod_AA_list\") [0];" << endl;
	os << "\t" << "for ( var i = 0; i < x.length ; i++ ) {" << endl;
	os << "\t\t" << "if ( x.options [i].selected ) {" << endl;
	os << "\t\t\t" << "addSelectOption ( \"" + prefix + "mod_AA\", x.options [i].text, \"\", \"\", \"\", \"\" );" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
void VariableModsSettingJavascript::printRemoveOption ( ostream& os ) const
{
	os << "function removeOption" << " ( form ) {" << endl;
	os << "\t" << "var x = document.getElementsByName(\"" + prefix + "mod_AA_list\") [0];" << endl;
	os << "\t" << "for ( var i = 0; i < x.length ; i++ ) {" << endl;
	os << "\t\t" << "if ( x.options [i].selected ) {" << endl;
	os << "\t\t\t" << "removeSelectOption ( \"" + prefix + "mod_AA\", x.options [i].text );" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
}

MSTagVariableModsSettingJavascript::MSTagVariableModsSettingJavascript ( FormValidatingJavascript* fvj ) :
	VariableModsSettingJavascript ( "msms_" ),
	fil ( "mod_AA", 1 ),
	//fimf (),
	mf (),
	motif ( fvj )
{
}
void MSTagVariableModsSettingJavascript::setValues ( const ParameterList* p )
{
	fil.setOptions ( p, "msms_max_modifications" );
}
void MSTagVariableModsSettingJavascript::printHTML ( ostream& os ) const
{
	VariableModsSettingJavascript::printHTML ( os );
	fil.printHTML ( os );
	//fimf.printHTML ( os );
	mf.printHTML ( os );
	motif.printHTML ( os );
}
void MSTagVariableModsSettingJavascript::printFunctions ( ostream& os ) const
{
	startJavascript ( os );
	printAddSelectOption ( os );
	printRemoveSelectOption ( os );
	printChangeVariableModType ( os );
	printAddOption ( os );
	printRemoveOption ( os );
	printRemoveMaxOptions ( os );
	printAddMaxOptions ( os );
	printSetMaxOptions ( os );
	endJavascript ( os );
}
void MSTagVariableModsSettingJavascript::printAddOption ( ostream& os ) const
{
	os << "function addOption" << " ( form ) {" << endl;
	os << "\t" << "var x = document.getElementsByName(\"" + prefix + "mod_AA_list\") [0];" << endl;
	os << "\t" << "var val2 = getSelectValue(form.mod_AA_limit);" << endl;
	os << "\t" << "var val3 = \"All\";" << endl;
	os << "\t" << "var val4 = getSelectValue(form.motif_offset);" << endl;
	os << "\t" << "var val5 = form.motif.value;" << endl;
	os << "\t" << "for ( var i = 0; i < x.length ; i++ ) {" << endl;
	os << "\t\t" << "if ( x.options [i].selected ) {" << endl;
	os << "\t\t\t" << "addSelectOption ( \"" + prefix + "mod_AA\", x.options [i].text, val2, val3, val4, val5 );" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
void MSTagVariableModsSettingJavascript::printRemoveMaxOptions ( ostream& os ) const
{
	os << "function removeMaxOptions" << " ( name ) {" << endl;
	os << "\t" << "var x = document.getElementsByName(name) [0];" << endl;
	os << "\t" << "var len1 = x.length;" << endl;
	os << "\t" << "for ( var i = 0 ; i < len1 ; i++ ) {" << endl;
	os << "\t\t" << "var text = x.options [i].text;" << endl;
	os << "\t\t" << "if ( text.substring ( 0, 3 ) == \"Max\" ) {" << endl;
	os << "\t\t\t" << "for ( var j = len1 - 1 ; j >= i ; j-- ) {" << endl;
	os << "\t\t\t\t" << "x.remove(j);" << endl;
	os << "\t\t\t" << "}" << endl;
	os << "\t\t\t" << "break;" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
void MSTagVariableModsSettingJavascript::printAddMaxOptions ( ostream& os ) const
{
	os << "function addMaxOptions" << " ( name, num ) {" << endl;
	os << "\t" << "var x = document.getElementsByName(name) [0];" << endl;
	os << "\t" << "for ( var i = 1 ; i <= num ; i++ ) {" << endl;
	os << "\t\t" << "var op = document.createElement(\"option\");" << endl;
	os << "\t\t" << "op.text = \"Max \" + i.toString ();" << endl;
	os << "\t\t" << "x.add(op);" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
void MSTagVariableModsSettingJavascript::printSetMaxOptions ( ostream& os ) const
{
	os << "function setMaxOptions" << " ( form ) {" << endl;
	os << "\t" << "var num = form.msms_max_modifications.value - 1;" << endl;
	os << "\t" << "removeMaxOptions ( \"mod_AA_limit\" );" << endl;
	os << "\t" << "addMaxOptions ( \"mod_AA_limit\", num );" << endl;
	os << "\t" << "for ( var i = 1 ; i <= 6 ; i++ ) {" << endl;
	os << "\t\t" << "removeMaxOptions ( \"mod_\" + i.toString () + \"_limit\" );" << endl;
	os << "\t\t" << "addMaxOptions ( \"mod_\" + i.toString () + \"_limit\", num );" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
VariableModsForm::VariableModsForm ( const VectorConstParameterListPtr& params, const string& prefix, FormValidatingJavascript* fvj ) :
	prefix ( prefix ),
	vmsj ( ( prefix == "msms_" ) ? new MSTagVariableModsSettingJavascript ( fvj ) : new VariableModsSettingJavascript )
{
	create ( params );
}
VariableModsForm::~VariableModsForm ()
{
	delete vmsj;
}
void VariableModsForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemModAAList::getName(prefix)]->setOptions ( p );
		formItemMap [FormItemModAAList::getName(prefix)]->setValue ( p );
		vmsj->setValues ( p );
	}
}
void VariableModsForm::createItems ()
{
	if ( prefix == "msms_" )
		formItemMap [FormItemModAAList::getName (prefix)] = new FormItemModAAList ( 6 );
	else																						// MS-Digest
		formItemMap [FormItemModAAList::getName ()]	= new FormItemModAAList ( false, true, 5 );

}
void VariableModsForm::printHTML ( ostream& os )
{
	formItemMap [FormItemModAAList::getName (prefix)]->printHTML ( os );
	vmsj->printHTML ( os );
}
int ExtraUserModsForm::numUserMods = 6;
ExtraUserModsForm::ExtraUserModsForm ( const VectorConstParameterListPtr& params, FormValidatingJavascript* fvj, bool massItem, int num ) :
	fvj ( fvj ),
	massItem ( massItem ),
	num ( num )
{
	create ( params );
}
void ExtraUserModsForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		for ( int i = 1 ; i <= numUserMods ; i++ ) {
			string sNum = gen_itoa ( i );
			formItemMap [FormItemLabel::getName("mod_"+sNum, num)]->setValue ( p );
			formItemMap [FormItemSpecificity::getName("mod_"+sNum, num)]->setValue ( p );
			formItemMap [FormItemComposition::getName("mod_"+sNum, num)]->setValue ( p );
			if ( massItem ) {
				formItemMap [FormItemAccurateMass::getName("mod_"+sNum, num)]->setValue ( p );
				formItemMap [FormItemLimit::getName("mod_"+sNum, num)]->setOptions ( p, "msms_max_modifications" );
				formItemMap [FormItemLimit::getName("mod_"+sNum, num)]->setValue ( p );
				formItemMap [FormItemMotifOffset::getName("mod_"+sNum, num)]->setValue ( p );
				formItemMap [FormItemMotif::getName("mod_"+sNum, num)]->setValue ( p );
			}
		}
	}
}
void ExtraUserModsForm::createItems ()
{
	for ( int i = 1 ; i <= numUserMods ; i++ ) {
		string sNum = gen_itoa ( i );
		formItemMap [FormItemLabel::getName("mod_"+sNum, num)]		= new FormItemLabel ( "Mod "+sNum, "mod_"+sNum, num, "" );
		formItemMap [FormItemSpecificity::getName("mod_"+sNum, num)]= new FormItemSpecificity ( "mod_"+sNum, num );
		formItemMap [FormItemComposition::getName("mod_"+sNum, num)]= new FormItemComposition ( "Mod "+sNum, "mod_"+sNum, num, "", fvj );
		if ( massItem ) {
			formItemMap [FormItemAccurateMass::getName("mod_"+sNum, num)]= new FormItemAccurateMass ( "mod_"+sNum, num, fvj );
			formItemMap [FormItemLimit::getName("mod_"+sNum, num)]= new FormItemLimit ( "mod_"+sNum, num );
			formItemMap [FormItemMotifOffset::getName("mod_"+sNum, num)]= new FormItemMotifOffset ( "mod_"+sNum, num );
			formItemMap [FormItemMotif::getName("mod_"+sNum, num)]= new FormItemMotif ( fvj, "mod_"+sNum, num );
		}
	}
}
void ExtraUserModsForm::printHTML ( ostream& os )
{
	ExpandableJavascriptBlock ejb ( "<a href=\"../html/instruct/allman.htm#user_defined_variable_mods\">User Defined Variable Modifications</a>" );
	ejb.printHeader ( os );
	for ( int i = 1 ; i <= numUserMods ; i++ ) {
		string sNum = gen_itoa ( i );
		formItemMap [FormItemLabel::getName("mod_"+sNum, num)]->printHTML ( os );
		formItemMap [FormItemSpecificity::getName("mod_"+sNum, num)]->printHTML ( os );
		if ( massItem ) os << "<br />" << endl;
		formItemMap [FormItemComposition::getName("mod_"+sNum, num)]->printHTML ( os );
		if ( massItem ) {
			formItemMap [FormItemAccurateMass::getName("mod_"+sNum, num)]->printHTML ( os );
		}
		os << "<br />" << endl;
		if ( massItem ) {
			formItemMap [FormItemLimit::getName("mod_"+sNum, num)]->printHTML ( os );
			formItemMap [FormItemMotifOffset::getName("mod_"+sNum, num)]->printHTML ( os );
			formItemMap [FormItemMotif::getName("mod_"+sNum, num)]->printHTML ( os );
			os << "<br />" << endl;
		}
	}
	ejb.printFooter ( os );
	showHiddenItems ( os );
}
void ExtraUserModsForm::replaceSubstrings ( string& s ) const
{
	for ( int i = 1 ; i <= numUserMods ; i++ ) {
		string sNum = gen_itoa ( i );
		s = genReplaceSubstrings ( s, FormItemLabel::getName("mod_"+sNum, num),				FormItemLabel::getName("mod_"+sNum, 1) );
		s = genReplaceSubstrings ( s, FormItemSpecificity::getName("mod_"+sNum, num),		FormItemSpecificity::getName("mod_"+sNum, 1) );
		s = genReplaceSubstrings ( s, FormItemComposition::getName("mod_"+sNum, num),		FormItemComposition::getName("mod_"+sNum, 1) );
		if ( massItem ) {
			s = genReplaceSubstrings ( s, FormItemAccurateMass::getName("mod_"+sNum, num),	FormItemAccurateMass::getName("mod_"+sNum, 1) );
			s = genReplaceSubstrings ( s, FormItemLimit::getName("mod_"+sNum, num),			FormItemLimit::getName("mod_"+sNum, 1) );
			s = genReplaceSubstrings ( s, FormItemMotifOffset::getName("mod_"+sNum, num),	FormItemMotifOffset::getName("mod_"+sNum, 1) );
			s = genReplaceSubstrings ( s, FormItemMotif::getName("mod_"+sNum, num),			FormItemMotif::getName("mod_"+sNum, 1) );
		}
	}
}
void ExtraUserModsForm::copyToCGI ( ostream& os, const ParameterList* params )
{
	for ( int i = 1 ; i <= numUserMods ; i++ ) {
		string sNum = gen_itoa ( i );
		string lab = params->getStringValue ( FormItemLabel::getName("mod_"+sNum, 1) );
		if ( !lab.empty () ) {
			params->copyToCGI ( os, FormItemLabel::getName("mod_"+sNum, 1) );
			params->copyToCGI ( os, FormItemSpecificity::getName("mod_"+sNum, 1) );
			params->copyToCGI ( os, FormItemComposition::getName("mod_"+sNum, 1) );
		}
	}
}

int CrosslinkingForm::numUserMods = 6;
CrosslinkingForm::CrosslinkingForm ( const VectorConstParameterListPtr& params, FormValidatingJavascript* fvj, bool tagForm, int num ) :
	fvj ( fvj ),
	tagForm ( tagForm ),
	num ( num )
{
	create ( params );
}
void CrosslinkingForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemLinkSearchType::getName (num)]->setValue ( p );
		if ( tagForm ) {
			formItemMap [FormItemMaxSavedTagHits::getName ()]->setValue ( p );
		}
		else {
			formItemMap [FormItemMaxLinkMolecules::getName ()]->setValue ( p );
		}
		formItemMap [FormItemLinkAA::getName (num)]->setValue ( p );
		formItemMap [FormItemComposition::getName("bridge", num)]->setValue ( p );
		if ( !tagForm ) {
			for ( int i = 1 ; i <= numUserMods ; i++ ) {
				string sNum = gen_itoa ( i );
				formItemMap [FormItemLabel::getName("mod_"+sNum)]->setValue ( p );
				formItemMap [FormItemAAModified::getName(sNum)]->setValue ( p );
				formItemMap [FormItemComposition::getName("mod_"+sNum)]->setValue ( p );
			}
		}
	}
}
void CrosslinkingForm::createItems ()
{
	if ( tagForm ) {
		formItemMap [FormItemLinkSearchType::getName (num)] = new FormItemLinkSearchType ( "No Link", "setLinkSearchDefaults( this.form )", num );
		formItemMap [FormItemMaxSavedTagHits::getName ()]	= new FormItemMaxSavedTagHits ( fvj );
	}
	else {
		formItemMap [FormItemLinkSearchType::getName (num)] = new FormItemLinkSearchType ( "", "setBridgeLinkSearchDefaults( this.form )", num );
		formItemMap [FormItemMaxLinkMolecules::getName ()] = new FormItemMaxLinkMolecules (fvj);
	}
	formItemMap [FormItemLinkAA::getName (num)]			= new FormItemLinkAA ( "setUserDefinedLinkAA( this.form )", num );
	formItemMap [FormItemComposition::getName("bridge", num)]= new FormItemComposition ( "Bridge", "bridge", num, "", fvj );
	if ( !tagForm ) {
		for ( int i = 1 ; i <= numUserMods ; i++ ) {
			string sNum = gen_itoa ( i );
			formItemMap [FormItemLabel::getName("mod_"+sNum)]		= new FormItemLabel ( "Mod "+sNum, "mod_"+sNum, 1, "" );
			formItemMap [FormItemAAModified::getName(sNum)]			= new FormItemAAModified ( sNum, i <= 3 ? "K" : "Protein N-term" );
			formItemMap [FormItemComposition::getName("mod_"+sNum)]	= new FormItemComposition ( "Mod "+sNum, "mod_"+sNum, 1, "", fvj );
		}
	}
}
void CrosslinkingForm::printHTML ( ostream& os )
{
	if ( tagForm ) {
		ExpandableJavascriptBlock ejb ( "Crosslinking" );
		ejb.printHeader ( os );
		formItemMap [FormItemLinkSearchType::getName (num)]->printHTML ( os );
		os << "<br />" << endl;
		formItemMap [FormItemMaxSavedTagHits::getName ()]->printHTML ( os );
		os << "<br />" << endl;
		printHTMLUDLP ( os );
		ejb.printFooter ( os );
	}
	else {
		formItemMap [FormItemLinkSearchType::getName (num)]->printHTML ( os );
		os << "<br />" << endl;
		formItemMap [FormItemMaxLinkMolecules::getName ()]->printHTML ( os );
		os << "<br />" << endl;
		printHTMLUDLP ( os );
	}
	showHiddenItems ( os );
}
void CrosslinkingForm::printHTMLUDLP ( ostream& os )
{
	ExpandableJavascriptBlock ejb ( "User Defined Link Parameters<br /><a href=\"../html/instruct/bridgeman.htm#user_link_params\">(see instructions)</a>" );
	ejb.printHeader ( os );
	formItemMap [FormItemLinkAA::getName (num)]->printHTML ( os );
	os << "<br />" << endl;
	formItemMap [FormItemComposition::getName("bridge", num)]->printHTML ( os );
	os << "<br />" << endl;
	if ( !tagForm ) {
		for ( int i = 1 ; i <= numUserMods ; i++ ) {
			string sNum = gen_itoa ( i );
			formItemMap [FormItemLabel::getName("mod_"+sNum)]->printHTML ( os );
			formItemMap [FormItemAAModified::getName(sNum)]->printHTML ( os );
			os << "<br />" << endl;
			formItemMap [FormItemComposition::getName("mod_"+sNum)]->printHTML ( os );
			os << "<br />" << endl;
		}
	}
	ejb.printFooter ( os );
}
MatrixModeForm::MatrixModeForm ( const VectorConstParameterListPtr& params, const string& prefix ) :
	prefix ( prefix ),
	names ( ModificationTable::getNames () )
{
	create ( params );
}
void MatrixModeForm::createItems ()
{
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		string n = names [i];
		formItemMap [getName (prefix) + gen_itoa (i)]	= new FormItemCheckbox ( n, "", getName (prefix), false, n );
	}
}
void MatrixModeForm::printHTML ( ostream& os )
{
	ExpandableJavascriptBlock* ejb = 0;
	
	if ( prefix == "msms" ) {
		ejb = new ExpandableJavascriptBlock ( "Matrix Modifications" );
		ejb->printHeader ( os );
	}
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		formItemMap [getName (prefix) + gen_itoa (i)]->printHTML ( os );
		if ( i+1 % 3 == 0 ) os << "<br />" << endl;
	}
	if ( prefix == "msms" ) {
		ejb->printFooter ( os );
		delete ejb;
	}
}
string MatrixModeForm::getName ( const string& ms )
{
	return ms + string ( "_search_type" );
}
void MatrixModeForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
			formItemMap [getName (prefix) + gen_itoa (i)]->setValue ( p );
		}
	}
}
#ifdef MYSQL_DATABASE
class FormItemProjectShowFilter : public FormItemSelect {
	static const char* options [];
public:
	FormItemProjectShowFilter ( const string& changeFunction );
	static string getName () { return "date_filter"; }
	static string getDefault () { return options [0]; }
};
const char* FormItemProjectShowFilter::options [] = { "All Projects", "Projects Created Between", "Projects Accessed Between", 0 };
FormItemProjectShowFilter::FormItemProjectShowFilter ( const string& changeFunction ) :
	FormItemSelect ( "Show", "", getName (), options, options [0], 1, changeFunction )
{
}

class FormItemStartYear : public FormItemSelect {
public:
	FormItemStartYear ( const string& firstYear, const string& label = "" );
	static string getName () { return "start_year"; }
	static string getDefault () { return genCurrentYear (); }
};
FormItemStartYear::FormItemStartYear ( const string& firstYear, const string& label ) :
	FormItemSelect ( label, "", getName (), getYearList ( firstYear, true ), genCurrentYear () )
{
}

class FormItemEndYear : public FormItemSelect {
public:
	FormItemEndYear ( const string& firstYear, const string& label = "" );
	static string getName () { return "end_year"; }
	static string getDefault () { return genCurrentYear (); }
};
FormItemEndYear::FormItemEndYear ( const string& firstYear, const string& label ) :
	FormItemSelect ( label, "", getName (), getYearList ( firstYear, true ), genCurrentYear () )
{
}

class FormItemStartMonth : public FormItemSelect {
	static const char* options [];
public:
	FormItemStartMonth ( const string& label = "" );
	static string getName () { return "start_month"; }
	static string getDefault () { return genCurrentIntMonth (); }
};
const char* FormItemStartMonth::options [] = { "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", 0 };
FormItemStartMonth::FormItemStartMonth ( const string& label ) :
	FormItemSelect ( label, "", getName (), options, genCurrentIntMonth () )
{
}
class FormItemEndMonth : public FormItemSelect {
	static const char* options [];
public:
	FormItemEndMonth ( const string& label = "" );
	static string getName () { return "end_month"; }
	static string getDefault () { return genCurrentIntMonth (); }
};
const char* FormItemEndMonth::options [] = { "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", 0 };
FormItemEndMonth::FormItemEndMonth ( const string& label ) :
	FormItemSelect ( label, "", getName (), options, genCurrentIntMonth () )
{
}

ProjectDateFilterForm::ProjectDateFilterForm ( const ParameterList* paramList )
{
	string dateFilter	= paramList->getStringValue ( FormItemProjectShowFilter::getName (),FormItemProjectShowFilter::getDefault () );
	string startYear	= paramList->getStringValue ( FormItemStartYear::getName (),		FormItemStartYear::getDefault () );
	string endYear		= paramList->getStringValue ( FormItemEndYear::getName (),			FormItemEndYear::getDefault () );
	string startMonth	= paramList->getStringValue ( FormItemStartMonth::getName (),		FormItemStartMonth::getDefault () );
	string endMonth		= paramList->getStringValue ( FormItemEndMonth::getName (),			FormItemEndMonth::getDefault () );
	pdf = new ProjectDateFilter ( dateFilter, startYear, endYear, startMonth, endMonth );
	VectorConstParameterListPtr pv;
	pv.push_back ( paramList );
	create ( pv );
}
ProjectDateFilterForm::~ProjectDateFilterForm ()
{
	delete pdf;
}
void ProjectDateFilterForm::createItems ()
{
	formItemMap [FormItemProjectShowFilter::getName ()]	= new FormItemProjectShowFilter ( "showDateFilterItems( this.form.date_filter )" );
	formItemMap [FormItemStartYear::getName ()]			= new FormItemStartYear ( "2007" );
	formItemMap [FormItemEndYear::getName ()]			= new FormItemEndYear ( "2007" );
	formItemMap [FormItemStartMonth::getName ()]		= new FormItemStartMonth ( "between" );
	formItemMap [FormItemEndMonth::getName ()]			= new FormItemEndMonth ( "and" );
}
void ProjectDateFilterForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemProjectShowFilter::getName ()]->setValue ( p );
		formItemMap [FormItemStartYear::getName ()]->setValue ( p );
		formItemMap [FormItemEndYear::getName ()]->setValue ( p );
		formItemMap [FormItemStartMonth::getName ()]->setValue ( p );
		formItemMap [FormItemEndMonth::getName ()]->setValue ( p );
	}
}
void ProjectDateFilterForm::printHTML ( ostream& os, bool showProjectsDiv )
{
	tableRowStartDiv ( os, "show_projects", showProjectsDiv );
		tableHeaderStart ( os, "", "center", true );
			formItemMap [FormItemProjectShowFilter::getName ()]->printHTML ( os );
			bool divInfo2 = true;
			if ( formItemMap [FormItemProjectShowFilter::getName ()]->getValue ( 0 ) == FormItemProjectShowFilter::getDefault () ) divInfo2 = false;
			divStart ( os, "dates", divInfo2 );
				formItemMap [FormItemStartMonth::getName ()]->printHTML ( os );
				formItemMap [FormItemStartYear::getName ()]->printHTML ( os );
				formItemMap [FormItemEndMonth::getName ()]->printHTML ( os );
				formItemMap [FormItemEndYear::getName ()]->printHTML ( os );
			divEnd ( os );
		tableHeaderEnd ( os );
	tableRowEnd ( os );
}
void ProjectDateFilterForm::printJavascriptFunctions ( ostream& os )
{
	os << "function showDateFilterItems ( item ) {" << endl;
	os << "\t" << "var val = getSelectValue(item);" << endl;
	os << "\t" << "if ( val == 'All Projects' ) {" << endl;
	os << "\t\t" << "hidediv ( 'dates' );" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "else {" << endl;
	os << "\t\t" << "showdiv ( 'dates' );" << endl;
	os << "\t" << "}" << endl;
	os << "}" << endl;
}
#endif

class FormItemLowMPlusH : public FormItemText {
public:
	FormItemLowMPlusH ( FormValidatingJavascript* fvj ) :
		FormItemText ( "M+H", "filterman.htm#precursor_m_plus_h", getName (), 6, 6, "600.0", fvj->addPositiveFloatingPointValidator (getName (), "Start M+H Range") ) {}
	static string getName () { return "low_m_plus_h"; }
};

class FormItemHighMPlusH : public FormItemText {
public:
	FormItemHighMPlusH ( FormValidatingJavascript* fvj ) :
		FormItemText ( "to", "", getName (), 6, 6, "3000.0", fvj->addPositiveFloatingPointValidator (getName (), "End M+H Range") ) {}
	static string getName () { return "high_m_plus_h"; }
};

class FormItemFullMPlusHRange : public FormItemCheckbox {
public:
	FormItemFullMPlusHRange () :
		FormItemCheckbox ( "All", "", getName (), true ) {}
	static string getName () { return "full_m_plus_h_range"; }
};

MPlusHForm::MPlusHForm ( FormValidatingJavascript* fvj, const VectorConstParameterListPtr& params ) :
	fvj ( fvj )
{
	create ( params );
}
void MPlusHForm::createItems ()
{
	formItemMap [FormItemLowMPlusH::getName ()]			= new FormItemLowMPlusH (fvj);
	formItemMap [FormItemHighMPlusH::getName ()]		= new FormItemHighMPlusH (fvj);
	formItemMap [FormItemFullMPlusHRange::getName ()]	= new FormItemFullMPlusHRange ();
}
void MPlusHForm::setValues ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];

		formItemMap [FormItemLowMPlusH::getName ()]->setValue ( p );
		formItemMap [FormItemHighMPlusH::getName ()]->setValue ( p );
		formItemMap [FormItemFullMPlusHRange::getName ()]->setValue ( p );
	}
}
void MPlusHForm::printHTML ( ostream& os )
{
	formItemMap [FormItemLowMPlusH::getName ()]->printHTML ( os );
	formItemMap [FormItemHighMPlusH::getName ()]->printHTML ( os );
	formItemMap [FormItemFullMPlusHRange::getName ()]->printHTML ( os );
}
FormItemPreviousAA::FormItemPreviousAA ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Previous AA", "allman.htm#previous_aa", getName (), 2, 2, "1", fvj->addPositiveNonZeroIntegerValidator ( getName (), "Previous AA" ) )
{
}
FormItemNextAA::FormItemNextAA ( FormValidatingJavascript* fvj ) :
	FormItemText ( "Next AA", "allman.htm#next_aa", getName (), 2, 2, "1", fvj->addPositiveNonZeroIntegerValidator ( getName (), "Next AA" ) )
{
}
