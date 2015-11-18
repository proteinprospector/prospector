/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_srch_form.h                                                *
*                                                                             *
*  Created    : December 7th 2004                                             *
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

#ifndef __lu_srch_form_h
#define __lu_srch_form_h

#include <ostream>
#include <lg_string.h>
#include <lgen_define.h>
#include <lu_html_form.h>
#include <lu_pros_form.h>

class FormValidatingJavascript;

class ResourceLinkTable {
protected:
	int numRows;
	int numColumns;
	VectorPairStringString resourceList;
	static PairStringString getPair ( const std::string& form, const std::string& name );
	static PairStringString getPair ( const std::string& form, const std::string& name, const std::vector <PairStringString>& p );
	void getNumRows ();
public:
	void printHTML ( std::ostream& os ) const;
};

class OldResourceLinkTable : public ResourceLinkTable {
	void getResourceList ( const std::string& s = "" );
public:
	OldResourceLinkTable ();
};

class SearchResultsForm : public ProspectorForm {
	bool set;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	SearchResultsForm ( const VectorConstParameterListPtr& params, bool set = false );
	virtual void printHTML ( std::ostream& os );
	void printHTMLFAIndex ( std::ostream& os );
};

class SaveResultsForm : public ProspectorForm {
	FormValidatingJavascript* fvj;
	bool defaultSave;
	bool tabDelimitedTextOption;

	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	SaveResultsForm ( const VectorConstParameterListPtr& params, FormValidatingJavascript* fvj, bool defaultSave = false, bool tabDelimitedTextOption = false );
	virtual void printHTML ( std::ostream& os );
};

class FormItemOutputType : public FormItemSelect {
	static const char* options [];
	static const char* options2 [];
	static const char* options3 [];
public:
	FormItemOutputType ( bool tabDelimitedTextOption = false );
	FormItemOutputType ( bool dummy1, bool dummy2 );
	static std::string getName () { return "output_type"; }
};

class FormItemDNAReadingFrame : public FormItemSelect {
	static const char* options [];
public:
	FormItemDNAReadingFrame ();
	static std::string getName () { return "dna_reading_frame"; }
};

class FormItemTaxonomy : public FormItemSelectMultiple {
public:
	FormItemTaxonomy ();
	static std::string getName () { return "species"; }
};

class FormItemTaxonomyRemove : public FormItemCheckbox {
public:
	FormItemTaxonomyRemove ();
	static std::string getName () { return "species_remove"; }
};

class FormItemTaxonomyNames : public FormItemTextArea {
public:
	FormItemTaxonomyNames ();
	static std::string getName () { return "species_names"; }
};

class ProteinMWForm : public ProspectorForm {
	FormValidatingJavascript* fvj;
	std::string ms;

	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	ProteinMWForm ( FormValidatingJavascript* fvj, const VectorConstParameterListPtr& params, const std::string& ms = "" );
	virtual void printHTML ( std::ostream& os );
};

class ProteinPIForm : public ProspectorForm {
	FormValidatingJavascript* fvj;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	ProteinPIForm ( FormValidatingJavascript* fvj, const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class FormItemNames : public FormItemTextArea {
public:
	FormItemNames ();
	static std::string getName () { return "names"; }
};

class FormItemAccessionNums : public FormItemTextArea {
public:
	FormItemAccessionNums ();
	static std::string getName () { return "accession_nums"; }
};

class FormItemAddAccessionNumbers : public FormItemTextArea {
public:
	FormItemAddAccessionNumbers ();
	static std::string getName () { return "add_accession_numbers"; }
};

class FormItemSaveParameters : public FormItemCheckbox {
public:
	FormItemSaveParameters ();
	static std::string getName () { return "save_params"; }
};

class MSToleranceForm : public ProspectorForm {
	FormValidatingJavascript* fvj;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSToleranceForm ( FormValidatingJavascript* fvj, const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSMSToleranceForm : public ProspectorForm {
	FormValidatingJavascript* fvj;
	bool pChrgItem;
	bool pChrgRangeItem;
	int num;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSMSToleranceForm ( FormValidatingJavascript* fvj, const VectorConstParameterListPtr& params, bool pChrgItem, bool pChrgRangeItem, int num = 1 );
	virtual void printHTML ( std::ostream& os );
};

class FormItemMassSystematicError : public FormItemText {
public:
	FormItemMassSystematicError ( FormValidatingJavascript* fvj, const std::string& ms );
	static std::string getName ( const std::string& ms )
		{ return ms + std::string ("parent_mass_systematic_error"); }
};

class FormItemMSMSPrecursorCharge : public FormItemSelect {
	static const char* options [];
public:
	FormItemMSMSPrecursorCharge ();
	static std::string getName () { return "msms_precursor_charge"; }
};

class FormItemMSMSPrecursorChargeRange : public FormItemSelect {
	static const char* options [];
public:
	FormItemMSMSPrecursorChargeRange ();
	static std::string getName () { return "msms_precursor_charge_range"; }
};

class SearchForm : public ProspectorForm {
	FormValidatingJavascript* fvj;
	std::string searchName;
	bool pattTypeSearch;
	bool protSearch;
	bool fitTypeSearch;
	bool msSearch;
	bool msmsSearch;
	int num;
	SearchResultsForm searchResultsForm;
	SaveResultsForm saveResultsForm;
	ProteinMWForm* proteinMWForm;
	ProteinMWForm* msProteinMWForm;
	ProteinMWForm* msmsProteinMWForm;
	ProteinPIForm proteinPIForm;

	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
	void printHTMLTopRight ( std::ostream& os );
	void printHTMLPreSearch ( std::ostream& os );
public:
	SearchForm ( const VectorConstParameterListPtr& params, const std::string& searchName, FormValidatingJavascript* fvj, int num = 1 );
	virtual void printHTML ( std::ostream& os );
	virtual void printCGI ( std::ostream& os ) const;
	virtual void printHTMLJavascriptHidden ( std::ostream& os ) const;
};

class FormItemSearchName : public FormItemText {
public:
	FormItemSearchName ( const std::string& value );
	static std::string getName () { return "search_name"; }
};

class FormItemReportTitle : public FormItemText {
public:
	FormItemReportTitle ( const std::string& value );
	static std::string getName () { return "report_title"; }
};

class FormItemVersion : public FormItemText {
public:
	FormItemVersion ();
	static std::string getName () { return "version"; }
};

class FormItemDatabase : public FormItemSelect {
public:
	FormItemDatabase ();
	FormItemDatabase ( bool userDefault );
	static std::string getName () { return "database"; }
	void setValue ( const ParameterList* p, const std::string& n = "" );
};

class FormItemMultipleDatabase : public FormItemSelectMultiple {
public:
	FormItemMultipleDatabase ();
	static std::string getName () { return "database"; }
	void setValue ( const ParameterList* p, const std::string& n = "" );
};

class FormItemNTerminusAALimit : public FormItemText {
public:
	FormItemNTerminusAALimit ( FormValidatingJavascript* fvj );
	static std::string getName () { return "n_term_aa_limit"; }
};

class FormItemEnzyme : public FormItemSelect {
public:
	FormItemEnzyme ( const std::string& v = "Trypsin", bool noNoEnzyme = false );
	static std::string getName () { return "enzyme"; }
};

class FormItemAllowNonSpecific : public FormItemSelect {
	static const char* options [];
public:
	FormItemAllowNonSpecific ();
	static std::string getName () { return "allow_non_specific"; }
};

class FormItemMissedCleavages : public FormItemText {
public:
	FormItemMissedCleavages ( FormValidatingJavascript* fvj, const std::string& val = "1" );
	static std::string getName () { return "missed_cleavages"; }
};

class FormItemS : public FormItemCheckbox {
public:
	FormItemS ( int num = 1 );
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "s" + n;
	}
};

class FormItemNTerm : public FormItemSelect {
public:
	FormItemNTerm ( int num = 1, bool showLabel = true );
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "nterm" + n;
	}
};

class FormItemCTerm : public FormItemSelect {
public:
	FormItemCTerm ( int num = 1, bool showLabel = true );
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "cterm" + n;
	}
};

class FormItemConstMod : public FormItemSelectMultiple {
public:
	FormItemConstMod ( bool showTerminalMods, bool nonSelected, int num = 1 );
	FormItemConstMod ();
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "const_mod" + n;
	}
};

class FormItemComment : public FormItemText {
public:
	FormItemComment ();
	static std::string getName () { return "comment"; }
};

class FormItemDisplayGraph : public FormItemCheckbox {
public:
	FormItemDisplayGraph ( bool set = true );
	static std::string getName () { return "display_graph"; }
};

class FormItemDetailedReport : public FormItemCheckbox {
public:
	FormItemDetailedReport ();
	static std::string getName () { return "detailed_report"; }
};

class FormItemModAA : public FormItemSelectMultiple {
	static const char* optionsFit [];
	static const char* optionsDigest [];
public:
	FormItemModAA ( const std::string& type );
	static std::string getName () { return "mod_AA"; }
};

class FormItemUserName : public FormItemSelect {
public:
	FormItemUserName ( const std::string& num );
	static std::string getName ( const std::string& num ) { return std::string ( "user" ) + num + std::string ( "_name" ); }
};

class FormItemMaxReportedHits : public FormItemText {
public:
	FormItemMaxReportedHits ( const std::string& ms, const std::string& value, FormValidatingJavascript* fvj, bool allowBlank );
	static std::string getName ( const std::string& ms = "" ) { return ms + std::string ( "max_reported_hits" ); }
};

class FormItemMinMatches : public FormItemText {
public:
	FormItemMinMatches ( FormValidatingJavascript* fvj );
	static std::string getName () { return "min_matches"; }
};

class FormItemMowseOn : public FormItemCheckbox {
public:
	FormItemMowseOn ();
	static std::string getName () { return "mowse_on"; }
};

class FormItemMowsePfactor : public FormItemText {
public:
	FormItemMowsePfactor ( FormValidatingJavascript* fvj );
	static std::string getName () { return "mowse_pfactor"; }
};

#ifdef CHEM_SCORE
class FormItemChemScore : public FormItemCheckbox {
public:
	FormItemChemScore ();
	static std::string getName () { return "chem_score"; }
};

class FormItemMetOxFactor : public FormItemText {
public:
	FormItemMetOxFactor ( FormValidatingJavascript* fvj );
	static std::string getName () { return "met_ox_factor"; }
};
#endif

class FormItemMaxHits : public FormItemText {
public:
	FormItemMaxHits ( const std::string& val = "9999999" );
	static std::string getName () { return "max_hits"; }
};

class FormItemMSMSMaxModifications : public FormItemText {
public:
	FormItemMSMSMaxModifications ( FormValidatingJavascript* fvj );
	static std::string getName () { return "msms_max_modifications"; }
};

class FormItemMSMSMaxPeptidePermutations : public FormItemText {
public:
	FormItemMSMSMaxPeptidePermutations ( FormValidatingJavascript* fvj );
	static std::string getName () { return "msms_max_peptide_permutations"; }
};

class FormItemModAAList : public FormItemSelectMultiple {
public:
	FormItemModAAList ( int siz );
	FormItemModAAList ( bool msms, bool digest, int siz );
	static std::string getName ( const std::string& s = "" ) { return s + std::string ( "mod_AA" ); }
};

class FormItemParentMassConvert : public FormItemSelect {
	static const char* optionsMS [];
	static const char* optionsMSMS [];
public:
	FormItemParentMassConvert ( const std::string& type, int num = 1 );
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "parent_mass_convert" + n;
	}
};

class FormItemMSParentMassTolerance : public FormItemText {
public:
	FormItemMSParentMassTolerance ();
	static std::string getName () { return "ms_parent_mass_tolerance"; }
};

class FormItemMSParentMassToleranceUnits : public FormItemSelect {
	static const char* options [];
public:
	FormItemMSParentMassToleranceUnits ( const std::string& v = "ppm" );
	static std::string getName () { return "ms_parent_mass_tolerance_units"; }
};

class FormItemMSMSParentMassTolerance : public FormItemText {
public:
	FormItemMSMSParentMassTolerance ( FormValidatingJavascript* fvj, const std::string& value = "200" );
	static std::string getName () { return "msms_parent_mass_tolerance"; }
};

class FormItemMSMSParentMassToleranceUnits : public FormItemSelect {
	static const char* options [];
public:
	FormItemMSMSParentMassToleranceUnits ( const std::string& v = "ppm" );
	static std::string getName () { return "msms_parent_mass_tolerance_units"; }
};

class FormItemFragmentMassesTolerance : public FormItemText {
public:
	FormItemFragmentMassesTolerance ( FormValidatingJavascript* fvj, const std::string& v = "300", int num = 1 );
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "fragment_masses_tolerance" + n;
	}
};

class FormItemFragmentMassesToleranceUnits : public FormItemSelect {
	static const char* options [];
public:
	FormItemFragmentMassesToleranceUnits ( const std::string& v = "ppm", int num = 1 );
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "fragment_masses_tolerance_units" + n;
	}
};

class FormItemMaxSavedTagHits : public FormItemText {
public:
	FormItemMaxSavedTagHits ( FormValidatingJavascript* fvj, const std::string& value = "1000" );
	static std::string getName () { return "max_saved_tag_hits"; }
};

class FormItemUseInstrumentIonTypes : public FormItemCheckbox {
public:
	FormItemUseInstrumentIonTypes ( bool val = true, int num = 1 );
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "use_instrument_ion_types" + n;
	}
};

class FormItemCreateParams : public FormItemCheckbox {
public:
	FormItemCreateParams ();
	static std::string getName () { return "create_params"; }
};

class FormItemScriptFilename : public FormItemText {
public:
	FormItemScriptFilename ();
	static std::string getName () { return "script_filename"; }
};

class FormItemParentContaminantMasses : public FormItemTextArea {
public:
	FormItemParentContaminantMasses ();
	static std::string getName () { return "parent_contaminant_masses"; }
};

class FormItemCompMaskType : public FormItemSelect {
	static const char* options [];
public:
	FormItemCompMaskType ( const std::string& defaultLogic = "AND" );
	static std::string getName () { return "comp_mask_type"; }
};

class FormItemUserProteinSequence : public FormItemTextArea {
public:
	FormItemUserProteinSequence ();
	FormItemUserProteinSequence ( bool databaseSequence, const std::string& sequence = "" );
	static std::string getName () { return "user_protein_sequence"; }
};

class FormItemIonType : public FormItemCheckbox {
public:
	FormItemIonType ( const std::string& label, const std::string& value, bool checked );
	static std::string getName () { return "it"; }
};

class FormItemAlternative : public FormItemCheckbox {
public:
	FormItemAlternative ();
	static std::string getName () { return "alternative"; }
};

class FormItemDiscriminating : public FormItemCheckbox {
public:
	FormItemDiscriminating ();
	static std::string getName () { return "discriminating"; }
};

class FormItemMaxLosses : public FormItemSelect {
	static const char* options [];
public:
	FormItemMaxLosses ();
	static std::string getName () { return "max_losses"; }
};

class FormItemMultiZInternal : public FormItemCheckbox {
public:
	FormItemMultiZInternal ();
	static std::string getName () { return "multi_z_internal"; }
};

class FormItemCalibrate : public FormItemCheckbox {
public:
	FormItemCalibrate ();
	static std::string getName () { return "calibrate"; }
};

class FormItemCalTolerance : public FormItemText {
public:
	FormItemCalTolerance ( const std::string& units, const std::string& val );
	static std::string getName () { return "cal_tolerance"; }
};

class FormItemDataPlotted : public FormItemSelect {
	static const char* options [];
public:
	FormItemDataPlotted ( const std::string& method = "Centroid" );
	static std::string getName () { return "data_plotted"; }
};

class FormItemHideProteinSequence : public FormItemCheckbox {
public:
	FormItemHideProteinSequence ( bool flag );
	static std::string getName () { return "hide_protein_sequence"; }
};

class FormItemExpectationCalculationMethod : public FormItemSelect {
	static const char* options [];
public:
	FormItemExpectationCalculationMethod ( const std::string& method );
	static std::string getName () { return "expect_calc_method"; }
};

class FormItemMSMSPkFilter : public FormItemSelect {
	static const char* options [];
	static const char* options2 [];
public:
	FormItemMSMSPkFilter ( bool allowUnprocessed, int num = 1 );
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "msms_pk_filter" + n;
	}
};

class FormItemMaxMSMSPeaks : public FormItemText {
public:
	FormItemMaxMSMSPeaks ( FormValidatingJavascript* fvj, bool allowBlank, int num = 1 );
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "msms_max_peaks" + n;
	}
};

class FormItemMSMSMinPrecursorMass : public FormItemText {
public:
	FormItemMSMSMinPrecursorMass ( FormValidatingJavascript* fvj, bool allowBlank );
	static std::string getName () { return "msms_min_precursor_mass"; }
};

class FormItemRawType : public FormItemSelect {
	static const char* options [];
public:
	FormItemRawType ( const std::string& rawType = options [0] );
	static std::string getName () { return "raw_type"; }
};

class FormItemIsotopePurity : public FormItemText {
public:
	FormItemIsotopePurity ( const std::string& element, int isotope, FormValidatingJavascript* fvj, double value = 100.0 );
	static std::string getName ( const std::string& element, int isotope ) { return "percent_" + element + gen_itoa ( isotope ); }
};

class FormItemReporterIonWindow : public FormItemText {
public:
	FormItemReporterIonWindow ( FormValidatingJavascript* fvj, const std::string& value = "0.4" );
	static std::string getName () { return "reporter_ion_window"; }
};

class FormItemLinkSearchType : public FormItemSelect {
public:
	FormItemLinkSearchType ( const std::string& val = "", const std::string& chargeFunction = "", int num = 1 );
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "link_search_type" + n;
	}
};

class FormItemMaxLinkMolecules : public FormItemText {
public:
	FormItemMaxLinkMolecules ( FormValidatingJavascript* fvj );
	static std::string getName () { return "max_link_molecules"; }
};

class FormItemLinkAA : public FormItemSelect {
public:
	FormItemLinkAA ( const std::string& chargeFunction = "", int num = 1 );
	static std::string getName ( int num = 1 )
	{
		std::string num1 = ( num == 1 ) ? "" : gen_itoa ( num );
		return "link_aa" + num1;
	}
};

class FormItemComposition : public FormItemText {
public:
	FormItemComposition ( const std::string& label, const std::string& n, int num, const std::string& v, FormValidatingJavascript* fvj );
	static std::string getName ( const std::string& n, int num = 1 )
	{
		std::string num1 = ( num == 1 ) ? "" : gen_itoa ( num );
		return n + std::string ( "_composition" ) + num1;
	}
};

class FormItemAccurateMass : public FormItemText {
public:
	FormItemAccurateMass ( const std::string& n, int num, FormValidatingJavascript* fvj );
	static std::string getName ( const std::string& n, int num = 1 )
	{
		std::string num1 = ( num == 1 ) ? "" : gen_itoa ( num );
		return n + std::string ( "_accurate_mass" ) + num1;
	}
};

class FormItemLabel : public FormItemText {
public:
	FormItemLabel ( const std::string& label, const std::string& n, int num, const std::string& v );
	static std::string getName ( const std::string& n, int num = 1 )
	{
		std::string num1 = ( num == 1 ) ? "" : gen_itoa ( num );
		return n + std::string ( "_label" ) + num1;
	}
};

class FormItemAAModified : public FormItemSelect {
public:
	FormItemAAModified ( const std::string& n, const std::string& v );
	static std::string getName ( const std::string& n ) { return "aa_modified_" + n; }
};

class FormItemSpecificity : public FormItemSelect {
	static const char* modOptions [];
public:
	FormItemSpecificity ( const std::string& n, int num );
	static std::string getName ( const std::string& n, int num = 1 )
	{
		std::string num1 = ( num == 1 ) ? "" : gen_itoa ( num );
		return n + std::string ( "_specificity" ) + num1;
	}
};

class FormItemModFile : public FormItemSelect {
	static const char* modFileOptions [];
public:
	FormItemModFile ();
	static std::string getName () { return "mod_AA_file"; }
};

class FormItemMotifOffset : public FormItemSelect {
	static const char* offsetOptions [];
public:
	FormItemMotifOffset ( const std::string& n = "", int num = 1 );
	static std::string getName ( const std::string& n = "", int num = 1 )
	{
		if ( n.empty () ) return "motif_offset";
		else {
			std::string num1 = ( num == 1 ) ? "" : gen_itoa ( num );
			return n + std::string ( "_motif_offset" ) + num1;
		}
	}
};

class FormItemMotif : public FormItemText {
public:
	FormItemMotif ( FormValidatingJavascript* fvj,  const std::string& n = "", int num = 1 );
	static std::string getName ( const std::string& n = "", int num = 1 ) {
		if ( n.empty () ) return "motif";
		else {
			std::string num1 = ( num == 1 ) ? "" : gen_itoa ( num );
			return n + std::string ( "_motif" ) + num1;
		}
	}
};

class FormItemLimit : public FormItemSelect {
	static const char* limitOptions [];
public:
	FormItemLimit ( const std::string& n, int num );
	static std::string getName ( const std::string& n, int num = 1 )
	{
		std::string num1 = ( num == 1 ) ? "" : gen_itoa ( num );
		return n + std::string ( "_limit" ) + num1;
	}
	void setOptions ( const ParameterList* p, const std::string& n = "" );
};

class UserAAForm : public ProspectorForm {
	FormValidatingJavascript* fvj;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	UserAAForm ( const VectorConstParameterListPtr& params, FormValidatingJavascript* fvj );
	virtual void printHTML ( std::ostream& os );
};

class CompIonForm : public ProspectorForm {
protected:
	std::string aaList;
	virtual void setValues ( const VectorConstParameterListPtr& params );
	virtual void createItems ();
	virtual std::string getName () const = 0;
public:
	CompIonForm ( const VectorConstParameterListPtr& params, const std::string& aaList );
	virtual void printHTML ( std::ostream& os ) = 0;
};

class PresentCompIonForm : public CompIonForm {
	std::string prefix;
	std::string getName () const;
	void setValues ( const VectorConstParameterListPtr& params );
public:
	PresentCompIonForm ( const VectorConstParameterListPtr& params, const std::string& aaList, const std::string& prefix = "" );
	virtual void printHTML ( std::ostream& os );
	void printHTML2Rows ( std::ostream& os );
};

class MSCompExtraIonForm : public CompIonForm {
	std::string getName () const;
public:
	MSCompExtraIonForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSCompAbsentIonForm : public CompIonForm {
	std::string getName () const;
public:
	MSCompAbsentIonForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSCompModifiedIonForm : public CompIonForm {
	std::string getName () const;
public:
	MSCompModifiedIonForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class ImmoniumIonForm : public ProspectorForm {
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
	std::string getName () const;
public:
	ImmoniumIonForm ( const VectorConstParameterListPtr& params );
	void setMSCompImmoniumDefaults ();
	virtual void printHTML ( std::ostream& os );
};

class HeaderForm : public ProspectorForm {
	std::string searchName;
	std::string reportTitle;
	void setValues ( const VectorConstParameterListPtr& params ) {}
	void createItems ();
public:
	HeaderForm ( const VectorConstParameterListPtr& params, const std::string& searchName, const std::string& reportTitle );
	virtual void printHTML ( std::ostream& os );
};

void msBridgeJavascriptFunctions ( std::ostream& os );
void massModCrosslinkingJavascriptFunctions ( std::ostream& os, int n = 1 );

class MassModificationForm : public ProspectorForm {
	static const int numMotifs;
	PresentCompIonForm compIonForm;
	FormValidatingJavascript* fvj;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
	void printSetModRangeTypeDefaults ( std::ostream& os );
	void setMMVisualizationFlags ( const std::string& val, bool& div_mm_1, bool& div_mm_2, bool& div_mm_3, bool& div_mm_4 ) const;
public:
	MassModificationForm ( const VectorConstParameterListPtr& params, FormValidatingJavascript* fvj );
	virtual void printHTML ( std::ostream& os );
	virtual void printCGI ( std::ostream& os ) const;
	virtual void printHTMLJavascriptHidden ( std::ostream& os ) const;
	bool getDivCL ();
	static int getNumMotifs () { return numMotifs; }
};

class ExpandableJavascriptBlock;

class VariableModsSettingJavascript;

class VariableModsForm : public ProspectorForm {
	std::string prefix;
	VariableModsSettingJavascript* vmsj;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	VariableModsForm ( const VectorConstParameterListPtr& params, const std::string& suffix, FormValidatingJavascript* fvj );
	~VariableModsForm ();
	virtual void printHTML ( std::ostream& os );
};

class ExtraUserModsForm : public ProspectorForm {
	static int numUserMods;
	FormValidatingJavascript* fvj;
	ExpandableJavascriptBlock* ejb;
	bool massItem;
	int num;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	ExtraUserModsForm ( const VectorConstParameterListPtr& params, FormValidatingJavascript* fvj, bool massItem, int num = 1 );
	virtual void printHTML ( std::ostream& os );
	static int getNumUserMods () { return numUserMods; }
	void replaceSubstrings ( std::string& s ) const;
	static void copyToCGI ( std::ostream& os, const ParameterList* params );
};

class CrosslinkingForm : public ProspectorForm {
	static int numUserMods;
	FormValidatingJavascript* fvj;
	bool tagForm;
	ExpandableJavascriptBlock* ejb;
	int num;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
	void printHTMLUDLP ( std::ostream& os );
public:
	CrosslinkingForm ( const VectorConstParameterListPtr& params, FormValidatingJavascript* fvj, bool tagForm = false, int num = 1 );
	virtual void printHTML ( std::ostream& os );
	static int getNumUserMods () { return numUserMods; }
};

class MatrixModeForm : public ProspectorForm {
	std::string prefix;
	StringVector names;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
	static std::string getName ( const std::string& ms );
public:
	MatrixModeForm ( const VectorConstParameterListPtr& params, const std::string& prefix );
	void printHTML ( std::ostream& os );
};

class FormItemUploadPeakList : public FormItemFile {
public:
	FormItemUploadPeakList () :
		FormItemFile ( "", "", getName (), 100, 256, "" ) {}
	static std::string getName () { return "upload_temp_peak_list"; }	// temp in the name indicates that the upload is written to the temp directory
};

#ifdef MYSQL_DATABASE
class ProjectDateFilter;

class ProjectDateFilterForm : public ProspectorForm {
	ProjectDateFilter* pdf;
public:
	ProjectDateFilterForm ( const ParameterList* paramList );
	~ProjectDateFilterForm ();
	void createItems ();
	void setValues ( const VectorConstParameterListPtr& params );
	void printHTML ( std::ostream& os ) {}
	void printHTML ( std::ostream& os, bool showProjectsDiv );
	void printJavascriptFunctions ( std::ostream& os );
	ProjectDateFilter* getProjectDateFilter () const { return pdf; }
};
#endif

class MPlusHForm : public ProspectorForm {
	FormValidatingJavascript* fvj;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MPlusHForm ( FormValidatingJavascript* fvj, const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class FormItemPreviousAA : public FormItemText {
public:
	FormItemPreviousAA ( FormValidatingJavascript* fvj );
	static std::string getName () { return "previous_aa"; }
};

class FormItemNextAA : public FormItemText {
public:
	FormItemNextAA ( FormValidatingJavascript* fvj );
	static std::string getName () { return "next_aa"; }
};

#endif /* ! __lu_srch_form_h */
