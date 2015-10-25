/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prod_par.cpp                                               *
*                                                                             *
*  Created    : June 18th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#include <lg_string.h>
#include <lgen_error.h>
#include <lu_aa_info.h>
#include <lu_ambiguity.h>
#include <lu_prod_par.h>
#include <lu_usermod.h>
#include <lu_cgi_val.h>
#include <lu_get_link.h>
#include <lu_mut_mtrx.h>
#include <lu_iso_par.h>
#include <lu_mass_seq.h>
#include <lu_mass_frag.h>
#include <lu_getfil.h>
#include <lu_param_list.h>
#include <lu_const_mod.h>
#include <lu_inst.h>
#include <lu_srch_form.h>
using std::showpos;
using std::noshowpos;
using std::string;
using std::ostream;
using std::endl;
using std::stringstream;
using std::vector;
using std::make_pair;
using std::pair;

static void printHTMLSeqHighlightAA ( ostream& os, const string& sequence, char aa, bool javascript )
{
	for ( StringSizeType i = 0 ; i < sequence.length () ; i++ ) {
		if ( sequence [i] == aa ) {
			if ( javascript )
				os << "<font size=\\\"+1\\\"><b>";
			else
				os << "<font size=\"+1\"><b>";
		}
		os << sequence [i];
		if ( sequence [i] == aa ) os << "</b></font>";
	}
}
string MSProductParameters::MAX_CHARGE = "max_charge";

string MSProductParameters::MAX_CHARGE_DEFAULT = "No Limit";

MSProductParameters::MSProductParameters ( ParameterList* params ) :
	MSProgramParameters	( params ),

	specID				( params ),
	dataSetInfo			( params, true ),
	parentMassTolerance	( "msms_parent_mass", params ),
	productMassTolerance( "fragment_masses", params ),

	maxCharge			( initMaxCharge ( params->getStringValue	( MAX_CHARGE, MAX_CHARGE_DEFAULT ) ) ),
	calibrate			( params->getBoolValue	( "calibrate", false ) ),
	calTolerance		( params->getDoubleValue( "cal_tolerance" ) ),
	dataPlotted			( params->getStringValue( "data_plotted", "Centroid" ) ),
	alternative			( params->getBoolValue	( "alternative", false ) ),
	discriminating		( params->getBoolValue	( "discriminating", false ) ),
	linkInfo			( new LinkInfo ( params ) ),
	instrumentName		( params->getStringValue ( "instrument_name", "" ) ),
	biemannParams		( 0 )
{
	ConstMod::setProductFlag ();
	Usermod::setProductFlag ();

	Usermod::initialiseAllUsermodAAInfo ();						// These are independent of the sequence
	ExtraUserMods::instance ().addUserMods ( params );
	ElementalFormulaInfo e1 ( "user_aa", params );
	userAAElemForms.push_back ( e1.getElementalFormula () );
	ElementalFormulaInfo e2 ( "user_aa_2", params );
	userAAElemForms.push_back ( e2.getElementalFormula () );
	ElementalFormulaInfo e3 ( "user_aa_3", params );
	userAAElemForms.push_back ( e3.getElementalFormula () );
	ElementalFormulaInfo e4 ( "user_aa_4", params );
	userAAElemForms.push_back ( e4.getElementalFormula () );
	initialise_amino_acid_weights ( MapStringConstModPtr (), userAAElemForms, true );
	initialiseInstrumentName ( instrumentName );
	massInfo = new MassInfo ( params );

	string mods = params->getStringValue ( "mods" );
	if ( !mods.empty () ) {
		params->removeName ( "mods" );
		ModificationAmbiguity ma ( mods );
		string startS, startNT, startCT, startNL;
		getSequenceElements ( params, "", startS, startNT, startCT, startNL );
		for ( int i = 0 ; ; i++ ) {
			string num = ( i == 0 ) ? "" : gen_itoa ( i+1 );
			string s = startS;
			string nt = startNT;
			string ct = startCT;
			string nl = startNL;
			sFlag.push_back ( i < 6 );
			constMods.resize ( i+1 );
			bool flag = ma.getNextSequence ( s, nt, ct, nl );
			params->addOrReplaceName ( "sequence" + num, s );
			params->addOrReplaceName ( "nterm" + num, nt );
			params->addOrReplaceName ( "cterm" + num, ct );
			params->addOrReplaceName ( "nloss" + num, nl );
			params->addOrReplaceName ( "s" + num, ( sFlag.back () == true ) ? "1" : "0" );
			nTermName.push_back ( nt );
			cTermName.push_back ( ct );
			nLossName.push_back ( nl );
			if ( !nTermName [i].empty () ) constMods [i]["n"] = new ConstMod ( nTermName [i] + " (N-term)" );
			if ( !cTermName [i].empty () ) constMods [i]["c"] = new ConstMod ( cTermName [i] + " (C-term)" );
			sequence.push_back ( initSequence ( s ) );
			if ( i == 0 ) biemannParams = new BiemannParameters ( params );
			if ( !flag ) break;
		}
	}
	else {
		for ( int i = 0 ; ; i++ ) {
			string num = ( i == 0 ) ? "" : gen_itoa ( i+1 );
			string s, nT, cT, nL;
			if ( !getSequenceElements ( params, num, s, nT, cT, nL ) ) break;
			nTermName.push_back ( nT );
			cTermName.push_back ( cT );
			nLossName.push_back ( nL );
			sFlag.push_back ( params->getBoolValue ( "s" + num, false ) );
			constMods.resize ( i+1 );
			if ( !nTermName [i].empty () ) constMods [i]["n"] = new ConstMod ( nTermName [i] + " (N-term)" );
			if ( !cTermName [i].empty () ) constMods [i]["c"] = new ConstMod ( cTermName [i] + " (C-term)" );
			sequence.push_back ( initSequence ( s ) );
			if ( i == 0 ) biemannParams = new BiemannParameters ( params );
		}
	}
	msmsPeakFilterOptions = new MSMSPeakFilterOptions ( params );
}
MSProductParameters::~MSProductParameters ()
{
	delete msmsPeakFilterOptions;
	delete biemannParams;
	delete massInfo;
}
bool MSProductParameters::getSequenceElements ( ParameterList* params, const string& num, string& sequence, string& nTermName, string& cTermName, string& neutralLossName )
{
	sequence = params->getStrippedStringValue ( "sequence" + num, "" );
	if ( sequence.empty () ) return false;
	nTermName = params->getStringValue ( "nterm" + num, "" );
	cTermName = params->getStringValue ( "cterm" + num, "" );
	neutralLossName = params->getStringValue ( "nloss" + num, "" );
	if ( nTermName.empty () && cTermName.empty () && neutralLossName.empty () ) {
		static StringVector nTermNames = ConstMod::getNTermNames ();
		static StringVector cTermNames = ConstMod::getCTermNames ();
		static StringVector nLossNames = ConstMod::getNLossNames ();
		IntVector plusIdx;
		IntVector minusIdx;
		for ( StringSizeType i = 0 ; i < sequence.length () ; i++ ) {
			char p = sequence [i];
			if ( p == '+' )		plusIdx.push_back ( i );
			else if ( p == '-' )minusIdx.push_back ( i );
			else if ( p == '(' ) {
				int bracket = 0;
				for ( ; i < sequence.length () ; i++ ) {
					char a = sequence [i];
					if ( a == '(' ) bracket++;
					if ( a == ')' ) bracket--;
					if ( bracket == 0 ) break;
				}
			}
		}
		if ( !plusIdx.empty () ) {
			for ( int i = 0 ; i < nLossNames.size () ; i++ ) {
				string nL = nLossNames [i];
				if ( isSuffix ( sequence, "+" + nL ) ) {
					sequence = sequence.substr ( 0, sequence.length () - 1 - nL.length () );
					neutralLossName = nL;
					break;
				}
			}
			if ( neutralLossName.empty () ) {
				string possNLoss = sequence.substr ( plusIdx.back () + 1 );	// Add 1 to the last +
				if ( genStringIsFloat ( possNLoss ) ) {						// Check if number
					sequence = sequence.substr ( 0, plusIdx.back () );
					neutralLossName = possNLoss;
				}
			}
			while ( !minusIdx.empty () && minusIdx.back () >= sequence.length () ) minusIdx.pop_back ();
		}
		if ( !minusIdx.empty () ) {
			for ( int i = 0 ; i < cTermNames.size () ; i++ ) {
				string nC = cTermNames [i];
				if ( nC.empty () ) continue;
				if ( isSuffix ( sequence, "-" + nC ) ) {
					sequence = sequence.substr ( 0, sequence.length () - 1 - nC.length () );
					cTermName = nC;
					break;
				}
			}
			if ( cTermName.empty () ) {
				int idx = minusIdx.back ();
				if ( sequence [idx-1] != '-' ) idx += 1;
				string possCTerm = sequence.substr ( idx );	// Add 1 to the last -
				if ( genStringIsFloat ( possCTerm ) ) {							// Check if number
					sequence = sequence.substr ( 0, idx - 1 );
					cTermName = possCTerm;
				}
			}
			while ( minusIdx.back () >= sequence.length () ) minusIdx.pop_back ();
		}
		if ( !minusIdx.empty () ) {
			for ( int j = 0 ; j < nTermNames.size () ; j++ ) {
				string nN = nTermNames [j];
				if ( nN.empty () ) continue;
				if ( isPrefix ( sequence, nN + "-" ) ) {
					sequence = sequence.substr ( nN.length () + 1 );
					nTermName = nN;
					break;
				}
			}
			if ( nTermName.empty () ) {
				int idx = ( sequence [0] == '-' ) ? minusIdx [1] : minusIdx [0];
				string possNTerm = sequence.substr ( 0, idx );
				if ( genStringIsFloat ( possNTerm ) ) {							// Check if number
					sequence = sequence.substr ( idx + 1 );
					nTermName = possNTerm;
				}
			}
		}
		params->addOrReplaceName ( "sequence" + num, sequence );
		params->addOrReplaceName ( "nterm" + num, nTermName );
		params->addOrReplaceName ( "cterm" + num, cTermName );
		params->addOrReplaceName ( "nloss" + num, neutralLossName );
	}
	return true;
}
int MSProductParameters::initMaxCharge ( const string& maxChargeString )
{
	return ( maxChargeString == "No Limit" ) ? 0 : atoi ( maxChargeString.c_str () );
}
void MSProductParameters::printHTML ( ostream& os )
{
	if ( !sequence.empty () && !sequence [0].empty () ) {

		int charge = 1;
		if ( dataSetInfo.getNumDataSets () > 0 ) {
			charge = dataSetInfo.getDataSet ( 0 )->getPrecursorCharge ();
		}
		int num = 1;
		for ( StringVectorVectorSizeType i = 0 ; i < sequence.size () ; i++ ) {
			if ( sFlag [i] ) {
				int colour = sequence.size () > 1 ? num : 0;
				if ( num > 6 ) {
					ErrorHandler::genError ()->error ( "A maximum of 6 sequences can be displayed at one time.\n" );
				}
				printSequence ( os, sequence [i], nTermName [i], cTermName [i], nLossName [i], charge, colour );
				num++;
			}
		}
		for ( ElementalFormulaInfoVectorSizeType j = 0 ; j < userAAElemForms.size () ; j++ ) {
			string f = userAAElemForms [j].getFormula ();
			if ( !f.empty () ) {
				os << "User AA Formula " << j+1 << ": <b>" << f << "</b><br />" << endl;
			}
		}
	}
	else
		os << "<div class=\"msprod_sequence\">No sequence specified.</div>" << endl;
}
void MSProductParameters::printSequence ( ostream& os, const StringVector& seq, const string& nName, const string& cName, const string& lName, int charge, int colour ) const
{
	static GenNameValueStream nvs ( MsparamsDir::instance ().getParamPath ( "unimod.txt" ) );
	static string mainURL = nvs.getStringValue ( "main_url" );
	static string startRange = nvs.getStringValue ( "start_range" );
	static string endRange = nvs.getStringValue ( "end_range" );

	os << "<div class=\"msprod_sequence";
	if ( colour != 0 ) os << "_" << colour;
	os << "\">" << endl;

		if ( !nName.empty () ) {
			printPSILink ( os, nName, true, mainURL, startRange, endRange );
			os << '-';
		}
		for ( StringVectorSizeType i = 0 ; i < seq.size () ; i++ ) {
			string s = seq [i];
			if ( s.length () == 1 ) os << s;
			else printPSILink ( os, s, false, mainURL, startRange, endRange );
		}
		if ( !cName.empty () ) {
			os << '-';
			printPSILink ( os, cName, true, mainURL, startRange, endRange );
		}
		if ( !lName.empty () ) {
			os << '+';
			printPSILink ( os, lName, true, mainURL, startRange, endRange );
		}
		if ( charge != 1 ) {
			os << "<sup>" << showpos << charge << noshowpos << "</sup>";
		}
		os << endl;

	os << "</div>" << endl;
}
void MSProductParameters::printPSILink ( ostream& os, const string& s, bool term, const string& mainURL, const string& startRange, const string& endRange )
{
	int start = term ? 0 : 2;
	int len = term ? s.length () : s.length () - 3;
	if ( genIsNumberStart ( s [start] ) ) {
		string massStr = s.substr ( start, len );
		double mass = atof ( massStr.c_str () );
		double startMass = mass - 1.0;
		double endMass = mass + 1.0;
		if ( !term ) os << s [0] << s [1];
		os << "<a";
		os << " target=\"_blank\"";
		os << " href=\"";
		os << mainURL;
		os << startRange;
		genPrint ( os, startMass, 4 );
		os << "&";
		os << endRange;
		genPrint ( os, endMass, 4 );
		os << "\"";
		os << ">";
		os << massStr;
		os << "</a>";
		if ( !term ) os << ')';
	}
	else os << s;
}
int MSProductLink::num = 0;
MSProductLink::MSProductLink ( const string& programName, const MSMSDataSetInfo* dataSetInfo, int dataSet ) :
	dataSetInfo ( dataSetInfo ),
	dataSet ( dataSet ),
	pTol ( 0 )			// Currently used as an indicator for if called from Search Compare
{
}
MSProductLink::MSProductLink ( const string& instrumentName, const string& searchKey, const SpecID& specID, const Tolerance* pTol, const Tolerance* tol ) :
	instrumentName ( instrumentName ),
	searchKey ( searchKey ),
	specID ( specID ),
	pTol ( pTol ),
	tol ( tol )
{
}
MSProductLink::MSProductLink ( const string& instrumentName, const string& dataFilename, const Tolerance* pTol, const Tolerance* tol ) :
	instrumentName ( instrumentName ),
	dataFilename ( dataFilename ),
	pTol ( pTol ),
	tol ( tol )
{
	num++;
	index = num;
}
void MSProductLink::start ( ostream& os, bool url ) const
{
	if ( url ) {
		os << ParameterList::getServer () + ProgramLink::getURLStart ( "mssearch" );
		os << "?";
		putCGI ( os );
	}
	else
		ProgramLink::openLink ( os, getLinkName (), -1 );
}
void MSProductLink::end ( ostream& os, bool url, const string& sequence, const string& sequence2 ) const
{
	if ( !url ) {
		os << "\\\">";
		os << sequence;
		if ( !sequence2.empty () ) {
			os << "<br />" << sequence2;
		}
		ProgramLink::closeLink ( os );
	}
}
void MSProductLink::putCGI ( ostream& os ) const
{
	const ParameterList* params = ProgramLink::getParams ();
	printCGIString ( os, "search_name", "msproduct" );
	printCGIString ( os, "output_type", "HTML" );
	printCGIString ( os, "report_title", "MS-Product" );
	printCGIString ( os, "version", Version::instance ().getVersion () );

	if ( pTol ) {
		if ( !searchKey.empty () ) {								// From Search Compare
			printCGIString ( os, "data_source", "List of Files" );
			printCGIString ( os, "search_key", searchKey );
			specID.putCGI ( os );
			MassInfo massInfo ( params );
			massInfo.putCGIFragment ( os );
			printCGI ( os, "use_instrument_ion_types", true );
		}
		else {														// MS-Viewer
			printCGIString ( os, "data_source", "Data From File" );
			printCGIString ( os, "data_filename", dataFilename );
			params->copyToCGI ( os, "use_instrument_ion_types" );
			if ( params->getBoolValue ( "use_instrument_ion_types" ) == false ) {
				params->copyToCGI ( os, "it" );
			}
			printCGI ( os, "msms_min_precursor_mass", 0.0 );
		}
		printCGIString ( os, "instrument_name", instrumentName );
		printCGI ( os, "display_graph", true );
		pTol->putCGI ( os, "msms_parent_mass" );
		tol->putCGI ( os, "fragment_masses" );
	}
	else {
		if ( dataSetInfo ) {	// MS-Tag
			dataSetInfo->putProductCGI ( os, dataSet );
			printCGI ( os, "display_graph", true );
			ToleranceInfo precursorTolerance ( "msms_parent_mass", params );
			precursorTolerance.putCGI ( os );
			ToleranceInfo productMassTolerance ( "fragment_masses", params );
			productMassTolerance.putCGI ( os );
			params->copyToCGI ( os, "instrument_name" );
			params->copyToCGI ( os, "msms_precursor_charge" );
			if ( params->getStringValue ( "search_name" ) == "msseq" ) {
				printCGI ( os, "use_instrument_ion_types", true );
			}
			else {
				params->copyToCGI ( os, "use_instrument_ion_types" );
				if ( params->getBoolValue ( "use_instrument_ion_types" ) == false ) {
					params->copyToCGI ( os, "it" );
				}
			}
		}
		else {					// MS-Fit/MS-Digest
			params->copyToCGI ( os, "instrument_name" );
			printCGI ( os, "use_instrument_ion_types", true );
		}
		MassInfo massInfo ( params );
		massInfo.putCGIFragment ( os );
	}
	params->copyToCGI ( os, "user_aa_composition" );
	params->copyToCGI ( os, "user_aa_2_composition" );
	params->copyToCGI ( os, "user_aa_3_composition" );
	params->copyToCGI ( os, "user_aa_4_composition" );
	ExtraUserMods::instance ().copyToCGI ( os );
	MSMSPeakFilterOptions::copyToCGI ( os, params );
	IsotopePurity::copyToCGI ( os, params );
}
void MSProductLink::write1 ( ostream& os, const PeptideSequence& ps, bool hideLinks, const int maxCharge, char aa ) const
{
	if ( hideLinks ) {
		printHTMLSeqHighlightAA ( os, ps.getPeptideString (), aa, false );
	}
	else {
		ProgramLink::openLink ( os, getLinkName (), -1 );
		if ( maxCharge == 0 )printCGIString ( os, "max_charge", "No Limit" );
		else				 printCGI ( os, "max_charge", maxCharge );
		printCGI ( os, "s", 1 );
		ps.putPeptideStringCGI ( os );
		os << "\\\">";
		printHTMLSeqHighlightAA ( os, ps.getPeptideString (), aa, true );
		ProgramLink::closeLink ( os );
	}
}
void MSProductLink::write2 ( ostream& os, const string& sequence, const PeptideSequence& ps, bool hideLinks, const int maxCharge, double err ) const
{
	if ( hideLinks ) {
		ps.printHTML ( os, sequence, err );
	}
	else {
		ProgramLink::openLink ( os, getLinkName (), -1 );
		if ( maxCharge == 0 )printCGIString ( os, "max_charge", "No Limit" );
		else				 printCGI ( os, "max_charge", maxCharge );
		printCGI ( os, "s", 1 );
		ps.putCGI ( os, err );
		os << "\\\">";
		ps.printHTML ( os, sequence, err );
		ProgramLink::closeLink ( os );
	}
}
void MSProductLink::write2 ( ostream& os, const CharVector& previousAA, const StringVector& sequence, const CharVector& nextAA, const vector <PeptideSequence>& ps, bool hideLinks, const int maxCharge, const LinkInfo* linkInfo ) const
{
	if ( !hideLinks ) {
		ProgramLink::openLink ( os, getLinkName (), -1 );
		if ( maxCharge == 0 )printCGIString ( os, "max_charge", "No Limit" );
		else				 printCGI ( os, "max_charge", maxCharge );
		printCGIString ( os, "count_pos_z", "Ignore Basic AA" );
		string linkSearchType = linkInfo->getName ();
		printCGIString ( os, "link_search_type", linkSearchType );
		if ( linkSearchType == "User Defined Link" ) {
			printCGIString ( os, "bridge_composition", linkInfo->getBridgeFormula () );
		}
		for ( int i = 0 ; i < sequence.size () ; i++ ) {
			printCGI ( os, "s" + ( i == 0 ? "" : gen_itoa ( i+1 ) ), true );
			ps [i].putCGI ( os, 0.0, i+1 );
		}
		os << "\\\">";
	}
	for ( int j = 0 ; j < sequence.size () ; j++ ) {
		os << "(" << previousAA [j] << ")";
		ps [j].printHTML ( os, sequence [j], 0.0 );
		os << "(" << nextAA [j] << ")";
		if ( j != sequence.size () - 1 ) os << "<br />";
	}
	if ( !hideLinks ) {
		ProgramLink::closeLink ( os );
	}
}
void MSProductLink::write3 ( ostream& os, const string& label, int maxCharge ) const
{
	ProgramLink::openLink ( os, getLinkName (), -1 );
	printCGI ( os, "max_charge", maxCharge );
	end ( os, false, label );
}
void MSProductLink::write4 ( ostream& os, const string& sequence, const string& mods, int maxCharge ) const
{
	ProgramLink::openLink ( os, getLinkName (), -1 );
	if ( maxCharge == 0 )printCGIString ( os, "max_charge", "No Limit" );
	else				 printCGI ( os, "max_charge", maxCharge );
	printCGIString ( os, "sequence", sequence );
	printCGI ( os, "s", 1 );
	if ( !mods.empty () )	printCGIString ( os, "mods", mods );
	printCGI ( os, "alternative", true );
	end ( os, false, sequence );
}
void MSProductLink::write4 ( ostream& os, const string& dbPeptide, const string& sequence, const string& nTermName, const string& cTermName, const string& neutralLossName, const string& mods, int maxCharge ) const
{
	ProgramLink::openLink ( os, getLinkName (), -1 );
	if ( maxCharge == 0 )printCGIString ( os, "max_charge", "No Limit" );
	else				 printCGI ( os, "max_charge", maxCharge );
	printCGIString ( os, "sequence", dbPeptide );
	printCGI ( os, "s", 1 );
	if ( !mods.empty () )	printCGIString ( os, "mods", mods );
	printCGI ( os, "alternative", true );
	os << "\\\">";
	if ( !nTermName.empty () ) os << nTermName << '-';
	os << sequence;
	if ( !cTermName.empty () ) os << '-' << cTermName;
	if ( !neutralLossName.empty () ) os << '+' << neutralLossName;
	ProgramLink::closeLink ( os );
}
void MSProductLink::write4 ( ostream& os, const StringVector& previousAA, const StringVector& sequence, const StringVector& nTermName, const StringVector& cTermName, const StringVector& neutralLossName, const StringVector& nextAA, bool hideLinks, const int maxCharge, const LinkInfo* linkInfo ) const
{
	if ( !hideLinks ) {
		ProgramLink::openLink ( os, getLinkName (), -1 );
		if ( maxCharge == 0 )printCGIString ( os, "max_charge", "No Limit" );
		else				 printCGI ( os, "max_charge", maxCharge );
		printCGIString ( os, "count_pos_z", "Ignore Basic AA" );
		string linkSearchType = linkInfo->getName ();
		printCGIString ( os, "link_search_type", linkSearchType );
		if ( linkSearchType == "User Defined Link" ) {
			printCGIString ( os, "bridge_composition", linkInfo->getBridgeFormula () );
		}
		for ( int i = 0 ; i < sequence.size () ; i++ ) {
			int n = i+1;
			string num = ( n == 1 ) ? "" : gen_itoa ( n ); 
			printCGI ( os, "s" + num, true );
			printCGIString ( os, "sequence" + num, sequence [i] );
			if ( !nTermName [i].empty () )		printCGIString ( os, "nterm" + num, nTermName [i] );
			if ( !cTermName [i].empty () )		printCGIString ( os, "cterm" + num, cTermName [i] );
			if ( !neutralLossName [i].empty () )printCGIString ( os, "nloss" + num, neutralLossName [i] );
		}
		os << "\\\">";
	}
	for ( int j = 0 ; j < sequence.size () ; j++ ) {
		if ( !previousAA [j].empty () ) os << "(" << previousAA [j] << ")";
		if ( !nTermName [j].empty () ) os << nTermName [j] << '-';
		os << sequence [j];
		if ( !cTermName [j].empty () ) os << '-' << cTermName [j];
		if ( !neutralLossName [j].empty () ) os << '+' << neutralLossName [j];
		if ( !nextAA [j].empty () ) os << "(" << nextAA [j] << ")";
		if ( j != sequence.size () - 1 ) os << "<br />";
	}
	if ( !hideLinks ) {
		ProgramLink::closeLink ( os );
	}
}
void MSProductLink::write5 ( ostream& os, const SpecID& spID, const string& sequence, const string& mods, int maxCharge, bool url ) const
{
	start ( os, url );

	spID.putCGI ( os );
	if ( maxCharge == 0 )printCGIString ( os, "max_charge", "No Limit" );
	else				 printCGI ( os, "max_charge", maxCharge );
	printCGI ( os, "msms_precursor_charge", maxCharge );
	printCGIString ( os, "sequence", sequence );
	printCGI ( os, "s", 1 );
	if ( !mods.empty () )	printCGIString ( os, "mods", mods );
	printCGI ( os, "alternative", true );
	end ( os, url, sequence );
}
void MSProductLink::write5 ( ostream& os, const SpecID& spID, const string& sequence1, const string& sequence2, int maxCharge, bool url, const string& linkSearchType, const string& bridgeComposition ) const
{
	start ( os, url );
	spID.putCGI ( os );
	if ( maxCharge == 0 )printCGIString ( os, "max_charge", "No Limit" );
	else				 printCGI ( os, "max_charge", maxCharge );
	printCGI ( os, "msms_precursor_charge", maxCharge );
	printCGIString ( os, "sequence", sequence1 );
	printCGI ( os, "s", 1 );
	printCGIString ( os, "sequence2", sequence2 );
	printCGI ( os, "s2", 1 );
	printCGIString ( os, "count_pos_z", "Ignore Basic AA" );
	printCGIString ( os, "link_search_type", linkSearchType );
	if ( linkSearchType == "User Defined Link" ) {
		printCGIString ( os, "bridge_composition", bridgeComposition );
	}
	end ( os, url, sequence1, sequence2 );
}
void MSProductLink::write5 ( ostream& os, const string& title, const string& sequence, const string& mods, int maxCharge, bool url ) const
{
	start ( os, url );
	printCGIString ( os, "title", title );
	writeLink ( os, sequence, mods, maxCharge, url );
}
void MSProductLink::write6 ( ostream& os, const string& index, const string& sequence, const string& mods, int maxCharge, bool url ) const
{
	start ( os, url );
	printCGIString ( os, "index", index );
	writeLink ( os, sequence, mods, maxCharge, url );
}
void MSProductLink::write7 ( ostream& os, const string& mOverZ, const string& sequence, const string& mods, int maxCharge, bool url ) const
{
	start ( os, url );
	printCGIString ( os, "m_over_z", mOverZ );
	writeLink ( os, sequence, mods, maxCharge, url );
}
void MSProductLink::write8 ( ostream& os, const string& scanNumber, const string& sequence, const string& mods, int maxCharge, bool url ) const
{
	start ( os, url );
	printCGIString ( os, "scan_number", scanNumber );
	writeLink ( os, sequence, mods, maxCharge, url );
}
void MSProductLink::writeLink ( ostream& os, const string& sequence, const string& mods, int maxCharge, bool url ) const
{
	if ( maxCharge == 0 )printCGIString ( os, "max_charge", "No Limit" );
	else				 printCGI ( os, "max_charge", maxCharge );
	printCGI ( os, "msms_precursor_charge", maxCharge );
	printCGIString ( os, "sequence", sequence );
	printCGI ( os, "s", 1 );
	if ( !mods.empty () )	printCGIString ( os, "mods", mods );
	printCGI ( os, "alternative", true );
	end ( os, url, sequence );
}
void MSProductLink::printHTML ( ostream& os ) const
{
	os << getLinkName () << "=\"";
	os << ProgramLink::getURLStart ( "mssearch" );
	os << "?";
	putCGI ( os );
	os << "\";\n";
}
