/******************************************************************************
*                                                                             *
*  Program    : searchCompare                                                 *
*                                                                             *
*  Filename   : sc_search_res.cpp                                             *
*                                                                             *
*  Created    : March 27th 2003                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif
#define SEARCH_RES_MAIN
#include <cmath>
#include <lg_stdlib.h>
#include <lg_time.h>
#include <lgen_file.h>
#include <lgen_math.h>
#include <lu_getfil.h>
#include <lu_ambiguity.h>
#include <lu_cgi_val.h>
#include <lu_app_gr.h>
#include <lu_check_db.h>
#include <lu_prod_par.h>
#include <lu_get_link.h>
#include <lu_html.h>
#include <lu_html_form.h>
#include <lu_quan_ratio.h>
#include <lu_pi.h>
#include <lu_sim_ent.h>
#include <lu_delim.h>
#include <lu_table.h>
#include <lu_proj_file.h>
#include <lu_usermod.h>
#include <lu_sctag_link.h>
#include <lu_xml.h>
#include <lu_mat_score.h>
#include <lu_file_type.h>
#include <lu_mass_seq.h>
#include <lu_fasta.h>
#include <lu_db_entry.h>
#include <lu_inst.h>
#include <lu_param_list.h>
#include <lu_fas_enz.h>
#include <lu_species.h>
#include <lu_srch_form.h>
#include <sc_anum_res.h>
#include <sc_search_res.h>
#include <sc_mzidentml.h>
#include <sc_quan.h>
#include <sc_xlink.h>
using namespace FileTypes;
using std::vector;
using std::string;
using std::unique;
using std::stable_sort;
using std::cout;
using std::ostream;
using std::istream;
using std::endl;
using std::pair;
using std::remove_if;
using std::runtime_error;
using std::ostringstream;
using std::copy;
using std::back_inserter;
using std::merge;
using std::make_pair;

IntVector ProteinInfo::dbIdx;
MapStringToInt ProteinInfo::dbNameMap;
vector <FastaServer*> ProteinInfo::fs;
vector <AccessionNumberMap*> ProteinInfo::am;
StringVector ProteinInfo::dbase;
vector <AccessionNumberLinkInfo> ProteinInfo::anli;
vector <AccessionNumberLinkInfo> ProteinInfo::upidli;
StringVector ProteinInfo::tempDirs;

IntVector ProteinInfo::initialise ( const ParameterList* params )
{
	IntVector iv;
	StringVector sv = params->getStringVectorValue ( "database" );
	StringVector dbs;
	DBSearchFlags dbFlags ( sv );
	for ( int a = 0 ; a < sv.size () ; a++ ) {
		string d = sv [a];
		if ( d == "User Protein" ) {
			dbs.push_back ( "User Protein" );
			if ( dbFlags.getRandomFlag () ) dbs.push_back ( "User Protein Random" );
			if ( dbFlags.getReverseFlag () )dbs.push_back ( "User Protein Reverse" );
		}
		else {
			PairStringString pss;
			if ( getConcatDBPair ( d, pss ) ) {
				dbs.push_back ( pss.first );
				dbs.push_back ( pss.second );
			}
			else
				dbs.push_back ( d );
		}
	}
	for ( int i = 0 ; i < dbs.size () ; i++ ) {
		string db = dbs [i];
		bool createDB = false;
		if ( db == "User Protein" ) {
			SingleEntryParameters sep ( params );
			PairStringBool psb = sep.createTemporaryDatabase ( false, false );
			db = psb.first;
			tempDirs.push_back ( genDirectoryFromPath ( db ) );
			createDB = psb.second;
		}
		else if ( db != "User Protein Random" && db != "User Protein Reverse" ) {
			if ( !genFileExists ( SeqdbDir::instance ().getDatabasePath ( db ) ) ) {
				db = getBestSubstituteDatabase ( db );
			}
		}
		MapStringToIntConstIterator cur = dbNameMap.find ( db );
		if ( cur == dbNameMap.end () ) {
			dbIdx.push_back ( fs.size () );
			dbNameMap [db] = fs.size ();
			if ( db.empty () || db == "User Protein Random" || db == "User Protein Reverse" ) {
				fs.push_back ( 0 );
				am.push_back ( getAccessionNumberMap ( db ) );
			}
			else {
				fs.push_back ( new FastaServer ( db, createDB ) );
				am.push_back ( getAccessionNumberMap ( fs.back ()->getFileName () ) );
			}
			dbase.push_back ( db );
			if ( reportLinks ) {
				anli.push_back ( AccessionNumberLinkInfo () );
				upidli.push_back ( AccessionNumberLinkInfo ( "upidlinks.txt" ) );
			}
		}
		else
			dbIdx.push_back ( (*cur).second );
		iv.push_back ( dbIdx.back () + 1 );
	}
	return iv;
}
void ProteinInfo::initialiseAccessionNumberLink ( ostream& os )
{
	if ( reportLinks ) {
		startJavascript ( os );
		for ( StringVectorSizeType i = 0 ; i < dbase.size () ; i++ ) {
			anli [i].printHTML ( os, dbase [i] );
			upidli [i].printHTML ( os, dbase [i] );
		}
		endJavascript ( os );
	}
}

bool ProteinInfo::reportUniprotID = true;
bool ProteinInfo::reportGeneName = true;
bool ProteinInfo::reportAccession = true;
bool ProteinInfo::reportVersion = false;
bool ProteinInfo::reportIndex = false;
bool ProteinInfo::reportLength = false;
bool ProteinInfo::reportMW = true;
bool ProteinInfo::reportPI = true;
bool ProteinInfo::reportSpecies = true;
bool ProteinInfo::reportName = true;
bool ProteinInfo::reportLinks = false;
vector <TaxonomyMatch*> ProteinInfo::taxMatch;
bool ProteinInfo::taxCheck = false;
int ProteinInfo::preferredMatchLength = 0;

ProteinInfo::ProteinInfo () :
	proteinMW ( 0.0 ),
	proteinPI ( 0.0 )
{
}
ProteinInfo::ProteinInfo ( const string& aNum ) :
	aNum ( aNum ),
	dnaReadingFrame ( -1 ),
	openReadingFrame ( -1 )
{
	setAnumInfo ( aNum );
	setFields ();
}
void ProteinInfo::setAnumInfo ( const string& str )
{
	string::size_type start = 0;
	string::size_type end = 0;
	string s1 = genNextString ( str, "$", start, end );
	int idx = 1;
	if ( s1.empty () ) {		// No $: this is an accession number
		acc = aNum;
	}
	else {
		string s2 = genNextString ( str, "$", start, end );
		if ( s2.empty () ) {	// Single $: this is database # then accession number
			idx = atoi ( s1.c_str () );
			acc = str.substr ( start );
		}
		else {
			string s3 = genNextString ( str, "$", start, end );
			if ( s3.empty () ) {// Two $: this is dna database
				acc = s1;
				dnaReadingFrame = atoi ( s2.c_str () );
				openReadingFrame = atoi ( str.substr ( start ).c_str () );
			}
			else {				// Three $: this is dna database
				idx = atoi ( s1.c_str () );
				acc = s2;
				dnaReadingFrame = atoi ( s3.c_str () );
				openReadingFrame = atoi ( str.substr ( start ).c_str () );
			}
		}
	}
	dbaseIndex = idx - 1;
}
void ProteinInfo::setFields ()
{
	if ( dbaseIndex < am.size () )
		indexNum = am [dbaseIndex]->getIndexNumber ( acc );
	else
		indexNum = -1;
	if ( indexNum != -1 ) {
		for ( fs [dbaseIndex]->firstLine ( indexNum ) ; fs [dbaseIndex]->isDoneLine () ; fs [dbaseIndex]->nextLine () ) {
			accessionNumbers.push_back ( fs [dbaseIndex]->getLineAccessionNumber () );
			species.push_back ( fs [dbaseIndex]->getLineSpecies () );
			names.push_back ( fs [dbaseIndex]->getLineName () );
			uniprotIDs.push_back ( fs [dbaseIndex]->getLineUniprotID () );

			organismNames.push_back ( fs [dbaseIndex]->getLineOrganismName () );
			geneNames.push_back ( fs [dbaseIndex]->getLineGeneName () );
			proteinExistences.push_back ( fs [dbaseIndex]->getLineProteinExistence () );
			sequenceVersion.push_back ( fs [dbaseIndex]->getLineSequenceVersion () );

			accessionInfo.push_back ( fs [dbaseIndex]->getLineAccessionInfo () );
		}
		char* seq = fs [dbaseIndex]->getProtein ( DatabaseEntry ( indexNum, dnaReadingFrame, openReadingFrame-1 ) );
		ProteinMW pmw ( seq );
		ProteinPI ppi ( seq );
		proteinMW = pmw.getMass ();
		proteinPI = ppi.getProteinPI ();
		length = strlen ( seq );
	}
	else {
		accessionNumbers.push_back ( acc );
		accessionInfo.push_back ( "" );
		species.push_back ( "" );
		names.push_back ( "" );
		uniprotIDs.push_back ( "" );

		organismNames.push_back ( "" );
		geneNames.push_back ( "" );
		proteinExistences.push_back ( "" );
		sequenceVersion.push_back ( "" );

		proteinMW = 0.0;
		proteinPI = 0.0;
		length = 0;
	}
	databaseIndex = dbaseIndex;
}
string ProteinInfo::getProteinSequence () const
{
	string s;
	if ( indexNum != -1 ) s = fs [databaseIndex]->getProtein ( DatabaseEntry ( indexNum, dnaReadingFrame, openReadingFrame-1 ) );
	return s;
}
bool ProteinInfo::isDecoy () const
{
	if ( isPrefix ( accessionNumbers [getPreferredSpeciesIndex ()], "-" ) || isFullyDecoyDatabase ( dbase [databaseIndex] ) )
		return true;
	else
		return false;
}
bool ProteinInfo::isDecoy ( const string& a )
{
	if ( isPrefix ( a, "-" ) ) {
		return true;
	}
	string::size_type start = 0;
	string::size_type end = 0;
	string s1 = genNextString ( a, "$", start, end );
	int idx = 1;
	if ( !s1.empty () ) {
		idx = atoi ( s1.c_str () );
	}
	if ( isFullyDecoyDatabase ( dbase [idx-1] ) ) {
		return true;
	}
	return false;
}
int ProteinInfo::getColspan ()
{
	int colspan = 0;
	if ( reportLength )	colspan++;
	if ( reportMW )		colspan++;
	if ( reportPI )		colspan++;
	if ( reportSpecies )colspan++;
	if ( reportName )	colspan++;
	return colspan;
}
void ProteinInfo::printHTMLHeader ( ostream& os, int rowspan )
{
	if ( reportLength )	tableHeader ( os, "Protein Length", "", "", false, 0, rowspan );
	if ( reportMW )		tableHeader ( os, "Protein MW", "", "", false, 0, rowspan );
	if ( reportPI )		tableHeader ( os, "Protein pI", "", "", false, 0, rowspan );
	if ( reportSpecies )tableHeader ( os, "Species", "", "", false, 0, rowspan );
	if ( reportName )	tableHeader ( os, "Protein Name", "", "", false, 0, rowspan );
}
void ProteinInfo::printHTMLANumHeader ( ostream& os, int rowspan )
{
	if ( reportAccession )	tableHeader ( os, "Acc #", "", "", false, 0, rowspan );
	if ( reportVersion )	tableHeader ( os, "Version", "", "", false, 0, rowspan );
	if ( reportIndex )		tableHeader ( os, "Index", "", "", false, 0, rowspan );
	if ( reportUniprotID )	tableHeader ( os, "UniProt ID", "", "", false, 0, rowspan );
	if ( reportGeneName )	tableHeader ( os, "Gene", "", "", false, 0, rowspan );
}
int ProteinInfo::getPreferredSpeciesIndex () const
{
	if ( species [0].empty () ) {
		return 0;
	}
	else {
		if ( taxCheck ) {
			return taxMatch [0]->getBestSpeciesIndex ( species );
		}
		else {
			return 0;
		}
	}
}
void ProteinInfo::printHTML ( ostream& os, bool empty ) const
{
	if ( empty || accessionNumbers.empty () ) {
		if ( reportLength )	tableEmptyCell ( os );
		if ( reportMW )		tableEmptyCell ( os );
		if ( reportPI )		tableEmptyCell ( os );
		if ( reportSpecies )tableEmptyCell ( os );
		if ( reportName )	tableEmptyCell ( os );
	}
	else {
		if ( reportLength )	tableCell ( os, length );
		if ( reportMW )		tableCell ( os, proteinMW, 1 );
		if ( reportPI )		tableCell ( os, proteinPI, 1 );
		if ( reportSpecies || reportName ) {
			int prefSpInd = getPreferredSpeciesIndex ();
			if ( reportSpecies )tableCell ( os, species [prefSpInd] );
			if ( reportName )	tableCell ( os, names [prefSpInd] );
		}
	}
}
void ProteinInfo::printHTMLANum ( ostream& os, bool empty ) const
{
	if ( empty || accessionNumbers.empty () ) {
		if ( reportAccession )	tableEmptyCell ( os );
		if ( reportVersion )	tableEmptyCell ( os );
		if ( reportIndex )		tableEmptyCell ( os );
		if ( reportUniprotID )	tableEmptyCell ( os );
		if ( reportGeneName )	tableEmptyCell ( os );
	}
	else {
		int prefSpInd = getPreferredSpeciesIndex ();
		if ( reportAccession ) {
			tableCellStart ( os, "", "", true );
			anli [databaseIndex].write2 ( os, accessionNumbers [prefSpInd], reportLinks && indexNum != -1 );
			if ( dnaReadingFrame != -1 ) os << "." << dnaReadingFrame;
			if ( openReadingFrame != -1 ) os << "." << openReadingFrame;
			tableCellEnd ( os );
		}
		if ( reportVersion )	tableCell ( os, sequenceVersion [prefSpInd] );
		if ( reportIndex )		tableCell ( os, indexNum );
		if ( reportUniprotID ) {
			tableCellStart ( os );
			upidli [databaseIndex].write2 ( os, uniprotIDs [prefSpInd], reportLinks && indexNum != -1 );
			tableCellEnd ( os );
		}
		if ( reportGeneName )	tableCell ( os, geneNames [prefSpInd] );
	}
}
void ProteinInfo::printHTMLLines ( ostream& os ) const
{
	if ( !accessionNumbers.empty () ) {
		for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
			os << "<b>Acc. #: </b>";
			anli [databaseIndex].write2 ( os, accessionNumbers [i], reportLinks && indexNum != -1 );
			os << " ";
			if ( !uniprotIDs [i].empty () ) {
				os << "<b>Uniprot ID: </b>";
				upidli [databaseIndex].write2 ( os, uniprotIDs [i], reportLinks && indexNum != -1 );
				os << " ";
			}
			os << " ";
			os << "<b>Species: </b>";
			os << species [i];
			os << " ";

			os << "<b>Name: </b>";
			os << names [i];
			os << "<br />" << endl;

			string str;
			if ( !organismNames [i].empty () )		str += "<b>Organism: </b>"			+ organismNames [i] + " ";
			if ( !geneNames [i].empty () )			str += "<b>Gene: </b>"				+ geneNames [i] + " ";
			if ( !proteinExistences [i].empty () )	str += "<b>Existence: </b>"			+ proteinExistences [i] + " ";
			if ( !sequenceVersion [i].empty () )	str += "<b>Version: </b>"			+ sequenceVersion [i] + " ";
			if ( !str.empty () ) os << str.substr ( 0, str.length () - 1 ) << "<br />" << endl;
		}

		string str2;
		if ( dnaReadingFrame != -1 )	str2 += "<b>DNA Reading Frame: </b>"	+ gen_itoa ( dnaReadingFrame ) + " ";
		if ( openReadingFrame != -1 )	str2 += "<b>Open Reading Frame:</b> "	+ gen_itoa ( openReadingFrame ) + " ";
		if ( !str2.empty () ) os << str2.substr ( 0, str2.length () - 1 ) << "<br />" << endl;

		os << "<b>Protein MW: </b>";
		genPrint ( os, proteinMW, 1 );
		os << " ";
		os << "<b>Protein pI: </b>";
		genPrint ( os, proteinPI, 1 );
		os << " ";
		os << "<b>Protein Length: </b>" << length;
		os << " ";
		os << "<b>Index: </b>" << indexNum;
		os << "<br />" << endl;
	}
}
void ProteinInfo::printDelimitedANumHeader ( ostream& os )
{
	if ( reportAccession )	delimitedHeader ( os, "Acc #" );
	if ( reportVersion )	delimitedHeader ( os, "Version" );
	if ( reportIndex )		delimitedHeader ( os, "Index" );
	if ( reportUniprotID )	delimitedHeader ( os, "UniProt ID" );
	if ( reportGeneName )	delimitedHeader ( os, "Gene" );
}
void ProteinInfo::printDelimitedHeader ( ostream& os )
{
	if ( reportLength )	delimitedHeader ( os, "Protein Length" );
	if ( reportMW )		delimitedHeader ( os, "Protein MW" );
	if ( reportPI )		delimitedHeader ( os, "Protein pI" );
	if ( reportSpecies )delimitedHeader ( os, "Species" );
	if ( reportName )	delimitedHeader ( os, "Protein Name" );
}
void ProteinInfo::printDelimited ( ostream& os ) const
{
	if ( accessionNumbers.empty () ) {
		if ( reportLength )	delimitedEmptyCell ( os );
		if ( reportMW )		delimitedEmptyCell ( os );
		if ( reportPI )		delimitedEmptyCell ( os );
		if ( reportSpecies )delimitedEmptyCell ( os );
		if ( reportName )	delimitedEmptyCell ( os );
	}
	else {
		if ( reportLength )	delimitedCell ( os, length );
		if ( reportMW )		delimitedCell ( os, proteinMW, 1 );
		if ( reportPI )		delimitedCell ( os, proteinPI, 1 );
		if ( reportSpecies || reportName ) {
			int prefSpInd = getPreferredSpeciesIndex ();
			if ( reportSpecies )delimitedCell ( os, species [prefSpInd] );
			if ( reportName )	delimitedCell ( os, names [prefSpInd] );
		}
	}
}
void ProteinInfo::printDelimitedANum ( ostream& os ) const
{
	if ( accessionNumbers.empty () ) {
		if ( reportAccession )	delimitedEmptyCell ( os );
		if ( reportVersion )	delimitedEmptyCell ( os );
		if ( reportIndex )		delimitedEmptyCell ( os );
		if ( reportUniprotID )	delimitedEmptyCell ( os );
		if ( reportGeneName )	delimitedEmptyCell ( os );
	}
	else {
		int prefSpInd = getPreferredSpeciesIndex ();
		string acc = accessionNumbers [prefSpInd];
		bool decoy = isPrefix ( acc, "-" );
		if ( reportAccession ) {
			ostringstream ost;
			if ( decoy || isFullyDecoyDatabase ( dbase [databaseIndex] ) )
				ost << "decoy";
			else
				ost << acc;
			if ( dnaReadingFrame != -1 ) ost << "$" << dnaReadingFrame;
			if ( openReadingFrame != -1 ) ost << "$" << openReadingFrame;

			delimitedCell ( os, ost.str () );
		}
		if ( reportVersion )	delimitedCell ( os, sequenceVersion [prefSpInd] );
		if ( reportIndex )		delimitedCell ( os, indexNum );
		if ( reportUniprotID ) {
			if ( decoy || isFullyDecoyDatabase ( dbase [databaseIndex] ) )
				delimitedCell ( os, "decoy" ); 
			else
				delimitedCell ( os, uniprotIDs [prefSpInd] );
		}
		if ( reportGeneName ) delimitedCell ( os, geneNames [prefSpInd] );
	}
}
void ProteinInfo::setTaxonomyMatch ( const StringVector& preferredSpecies )
{
	if ( !preferredSpecies.empty () ) {
		try {
			for ( StringVectorSizeType i = 0 ; i < preferredSpecies.size () ; i++ ) {
				taxMatch.push_back ( new TaxonomyMatch ( preferredSpecies [i] ) );
			}
			taxCheck = true;
		}
		catch ( runtime_error e ) {
			ErrorHandler::genError ()->error ( e );
		}
	}
}
void ProteinInfo::initProteinSequence () const
{
	if ( indexNum != -1 ) proteinSeq = fs [databaseIndex]->getProtein ( DatabaseEntry ( indexNum, dnaReadingFrame, openReadingFrame-1 ) );
}
string ProteinInfo::getPreviousAA ( int startAA, int numAA ) const
{
	if ( indexNum == -1 || startAA > length ) return "*";
	if ( proteinSeq.empty () ) initProteinSequence ();
	int numDisp = genMin ( numAA, startAA - 1 );
	string ret;
	if ( numDisp < numAA ) ret += string ( "-" );
	ret += proteinSeq.substr ( startAA-numDisp-1, numDisp );
	return ret;
}
string ProteinInfo::getNextAA ( int endAA, int numAA ) const
{
	if ( indexNum == -1 || endAA > length ) return "*";
	if ( proteinSeq.empty () ) initProteinSequence ();
	int numDisp = genMin ( numAA, length - endAA );
	string ret = proteinSeq.substr ( endAA, numDisp );
	if ( numDisp < numAA ) ret += string ( "-" );
	return ret;
}
time_t ProteinInfo::getDatabaseTime ()
{
	if ( fs.empty () )	return time ( 0 );
	else				return fs [0]->getDatabaseTime ();
}
string ProteinInfo::getDatabasePath ()
{
	if ( fs.empty () )	return "";
	else				return fs [0]->getDatabasePath ();
}
string ProteinInfo::getDatabaseName ( int i )
{
	if ( fs.empty () )	return "";
	else				return fs [i]->getFileName ();
}
string ProteinInfo::getDatabasePath ( int i )
{
	if ( i < fs.size () )	return fs [i]->getDatabasePath ();
	else					return "";
}
int ProteinInfo::fsSize ()
{
	return fs.size ();
}
bool ProteinInfo::getDNADatabase ()
{
	if ( fs.empty () )	return false;
	else				return fs [0]->getDNADatabase ();
}
int ProteinInfo::getNumEntries ()
{
	if ( fs.empty () )	return -1;
	else				return fs [0]->getNumEntries ();
}
int ProteinInfo::getNumEntries ( int i )
{
	if ( i < fs.size () )	return fs [i]->getNumEntries ();
	else					return -1;
}
PairIntString ProteinInfo::getANumPair ( const string& aNum )
{
	int pos = aNum.find ( '$' );
	return make_pair ( atoi ( aNum.substr ( 0, pos ).c_str () ) - 1, aNum.substr ( pos+1 ) );
}
StringVector ProteinInfo::getDBSearchList1 ()
{
	StringVector sv = dbase;
	for ( int i = 1 ; i < sv.size () ; i++ ) {
		string db = sv [i];
		int len = 0;
		if ( isSuffix ( db, ".reverse" ) )	len = 8;
		if ( isSuffix ( db, ".random" ) )	len = 7;
		if ( len && sv [i-1] == db.substr ( 0, db.length () - len ) ) sv [i-1] = db + ".concat";
	}
	return sv;
}
BoolDeque ProteinInfo::getDBSearchFlagList1 ()
{
	BoolDeque bd;
	for ( int i = 0 ; i < dbase.size () ; i++ ) {
		bool flag = true;
		string db = dbase [i];
		if ( isSuffix ( db, "UserProtein.fasta" ) )	flag = false;
		if ( db == "User Protein Reverse" )			flag = false;
		if ( db == "User Protein Random" )			flag = false;
		if ( isSuffix ( db, ".reverse" ) || isSuffix ( db, ".random" ) ) flag = false;
		bd.push_back ( flag );
	}
	return bd;
}
bool ProteinInfo::getNonDecoyFlag ( int idx )
{
	static BoolDeque nonDecoyDBFlag;
	if ( nonDecoyDBFlag.empty () ) {	// initialise if necessary
		for ( int i = 0 ; i < dbase.size () ; i++ ) {
			bool flag = true;
			string db = dbase [i];
			if ( db == "User Protein Reverse" )			flag = false;
			if ( db == "User Protein Random" )			flag = false;
			if ( isSuffix ( db, ".reverse" ) || isSuffix ( db, ".random" ) ) flag = false;
			nonDecoyDBFlag.push_back ( flag );
		}
	}
	return nonDecoyDBFlag [idx];
}
void ProteinInfo::getUserProteinIndex ( int& index, const string& db, const string& possDB, const string& diskDB )
{
	if ( isPrefix ( db, possDB ) ) {
		string aUPIdx = gen_strtrim ( db.substr ( possDB.length () ) );
		int upIdx = 1;
		if ( !aUPIdx.empty () && genStringIsInteger ( aUPIdx ) ) upIdx = atoi ( aUPIdx.c_str () );
		for ( StringVectorSizeType i = 0, j = 0 ; i < dbase.size () ; i++ ) {
			if ( isSuffix ( dbase [i], diskDB ) ) {
				j++;
				if ( j == upIdx ) index = i+1;
			}
		}
	}
}
int ProteinInfo::getDBIndex ( const string& db )
{
	int index = -1;
	
	for ( StringVectorSizeType i = 0 ; i < dbase.size () ; i++ ) {	//Look for an exact match
		if ( db == dbase [i] ) return i+1;
	}
	string db2;
	if ( isSuffix ( db, ".random.concat" ) ) db2 = db.substr ( 0, db.length () - 14 );
	if ( isSuffix ( db, ".reverse.concat" ) ) db2 = db.substr ( 0, db.length () - 15 );
	if ( !db2.empty () ) {	// Check if this is a split concat database
		for ( StringVectorSizeType i = 0 ; i < dbase.size () ; i++ ) {
			if ( db2 == dbase [i] ) return i+1;
		}
	}
	getUserProteinIndex ( index, db, "User Protein Random", "User Protein Random" );
	if ( index != -1 ) return index;
	getUserProteinIndex ( index, db, "User Protein Reverse", "User Protein Reverse" );
	if ( index != -1 ) return index;
	getUserProteinIndex ( index, db, "User Protein", "UserProtein.fasta" );
	if ( index != -1 ) return index;
	return index;
}
string ProteinInfo::getDatabaseMZIdentMLRef () const
{
	return "SDB_" + gen_itoa ( dbaseIndex+1 );
}
void ProteinInfo::addDatabaseMZIdentMLInfo ( VectorXMLOutputItemPtr& subItems )
{
	for ( SearchResultsPtrVectorSizeType i = 0 ; i < fs.size () ; i++ ) {
		VectorXMLOutputItemPtr subItems1A;
		VectorXMLOutputItemPtr subItems1A1;
		subItems1A1.push_back ( new MZIdentML_CVParam ( "MS:1001348", "FASTA format", "PSI-MS" ) );
		subItems1A.push_back ( new MZIdentML_FileFormat ( subItems1A1 ) );
		VectorXMLOutputItemPtr subItems1A2;
		subItems1A2.push_back ( new MZIdentML_userParam ( getDatabaseName ( i ) ) );
		subItems1A.push_back ( new MZIdentML_DatabaseName ( subItems1A2 ) );
		subItems1A.push_back ( new MZIdentML_CVParam ( "MS:1001073", "database type amino acid", "PSI-MS" ) );
		subItems.push_back ( new MZIdentML_SearchDatabase ( subItems1A, "SDB_" + gen_itoa ( i+1 ), getDatabasePath ( i ), "", getNumEntries ( i ) ) );
	}
}
void ProteinInfo::resetANumMap ()
{
	for ( int i = 0 ; i < am.size () ; i++ ) {
		am [i]->reset ();
	}
}
void ProteinInfo::deleteTempDirs ()
{
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		delete fs [i];
	}
	for ( StringVectorSizeType j = 0 ; j < tempDirs.size () ; j++ ) {
		genUnlinkDirectory ( tempDirs [j] );
	}
}
MapStringToDouble HitPeptide::msd;
HitPeptide::HitPeptide ( const string& nTerm, const string& peptide, const string& cTerm, const string& neutralLoss ) :
	nTerm ( nTerm ),
	peptide ( peptide ),
	cTerm ( cTerm ),
	neutralLoss ( neutralLoss ),
	mmod ( false )
{
	string s = gen_strstriptags2 ( peptide, '(', ')' );
	if ( s != peptide ) {
		databasePeptide = s;
		for ( StringSizeType i = 0 ; i < peptide.length () ; i++ ) {
			if ( peptide [i] == '(' ) {
				if ( genIsNumberStart ( peptide [i+1] ) ) {
					mmod = true;
					break;
				}
				else {
					int bracket = 0;
					for ( ; i < peptide.length () ; i++ ) {
						char a = peptide [i];
						if ( a == '(' ) bracket++;
						if ( a == ')' ) bracket--;
						if ( bracket == 0 ) break;
					}
				}
			}
		}
	}
	mmod = mmod || ( !nTerm.empty () && genIsNumberStart ( nTerm [0] ) );
	mmod = mmod || ( !cTerm.empty () && genIsNumberStart ( cTerm [0] ) );
	mmod = mmod || ( !neutralLoss.empty () && genIsNumberStart ( neutralLoss [0] ) );
}
string HitPeptide::getBlibPeptide () const
{
	string bPep;
	double nTermMass = 0.0;
	if ( !nTerm.empty () ) {
		if ( genIsNumberStart ( nTerm [0] ) )	nTermMass = atof ( nTerm.c_str () );
		else									nTermMass = getModMass ( nTerm );
	}
	double cTermMass = 0.0;
	if ( !cTerm.empty () ) {
		if ( genIsNumberStart ( cTerm [0] ) )	cTermMass = atof ( cTerm.c_str () );
		else									cTermMass = getModMass ( cTerm );
	}
	StringSizeType len = peptide.length ();
	for ( StringSizeType i = 0 ; i < len ; i++ ) {
		if ( peptide [i] == '(' ) {
			int bracket = 0;
			StringSizeType start = i+1;
			bool ntFlag = ( i == 1 );
			for ( ; i < len ; i++ ) {
				char a = peptide [i];
				if ( a == '(' ) bracket++;
				if ( a == ')' ) bracket--;
				if ( bracket == 0 ) break;
			}
			const string& mod = peptide.substr ( start, i-start );
			double mm = getModMass ( mod );
			if ( ntFlag ) mm += nTermMass;
			if ( i == len - 1 ) mm += cTermMass;
			bPep += "[" + gen_ftoa ( mm, "%+.1f" ) + "]";
		}
		else {
			double mm = 0.0;
			if ( i == 0 && ( i == len - 1 || peptide [i+1] != '(' ) ) mm += nTermMass;
			if ( i == len - 1 )	mm += cTermMass;
			bPep += peptide [i];
			if ( mm ) bPep += "[" + gen_ftoa ( mm, "%+.1f" ) + "]";
		}
	}
	return bPep;
}
void HitPeptide::getBlibMods ( VectorPairIntDouble& vpid ) const
{
	double nTermMass = 0.0;
	if ( !nTerm.empty () ) {
		if ( genIsNumberStart ( nTerm [0] ) )	nTermMass = atof ( nTerm.c_str () );
		else									nTermMass = getModMass ( nTerm );
	}
	double cTermMass = 0.0;
	if ( !cTerm.empty () ) {
		if ( genIsNumberStart ( cTerm [0] ) )	cTermMass = atof ( cTerm.c_str () );
		else									cTermMass = getModMass ( cTerm );
	}
	StringSizeType len = peptide.length ();
	int idx = 1;
	for ( StringSizeType i = 0 ; i < len ; i++ ) {
		if ( peptide [i] == '(' ) {
			int bracket = 0;
			StringSizeType start = i+1;
			bool ntFlag = ( i == 1 );
			for ( ; i < len ; i++ ) {
				char a = peptide [i];
				if ( a == '(' ) bracket++;
				if ( a == ')' ) bracket--;
				if ( bracket == 0 ) break;
			}
			const string& mod = peptide.substr ( start, i-start );
			double mm = getModMass ( mod );
			if ( ntFlag ) mm += nTermMass;
			if ( i == len - 1 ) mm += cTermMass;
			vpid.push_back ( make_pair ( idx-1, mm ) );
		}
		else {
			double mm = 0.0;
			if ( i == 0 && ( i == len - 1 || peptide [i+1] != '(' ) ) mm += nTermMass;
			if ( i == len - 1 )	mm += cTermMass;
			if ( mm ) vpid.push_back ( make_pair ( idx, mm ) );
			idx++;
		}
	}
}
double HitPeptide::getMModValue () const
{
	if ( mmod ) {
		int nMods = 0;
		double val = 0.0;
		for ( StringSizeType i = 0 ; i < peptide.length () ; i++ ) {
			if ( peptide [i] == '(' ) {
				if ( genIsNumberStart ( peptide [i+1] ) ) {
					val = atof ( peptide.substr ( i+1 ).c_str () );
					nMods++;
				}
				int bracket = 0;
				for ( ; i < peptide.length () ; i++ ) {
					char a = peptide [i];
					if ( a == '(' ) bracket++;
					if ( a == ')' ) bracket--;
					if ( bracket == 0 ) break;
				}
			}
		}
		if ( !nTerm.empty () && genIsNumberStart ( nTerm [0] ) ) {
			val = atof ( nTerm.c_str () );
			nMods++;
		}
		if ( !cTerm.empty () && genIsNumberStart ( cTerm [0] ) ) {
			val = atof ( cTerm.c_str () );
			nMods++;
		}
		if ( !neutralLoss.empty () && genIsNumberStart ( neutralLoss [0] ) ) {
			val = atof ( neutralLoss.c_str () );
			nMods++;
		}
		if ( nMods == 1 ) return val;
	}
	return 0.0;
}
int HitPeptide::getMModPosition () const
{
	if ( mmod ) {
		int nMods = 0;
		double val = 0.0;
		size_t len = 0;
		int idx = 0;
		for ( StringSizeType i = 0 ; i < peptide.length () ; i++ ) {
			if ( peptide [i] == '(' ) {
				if ( genIsNumberStart ( peptide [i+1] ) ) {
					idx = len;
					nMods++;
				}
				int bracket = 0;
				for ( ; i < peptide.length () ; i++ ) {
					char a = peptide [i];
					if ( a == '(' ) bracket++;
					if ( a == ')' ) bracket--;
					if ( bracket == 0 ) break;
				}
			}
			else len++;
		}
		if ( !nTerm.empty () && genIsNumberStart ( nTerm [0] ) ) {
			idx = 1;
			nMods++;
		}
		if ( !cTerm.empty () && genIsNumberStart ( cTerm [0] ) ) {
			idx = len;
			nMods++;
		}
		if ( !neutralLoss.empty () && genIsNumberStart ( neutralLoss [0] ) ) {
			idx = 0;
			nMods++;
		}
		if ( nMods == 1 ) return idx;
	}
	return 0;
}
string HitPeptide::getModString () const
{
	string str;
	if ( !nTerm.empty () ) str += nTerm + "@N-term; ";
	for ( StringSizeType i = 0, idx = 0 ; i < peptide.length () ; i++ ) {
		if ( peptide [i] == '(' ) {
			int bracket = 0;
			StringSizeType start = i+1;
			for ( ; i < peptide.length () ; i++ ) {
				char a = peptide [i];
				if ( a == '(' ) bracket++;
				if ( a == ')' ) bracket--;
				if ( bracket == 0 ) break;
			}
			str += peptide.substr ( start, i-start ) + "(" + peptide [start - 2] + ")@" + gen_itoa ( idx ) + "; ";
		}
		else idx++;
	}
	if ( !cTerm.empty () ) str += cTerm + "@C-term; ";
	if ( !neutralLoss.empty () ) str += neutralLoss + "@Neutral loss; ";
	if ( !str.empty () ) str = str.substr ( 0, str.length () - 2 ); // Delete trailing "; "
	return str;
}
void HitPeptide::getModMassesAndIndicies ( VectorPairIntDouble& vpid ) const
{
	for ( StringSizeType i = 0, idx = 0 ; i < peptide.length () ; i++ ) {
		if ( peptide [i] == '(' ) {
			int bracket = 0;
			StringSizeType start = i+1;
			for ( ; i < peptide.length () ; i++ ) {
				char a = peptide [i];
				if ( a == '(' ) bracket++;
				if ( a == ')' ) bracket--;
				if ( bracket == 0 ) break;
			}
			vpid.push_back ( make_pair ( idx, getModMass ( peptide.substr ( start, i-start ) ) ) );
		}
		else idx++;
	}
}
void HitPeptide::getModMassesIndiciesAndString ( VectorPairIntPairStringDouble& vpid ) const
{
	for ( StringSizeType i = 0, idx = 0 ; i < peptide.length () ; i++ ) {
		if ( peptide [i] == '(' ) {
			int bracket = 0;
			StringSizeType start = i+1;
			for ( ; i < peptide.length () ; i++ ) {
				char a = peptide [i];
				if ( a == '(' ) bracket++;
				if ( a == ')' ) bracket--;
				if ( bracket == 0 ) break;
			}
			const string& mod = peptide.substr ( start, i-start );
			vpid.push_back ( make_pair ( idx, make_pair ( mod, getModMass ( mod ) ) ) );
		}
		else idx++;
	}
}
double HitPeptide::getModMass ( const string& s ) const
{
	if ( !s.empty () ) {
		if ( genIsNumberStart ( s [0] ) ) return atof ( s.c_str () );
		else {
			MapStringToDoubleConstIterator cur = msd.find ( s );

			if ( cur != msd.end () ) return (*cur).second;
			else {
				Usermod u ( s, false );
				ElementalFormula ef = u.getElementalFormula ();
				double m = formula_to_monoisotopic_mass ( ef );
				msd [s] = m;
				return m;
			}
		}
	}
	else return 0.0;
}
string HitPeptide::getCommandLineNVPair ( int i ) const
{
	string num = ( i == 1 ) ? "" : gen_itoa ( i );
	string s;
	s += ::getCommandLineNVPair ( "s" + num, "1" );
	s += " ";
	s += ::getCommandLineNVPair ( "sequence" + num, peptide );
	s += " ";
	s += ::getCommandLineNVPair ( "nterm" + num, nTerm );
	s += " ";
	s += ::getCommandLineNVPair ( "cterm" + num, cTerm );
	s += " ";
	s += ::getCommandLineNVPair ( "nloss" + num, neutralLoss );
	s += " ";
	return s;
}
AACalculator* PeptidePosition::aaCalc = 0;
double PeptidePosition::defaultResolution;
int PeptidePosition::numIsotopePeaks = 1;
string PeptidePosition::compMaskType;
unsigned int PeptidePosition::compMask = 0;
MapStringToUInt PeptidePosition::compMaskMap;
StringVector PeptidePosition::instrument;
BoolDeque PeptidePosition::chargeReducedFragmentation;
StringVectorVector PeptidePosition::fractionNames;
BoolDeque PeptidePosition::spottingPlates;
BoolDeque PeptidePosition::spectrumNumber;
StringVectorVector PeptidePosition::rawTypes;
VectorConstParameterListPtr PeptidePosition::params;
StringVector PeptidePosition::searchKey;
vector <Tolerance*> PeptidePosition::parentTolerances;
StringVector PeptidePosition::sysErrorStr;
DoubleVectorVector PeptidePosition::offsets;
vector <Tolerance*> PeptidePosition::fragmentTolerances;
bool PeptidePosition::multipleErrorUnits = false;
bool PeptidePosition::multipleFractionNames = false;
bool PeptidePosition::spottingPlatesFlag = false;
bool PeptidePosition::spectrumNumberFlag = false;
bool PeptidePosition::quanMultiNormalFlag = false;
bool PeptidePosition::enzymeInit = false;
bool PeptidePosition::reportCheckboxes = false;
bool PeptidePosition::reportSearchNumber = true;
bool PeptidePosition::reportMPlusH = false;
bool PeptidePosition::reportMOverZ = false;
bool PeptidePosition::reportCharge = false;
bool PeptidePosition::reportIntensity = false;
bool PeptidePosition::reportMPlusHCalc = false;
bool PeptidePosition::reportMOverZCalc = false;
bool PeptidePosition::reportError = false;
int PeptidePosition::reportPreviousAA = 0;
bool PeptidePosition::reportDBPeptide = false;
int PeptidePosition::reportNextAA = 0;
string PeptidePosition::peptideModType = "Off";
bool PeptidePosition::reportProteinMods = false;
bool PeptidePosition::reportTime = false;
bool PeptidePosition::reportMSMSInfo = false;
bool PeptidePosition::reportStartAA = false;
bool PeptidePosition::reportEndAA = false;
bool PeptidePosition::reportElemComp = false;
bool PeptidePosition::reportMissedCleavages = false;
bool PeptidePosition::reportLength = false;
bool PeptidePosition::reportComposition = false;
bool PeptidePosition::reportMModValue = false;
bool PeptidePosition::reportLinks = false;
bool PeptidePosition::runMSProductFlag = false;
double PeptidePosition::rtIntervalStart = 0.0;
double PeptidePosition::rtIntervalEnd = 0.0;
const double PeptidePosition::INVALID_ERROR = std::numeric_limits<double>::infinity();
PeptidePosition::PeptidePosition ( const string& accessionNumber, const HitPeptide* hitPeptide, double error, const SpecID* spID, const MSMSSpectrumInfo* mmsi, int startAA, int searchIndex ) :
	accessionNumber ( accessionNumber ),
	hitPeptide ( hitPeptide ),
	error ( error ),
	startAA ( startAA ),
	specID ( spID ),
	mmsi ( mmsi ),
	firstOccurence ( true ), 
	quanRatio ( 0 ),
	searchIndex ( searchIndex )
{
	setSpectrumNumber ( spID, searchIndex );
}
PeptidePosition::PeptidePosition ( const SpecID* spID, const MSMSSpectrumInfo* mmsi, int searchIndex ) :
	specID ( spID ),
	mmsi ( mmsi ),
	hitPeptide ( 0 ),
	firstOccurence ( false ), 
	quanRatio ( 0 ),
	searchIndex ( searchIndex )
{
	setSpectrumNumber ( spID, searchIndex );
	error = INVALID_ERROR;
}
void PeptidePosition::setSpectrumNumber ( const SpecID* spID, int searchIndex )
{
	if ( spID ) {
		bool flag = specID->getSpecNum () != 1;
		if ( flag ) {
			int siz = spectrumNumber.size ();
			if ( searchIndex >= siz ) {
				spectrumNumber.resize ( searchIndex+1 );
				for ( int i = siz ; i <= searchIndex ; i++ ) {
					spectrumNumber [i] = false;
				}
			}
			spectrumNumber [searchIndex] = true;
			spectrumNumberFlag = true;
		}
	}
	if ( spectrumNumber.empty () ) spectrumNumber.push_back ( false );
}
unsigned int PeptidePosition::peptideToMask () const
{
	unsigned int mask = 0;
	const string& pep = getPeptide ();
	for ( StringSizeType i = 0 ; i < pep.length () ; i++ ) {
		string s;
		if ( i < pep.length () - 1 && pep [i+1] == '(' ) {
			string cur;
			cur += pep [i++];
			int bracket = 0;
			for ( ; i < pep.length () ; i++ ) {
				char a = pep [i];
				cur += a;
				if ( a == '(' ) bracket++;
				if ( a == ')' ) bracket--;
				if ( bracket == 0 ) break;
			}
			s = cur;
		}
		else {
			s += pep [i];
		}
		MapStringToUIntConstIterator cur = compMaskMap.find ( s );
		if ( cur != compMaskMap.end () ) mask |= (*cur).second;
	}
	const string& nT = getNTerm ();
	if ( !nT.empty () ) {
		MapStringToUIntConstIterator cur = compMaskMap.find ( "n" + nT );
		if ( cur != compMaskMap.end () ) mask |= (*cur).second;
	}
	const string& cT = getCTerm ();
	if ( !cT.empty () ) {
		MapStringToUIntConstIterator cur = compMaskMap.find ( "c" + cT );
		if ( cur != compMaskMap.end () ) mask |= (*cur).second;
	}
	const string& nL = getNeutralLoss ();
	if ( !nL.empty () ) {
		MapStringToUIntConstIterator cur = compMaskMap.find ( "." + nL );
		if ( cur != compMaskMap.end () ) mask |= (*cur).second;
	}
	return mask;
}
bool PeptidePosition::checkComposition () const
{
	if ( compMask ) {
		if ( compMaskType == "AND" ) {
			if ( ( compMask & peptideToMask () ) == compMask ) {
				return true;
			}
		}
		else if ( compMaskType == "OR" ) {
			if ( ( compMask & peptideToMask () ) != 0 ) {
				return true;
			}
		}
		return false;
	}
	else
		return false;
}
bool PeptidePosition::isDecoyHit () const
{
	if ( accessionNumber.empty () ) return false;
	else {
		PairIntString pis = ProteinInfo::getANumPair ( accessionNumber );
		int idx = pis.first;
		string an = pis.second;
		if ( an [0] != '-' && ProteinInfo::getNonDecoyFlag ( idx ) )
			return false;
		else
			return true;
	}
}
void PeptidePosition::setQuanResults () const
{
	QuanPeptide qp ( getPeptide (), getNTerm (), getCTerm (), getNeutralLoss () );
	quanRatio = PeptidePositionQuan::getQuanRatio ( getMOverZ (), getCharge (), qp, getSpecIDasID (), searchIndex );
}
int PeptidePosition::getColspan ( int searchNumber, bool tabDelim, bool showTimes )
{
	int colspan = 0;
	if ( !tabDelim && reportCheckboxes )	colspan++;
	if ( reportSearchNumber ) colspan++;
	if ( showTimes ) colspan += getColspanPeak1 ();
	if ( reportMPlusHCalc )	colspan++;
	if ( reportMOverZCalc )	colspan++;
	if ( reportError ) {
		if ( multipleErrorUnits && searchNumber == -1 ) colspan += 2;
		else colspan++;
	}
	if ( tabDelim && reportPreviousAA )	colspan++;
	if ( reportDBPeptide )	colspan++;
	if ( tabDelim && reportNextAA )		colspan++;
	if ( peptideModType != "Off" ) {
		if ( peptideModType == "All Mods (2 columns)" )	colspan += 2;
		else											colspan++;
	}
	if ( reportProteinMods )colspan++;
	if ( reportLength )		colspan++;
	if ( reportComposition )colspan++;
	if ( reportMModValue )	colspan++;
	if ( reportMissedCleavages )colspan++;
	if ( showTimes ) colspan += getColspanPeak2 ( searchNumber, tabDelim );
	if ( reportMSMSInfo )	colspan++;
	colspan += PeptidePositionQuan::getColspan ( searchNumber );
	if ( reportStartAA )	colspan++;
	if ( reportEndAA )		colspan++;
	if ( reportElemComp )	colspan++;
	return colspan;
}
int PeptidePosition::getColspanPeak1 ()
{
	int colspan = 0;
	if ( reportMPlusH )		colspan++;
	if ( reportMOverZ )		colspan++;
	if ( reportCharge )		colspan++;
	if ( reportIntensity )	colspan++;
	return colspan;
}
int PeptidePosition::getColspanPeak2 ( int searchNumber, bool tabDelim )
{
	int colspan = 0;
	if ( reportTime ) {
		if ( tabDelim ) colspan++;				// Fraction always reported
		else {
			if ( searchNumber == -1 ) {
				if ( multipleFractionNames ) colspan++;
			}
			else {
				if ( fractionNames [searchNumber].size () > 1 ) colspan++;
			}
		}
		if ( searchNumber == -1 )
			colspan += spottingPlatesFlag ? 3 : ( tabDelim && spectrumNumberFlag ) ? 2 : 1;
		else
			colspan += spottingPlates [searchNumber] ? 3 : ( tabDelim && spectrumNumber [searchNumber] ) ? 2 : 1;
	}
	return colspan;
}
void PeptidePosition::printHeaderHTML ( ostream& os, int searchNumber, const string& styleID, bool showTimes )
{
	if ( reportCheckboxes )	tableHeader ( os, "", styleID );
	if ( reportSearchNumber ) tableHeader ( os, "Search #", styleID );
	if ( showTimes ) printHeaderHTMLPeak1 ( os, styleID );
	if ( reportMPlusHCalc )	tableHeader ( os, "M+H Calc", styleID );
	if ( reportMOverZCalc )	tableHeader ( os, "m/z Calc", styleID );
	if ( reportError ) {
		if ( searchNumber == -1 ) {
			if ( multipleErrorUnits ) {
				tableHeader ( os, "Error", styleID );
				tableHeader ( os, "Units", styleID );
			}
			else
				tableHeader ( os, parentTolerances [0]->getUnitsString (), styleID );
		}
		else
			tableHeader ( os, parentTolerances [searchNumber]->getUnitsString (), styleID );
	}
	if ( reportDBPeptide )	tableHeader ( os, "DB Peptide", styleID );
	if ( peptideModType == "Mods In Peptide" )			tableHeader ( os, "Peptide", styleID );
	else if ( peptideModType == "Variable Mods Only" )	tableHeader ( os, "Variable Mods", styleID );
	else if ( peptideModType == "Constant Mods Only" )	tableHeader ( os, "Constant Mods", styleID );
	else if ( peptideModType == "All Mods (1 column)" )	tableHeader ( os, "Mods", styleID );
	else if ( peptideModType == "All Mods (2 columns)" ) {
		tableHeader ( os, "Constant Mods", styleID );
		tableHeader ( os, "Variable Mods", styleID );
	}
	if ( reportProteinMods )tableHeader ( os, "Protein Mods", styleID );
	if ( reportLength )		tableHeader ( os, "Length", styleID );
	if ( reportComposition )tableHeader ( os, "Composition", styleID );
	if ( reportMModValue )	tableHeader ( os, "M Mod", styleID );
	if ( reportMissedCleavages )tableHeader ( os, "M Cl", styleID );
	if ( showTimes ) printHeaderHTMLPeak2 ( os, searchNumber, styleID );
	if ( reportMSMSInfo )	SpecID::printTableMSMSInfoHeader ( os, styleID );
	PeptidePositionQuan::printHeaderHTML ( os, searchNumber, styleID );
	if ( reportStartAA )	tableHeader ( os, "Start", styleID );
	if ( reportEndAA )		tableHeader ( os, "End", styleID );
	if ( reportElemComp )	tableHeader ( os, "Elemental Composition", styleID );
}
void PeptidePosition::printHeaderHTMLPeak1 ( ostream& os, const string& styleID, int rowspan )
{
	if ( reportMPlusH )		tableHeader ( os, "M+H", styleID, "", false, 0, rowspan );
	if ( reportMOverZ )		tableHeader ( os, "m/z", styleID, "", false, 0, rowspan );
	if ( reportCharge )		tableHeader ( os, "z", styleID, "", false, 0, rowspan );
	if ( reportIntensity )	tableHeader ( os, "Intensity", styleID, "", false, 0, rowspan );
}
void PeptidePosition::printHeaderHTMLPeak2 ( ostream& os, int searchNumber, const string& styleID, int rowspan )
{
	if ( reportTime ) {
		if ( searchNumber == -1 ) {
			if ( multipleFractionNames ) tableHeader ( os, "Fraction", styleID, "", false, 0, rowspan );
			SpecID::printTableHeader ( os, spottingPlatesFlag, spottingPlatesFlag, styleID, rowspan );
		}
		else {
			if ( fractionNames [searchNumber].size () > 1 ) tableHeader ( os, "Fraction", styleID, "", false, 0, rowspan );
			SpecID::printTableHeader ( os, spottingPlates [searchNumber], spottingPlates [searchNumber], styleID, rowspan );
		}
	}
}
string PeptidePosition::getSequence () const
{
	return getNTerm () + '-' + getPeptide () + '-' + getCTerm () + '+' + getNeutralLoss ();
}
string PeptidePosition::getSequencePlusMods () const
{
	return getDBPeptide () + ';' + SCModInfo::getConstModsString ( getDBPeptide (), searchIndex ) + ';' + SCModInfo::getModsString ( searchIndex, specID, getDBPeptide () );
}
string PeptidePosition::getPrintedSequence () const
{
	string str;
	if ( !getNTerm ().empty () ) str += getNTerm () + '-';
	str += getPeptide ();
	if ( !getCTerm ().empty () ) str += '-' + getCTerm ();
	if ( !getNeutralLoss ().empty () ) str += '+' + getNeutralLoss ();
	return str;
}
int PeptidePosition::getNumTolTerm ( const ProteinInfo& proteinInfo ) const
{
	if ( enzymeInit ) {
		string prev = getPrevAA ( proteinInfo );
		if ( prev == "" ) return -1;
		prev = prev.substr ( prev.length () - 1 );
		string p = getDBPeptide ();
		int num = 0;
		if ( prev == "-" ) num++;
		else {
			const IntVector& ci = enzyme_fragmenter ( prev + p );
			if ( ci [0] == 0 ) num++;
		}
		string next = getNextAA ( proteinInfo );
		if ( next == "" ) return -1;
		next = next.substr ( 0, 1 );
		if ( next == "-" ) num++;
		else {
			p += next;
			const IntVector& ci = enzyme_fragmenter ( p );
			if ( ci.size () >= 2 && ci [ci.size ()-2] == p.length () - 2 ) num++;
		}
		return num;
	}
	else return 2;
}
string PeptidePosition::getMissedCleavages () const
{
	if ( enzymeInit ) return gen_itoa ( enzyme_fragmenter ( getDBPeptide () ).size () - 1 );
	else return "-";
}
string PeptidePosition::getElemComp () const
{
	ElementalFormula ef;
	try {
		if ( !aaCalc->calculateStrippedElementalCompositionWithTerminii ( getPeptide (), getNTerm (), getCTerm (), getNeutralLoss (), ef ) ) {
			return "";
		}
	}
	catch ( AACalculatorNoElementalComposition ) {
		return "";
	}
	return ef.getFormula ();
}
void PeptidePosition::printParamsHTML ( ostream& os, int searchNumber )
{
	ExpandableJavascriptBlock ejb ( "Original Search Parameters " + gen_itoa ( searchNumber + 1 ) );

	ejb.printHeader ( os );
	os << "<p>" << endl;
	params [searchNumber]->HTMLParameters ( os );
	os << "</p>" << endl;
	os << "<p>" << endl;
	StringVector& sv = fractionNames [searchNumber];
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		os << "Fraction " << i+1 << ": " << "<b>" << sv [i] << "</b>" << "<br />" << endl;
	}
	os << "</p>" << endl;
	ejb.printFooter ( os );
	os << "<br />" << endl;
}
void PeptidePosition::printHTMLEmpty ( ostream& os, int searchNumber, const string& styleID, bool showTimes )
{
	if ( reportCheckboxes )	tableEmptyCell ( os, styleID );
	if ( reportSearchNumber )tableEmptyCell ( os, styleID );
	if ( showTimes ) printHTMLEmptyPeak1 ( os, styleID );
	if ( reportMPlusHCalc )	tableEmptyCell ( os, styleID );
	if ( reportMOverZCalc )	tableEmptyCell ( os, styleID );
	if ( reportError )		tableEmptyCell ( os, styleID );
	if ( reportDBPeptide )	tableEmptyCell ( os, styleID );
	if ( peptideModType != "Off" ) {
		if ( peptideModType == "All Mods (2 columns)" )	tableEmptyNCells ( os, 2, styleID );
		else											tableEmptyCell ( os, styleID );
	}
	if ( reportProteinMods )tableEmptyCell ( os, styleID );
	if ( reportLength )		tableEmptyCell ( os, styleID );
	if ( reportComposition )tableEmptyCell ( os, styleID );
	if ( reportMModValue )	tableEmptyCell ( os, styleID );
	if ( reportMissedCleavages )tableEmptyCell ( os, styleID );
	if ( showTimes ) printHTMLEmptyPeak2 ( os, searchNumber, styleID );
	if ( reportMSMSInfo )	tableEmptyCell ( os, styleID );
	PeptidePositionQuan::printQuanBlankHTML ( os, searchNumber, styleID );
	if ( reportStartAA )	tableEmptyCell ( os, styleID );
	if ( reportEndAA )		tableEmptyCell ( os, styleID );
	if ( reportElemComp )	tableEmptyCell ( os, styleID );
}
void PeptidePosition::printHTMLEmptyPeak1 ( ostream& os, const string& styleID )
{
	if ( reportMPlusH )		tableEmptyCell ( os, styleID );
	if ( reportMOverZ )		tableEmptyCell ( os, styleID );
	if ( reportCharge )		tableEmptyCell ( os, styleID );
	if ( reportIntensity )	tableEmptyCell ( os, styleID );
}
void PeptidePosition::printHTMLEmptyPeak2 ( ostream& os, int searchNumber, const string& styleID )
{
	if ( reportTime ) {
		if ( fractionNames [searchNumber].size () > 1 ) tableEmptyCell ( os, styleID );
		SpecID::printTableEmptyCell ( os, spottingPlates [searchNumber], spottingPlates [searchNumber], styleID );
	}
}
void PeptidePosition::printDelimitedEmpty ( ostream& os, int searchNumber, bool showTimes )
{
	if ( reportSearchNumber ) delimitedEmptyCell ( os );
	if ( showTimes ) printDelimitedEmptyPeak1 ( os );
	if ( reportMPlusHCalc )	delimitedEmptyCell ( os );
	if ( reportMOverZCalc )	delimitedEmptyCell ( os );
	if ( reportError )		delimitedEmptyCell ( os );
	if ( reportPreviousAA )	delimitedEmptyCell ( os );
	if ( reportDBPeptide )	delimitedEmptyCell ( os );
	if ( reportNextAA )		delimitedEmptyCell ( os );
	if ( peptideModType != "Off" ) {
		if ( peptideModType == "All Mods (2 columns)" )	delimitedEmptyNCells ( os, 2 );
		else											delimitedEmptyCell ( os );
	}
	if ( reportProteinMods )delimitedEmptyCell ( os );
	if ( reportLength )		delimitedEmptyCell ( os );
	if ( reportComposition )delimitedEmptyCell ( os );
	if ( reportMModValue )	delimitedEmptyCell ( os );
	if ( reportMissedCleavages )delimitedEmptyCell ( os );
	if ( showTimes ) printDelimitedEmptyPeak2 ( os, searchNumber );
	if ( reportMSMSInfo )	delimitedEmptyCell ( os );
	PeptidePositionQuan::printQuanBlankDelimited ( os, searchNumber );
	if ( reportStartAA )	delimitedEmptyCell ( os );
	if ( reportEndAA )		delimitedEmptyCell ( os );
	if ( reportElemComp )	delimitedEmptyCell ( os );
}
void PeptidePosition::printDelimitedEmptyPeak1 ( ostream& os )
{
	if ( reportMPlusH )		delimitedEmptyCell ( os );
	if ( reportMOverZ )		delimitedEmptyCell ( os );
	if ( reportCharge )		delimitedEmptyCell ( os );
	if ( reportIntensity )	delimitedEmptyCell ( os );
}
void PeptidePosition::printDelimitedEmptyPeak2 ( ostream& os, int searchNumber )
{
	if ( reportTime ) {
		delimitedEmptyCell ( os );
		SpecID::printDelimitedEmptyCell ( os, spottingPlates [searchNumber], spottingPlates [searchNumber] || spectrumNumber [searchNumber] );
	}
}
void PeptidePosition::printHTML ( ostream& os, const ProteinInfo& proteinInfo, const string& styleID, const SCMSTagLink& smtl, const string& id, bool joint, bool showTimes ) const
{
	MSProductLink* productLink = 0;

	if ( reportCheckboxes ) {
		tableCellStart ( os, styleID );
			string val = gen_itoa ( searchIndex );
			val += '\t';
			val += id;
			val += '\t';
			val += getSpecID ();
			val += '\t';
			//val += accessionNumber;
			//val += '\t';
			val += getPeptide ();
			FormItemCheckbox fic ( "", "", "cb", false, val );
			fic.printHTML ( os );
		tableCellEnd ( os );
	}
	if ( reportSearchNumber ) tableCell ( os, searchIndex+1, false, styleID );
	if ( reportLinks ) {
		productLink = new MSProductLink ( instrument [searchIndex], searchKey [searchIndex], getSpecIDasID (), parentTolerances [searchIndex], fragmentTolerances [searchIndex] );
		startJavascript ( os );
		productLink->printHTML ( os );
		endJavascript ( os );
	}
	if ( showTimes ) printHTMLPeak1 ( os, styleID );
	if ( getPeptide ().empty () ) {
		if ( reportMPlusHCalc )	tableEmptyCell ( os, styleID );
		if ( reportMOverZCalc )	tableEmptyCell ( os, styleID );
		if ( reportError ) {
			tableEmptyCell ( os, styleID );
			if ( joint && multipleErrorUnits ) {
				tableEmptyCell ( os, styleID );
			}
		}
		if ( reportDBPeptide ) {
			tableCellStart ( os, styleID );
				string pks = "Pks: " + gen_itoa ( getNumPeaks () );
				if ( reportLinks )
					productLink->write3 ( os, pks, getMaxFragmentCharge ( searchIndex ) );
				else
					os << pks;
			tableCellEnd ( os );
		}
		if ( peptideModType != "Off" ) {
			if ( peptideModType == "All Mods (2 columns)" )	tableEmptyNCells ( os, 2, styleID );
			else											tableEmptyCell ( os, styleID );
		}
		if ( reportProteinMods )tableEmptyCell ( os, styleID );
		if ( reportLength )		tableCell ( os, "---", false, false, styleID );
		if ( reportComposition )tableCell ( os, "---", false, false, styleID );
		if ( reportMModValue )	tableCell ( os, "---", false, false, styleID );
		if ( reportMissedCleavages )tableCell ( os, "---", false, false, styleID );
	}
	else {
		if ( reportMPlusHCalc )	tableCell ( os, getMPlusHCalc ( searchIndex ), 4, false, styleID );
		if ( reportMOverZCalc )	tableCell ( os, getMOverZCalc ( searchIndex ), 4, false, styleID );
		if ( reportError ) {
			tableCellSigFig ( os, error, 2, false, styleID );
			if ( joint && multipleErrorUnits ) {
				tableCell ( os, parentTolerances [searchIndex]->getUnitsString (), false, false, styleID );
			}
		}
		if ( reportDBPeptide ) {
			tableDataStart ( os, styleID, "", true );
				if ( peptideModType == "Mods In Peptide" ) os << getDBPeptide ();
				else {
					if ( firstOccurence )	os << "<b>";
					if ( reportPreviousAA ) os << "(" << proteinInfo.getPreviousAA ( startAA, reportPreviousAA ) << ")";
					if ( reportLinks ) {
						productLink->write4 ( os, getDBPeptide (), SCModInfo::getAllModsString ( searchIndex, specID, getDBPeptide () ), getMaxFragmentCharge ( searchIndex ) );
					}
					else os << getDBPeptide ();
					if ( reportNextAA ) os << "(" << proteinInfo.getNextAA ( getEndAA (), reportNextAA ) << ")";
					if ( firstOccurence ) os << "</b>";
				}
			tableCellEnd ( os );
		}
		if ( peptideModType != "Off" ) {
			string cMod = SCModInfo::getConstModsString ( getDBPeptide (), searchIndex );
			string vMod = SCModInfo::getModsString ( searchIndex, specID, getDBPeptide () );
			if ( peptideModType == "Variable Mods Only" )		tableCell ( os, vMod, true, false, styleID );
			else if ( peptideModType == "Constant Mods Only" )	tableCell ( os, cMod, true, false, styleID );
			else if ( peptideModType == "Mods In Peptide" )	{
				tableDataStart ( os, styleID, "", true );
					if ( firstOccurence )	os << "<b>";
					if ( reportPreviousAA ) os << "(" << proteinInfo.getPreviousAA ( startAA, reportPreviousAA ) << ")";
					if ( reportLinks ) {
						productLink->write4 ( os, getDBPeptide (), getPeptide (), getNTerm (), getCTerm (), getNeutralLoss (), SCModInfo::getAllModsString ( searchIndex, specID, getDBPeptide () ), getMaxFragmentCharge ( searchIndex ) );
					}
					else os << getDBPeptide ();
					if ( reportNextAA ) os << "(" << proteinInfo.getNextAA ( getEndAA (), reportNextAA ) << ")";
					if ( firstOccurence ) os << "</b>";
				tableCellEnd ( os );
			}
			else if ( peptideModType == "All Mods (1 column)" )	{
				tableCell ( os, cMod.empty () ? vMod : ( vMod.empty () ? cMod : cMod + ';' + vMod ), true, false, styleID );
			}
			else if ( peptideModType == "All Mods (2 columns)" ) {
				tableCell ( os, cMod, true, false, styleID );
				tableCell ( os, vMod, true, false, styleID );
			}
		}
		if ( reportProteinMods )tableCell ( os, SCModInfo::getModsString ( searchIndex, specID, getDBPeptide (), startAA ), false, false, styleID );
		if ( reportLength )		tableCell ( os, getPeptideLength (), false, styleID );
		if ( reportComposition )tableCell ( os, checkComposition (), false, styleID );
		if ( reportMModValue )	{
			double m = getMModValue ();
			if ( m )tableCell ( os, m, 0, false, styleID );
			else	tableEmptyCell ( os, styleID );
		}
		if ( reportMissedCleavages )tableCell ( os, getMissedCleavages (), false, false, styleID );
	}
	if ( showTimes ) printHTMLPeak2 ( os, styleID, smtl, joint );
	if ( reportMSMSInfo ) specID->printTableMSMSInfoCell ( os, styleID );
	if ( getPeptide ().empty () || quanRatio == 0 )
		PeptidePositionQuan::printQuanBlankHTML ( os, searchIndex, styleID );
	else {
		PeptidePositionQuan::printHTML ( os, quanRatio, styleID, joint );
	}
	if ( getPeptide ().empty () ) {
		if ( reportStartAA )tableCell ( os, "---", false, false, styleID );
		if ( reportEndAA )	tableCell ( os, "---", false, false, styleID );
	}
	else {
		if ( reportStartAA )tableCell ( os, startAA, false, styleID );
		if ( reportEndAA )	tableCell ( os, getEndAA (), false, styleID );
	}
	if ( reportElemComp )	tableCell ( os, getElemComp (), true, false, styleID );
	delete productLink;
}
void PeptidePosition::printHTMLPeak1 ( ostream& os, const string& styleID ) const
{
	MSParentLink* parentLink = 0;
	if ( reportLinks ) {
		if ( !rawTypes [searchIndex][getFraction ()-1].empty () ) parentLink = new MSParentLink;
		startJavascript ( os );
		if ( !rawTypes [searchIndex][getFraction ()-1].empty () ) parentLink->printHTML ( os );
		endJavascript ( os );
	}
	if ( reportMPlusH )		tableCell ( os, getMPlusH (), 4, false, styleID );
	if ( reportMOverZ ) {
		if ( reportLinks && !rawTypes [searchIndex][getFraction ()-1].empty () ) {
			tableDataStart ( os, styleID, "", true );
				string units = parentTolerances [searchIndex]->getUnitsString ();
				parentLink->write ( os, getSpecIDasID (), getMOverZ (), getCharge (), rtIntervalStart, rtIntervalEnd, PeakFitData::getSNRThreshold (), getPeptide (), getNTerm (), getCTerm (), getNeutralLoss (), searchKey [searchIndex], sysErrorStr [searchIndex], units );
			tableCellEnd ( os );
		}
		else
			tableCell ( os, getMOverZ (), 4, false, styleID );
	}
	if ( reportCharge )		tableCell ( os, getCharge (), false, styleID );
	if ( reportIntensity )	tableCellSigFig ( os, getIntensity (), 3, false, styleID );
}
void PeptidePosition::printHTMLPeak2 ( ostream& os, const string& styleID, const SCMSTagLink& smtl, bool joint ) const
{
	if ( reportTime ) {
		if ( ( joint && multipleFractionNames ) || fractionNames [searchIndex].size () > 1 ) tableCell ( os, fractionNames [searchIndex][specID->getFraction ()-1], false, false, styleID );
		bool sPlates = joint ? spottingPlatesFlag : spottingPlates [searchIndex];
		if ( reportLinks )
			specID->printTableCell2 ( os, sPlates, sPlates, styleID, smtl, searchKey [searchIndex] );
		else
			specID->printTableCell ( os, sPlates, sPlates, styleID );
	}
}
void PeptidePosition::printHeaderDelimited ( ostream& os, int searchNumber, bool showTimes )
{
	if ( reportSearchNumber ) delimitedHeader ( os, "Search #" );
	if ( showTimes ) printHeaderDelimitedPeak1 ( os );
	if ( reportMPlusHCalc )	delimitedHeader ( os, "M+H Calc" );
	if ( reportMOverZCalc )	delimitedHeader ( os, "m/z Calc" );
	if ( reportError ) {
		if ( searchNumber == -1 ) {
			if ( multipleErrorUnits ) {
				delimitedHeader ( os, "Error" );
				delimitedHeader ( os, "Units" );
			}
			else
				delimitedHeader ( os, parentTolerances [0]->getUnitsString ()  );
		}
		else
			delimitedHeader ( os, parentTolerances [searchNumber]->getUnitsString ()  );
	}
	if ( reportPreviousAA )	delimitedHeader ( os, "Prev AA" );
	if ( reportDBPeptide )	delimitedHeader ( os, "DB Peptide" );
	if ( peptideModType == "Mods In Peptide" ) delimitedHeader ( os, "Peptide" );
	if ( reportNextAA )		delimitedHeader ( os, "Next AA" );
	if ( peptideModType == "Variable Mods Only" )	delimitedHeader ( os, "Variable Mods" );
	else if ( peptideModType == "Constant Mods Only" )	delimitedHeader ( os, "Constant Mods" );
	else if ( peptideModType == "All Mods (1 column)" )	delimitedHeader ( os, "Mods" );
	else if ( peptideModType == "All Mods (2 columns)" ) {
		delimitedHeader ( os, "Constant Mods" );
		delimitedHeader ( os, "Variable Mods" );
	}
	if ( reportProteinMods )delimitedHeader ( os, "Protein Mods" );
	if ( reportLength )		delimitedHeader ( os, "Length" );
	if ( reportComposition )delimitedHeader ( os, "Composition" );
	if ( reportMModValue )	delimitedHeader ( os, "M Mod" );
	if ( reportMissedCleavages )delimitedHeader ( os, "M Cl" );
	if ( showTimes ) printHeaderDelimitedPeak2 ( os, searchNumber );
	if ( reportMSMSInfo )	SpecID::printDelimitedMSMSInfoHeader ( os );
	PeptidePositionQuan::printHeaderDelimited ( os, searchNumber );
	if ( reportStartAA )	delimitedHeader ( os, "Start" );
	if ( reportEndAA )		delimitedHeader ( os, "End" );
	if ( reportElemComp )	delimitedHeader ( os, "Elemental Composition" );
}
void PeptidePosition::printHeaderDelimitedPeak1 ( ostream& os )
{
	if ( reportMPlusH )		delimitedHeader ( os, "M+H" );
	if ( reportMOverZ )		delimitedHeader ( os, "m/z" );
	if ( reportCharge )		delimitedHeader ( os, "z" );
	if ( reportIntensity )	delimitedHeader ( os, "Intensity" );
}
void PeptidePosition::printHeaderDelimitedPeak2 ( ostream& os, int searchNumber )
{
	if ( reportTime ) {
		delimitedHeader ( os, "Fraction" );
		bool sPlates = ( searchNumber == -1 ) ? spottingPlatesFlag : spottingPlates [searchNumber];
		bool sNum = ( searchNumber == -1 ) ? spectrumNumberFlag : spectrumNumber [searchNumber];
		SpecID::printDelimitedHeader ( os, sPlates, sPlates || sNum );
	}
}
void PeptidePosition::printDelimited ( ostream& os, const ProteinInfo& proteinInfo, bool joint, bool showTimes ) const
{
	if ( reportSearchNumber ) delimitedCell ( os, searchIndex+1 );
	if ( showTimes ) printDelimitedPeak1 ( os );
	if ( getPeptide ().empty () ) {
		if ( reportMPlusHCalc )	delimitedEmptyCell ( os );
		if ( reportMOverZCalc )	delimitedEmptyCell ( os );
		if ( reportError ) {
			if ( joint && multipleErrorUnits )
				delimitedEmptyNCells ( os, 2 );
			else
				delimitedEmptyCell ( os );
		}
		if ( reportPreviousAA )	delimitedEmptyCell ( os );
		if ( reportDBPeptide )	delimitedCell ( os, "Pks: " + gen_itoa ( getNumPeaks () ) );
		if ( reportNextAA )		delimitedEmptyCell ( os );
		if ( peptideModType != "Off" ) {
			if ( peptideModType == "All Mods (2 columns)" )	delimitedEmptyNCells ( os, 2 );
			else											delimitedEmptyCell ( os );
		}
		if ( reportProteinMods )delimitedEmptyCell ( os );
		if ( reportLength )		delimitedCell ( os, "---" );
		if ( reportComposition )delimitedCell ( os, "---" );
		if ( reportMModValue )	delimitedCell ( os, "---" );
		if ( reportMissedCleavages )delimitedCell ( os, "---" );
	}
	else {
		if ( reportMPlusHCalc )	delimitedCell ( os, getMPlusHCalc ( searchIndex ), 4 );
		if ( reportMOverZCalc )	delimitedCell ( os, getMOverZCalc ( searchIndex ), 4 );
		if ( reportError ) {
			delimitedCellSigFig ( os, error, 2 );
			if ( joint && multipleErrorUnits ) delimitedCell ( os, parentTolerances [searchIndex]->getUnitsString () );
		}
		if ( reportPreviousAA ) delimitedCell ( os, proteinInfo.getPreviousAA ( startAA, reportPreviousAA ) );
		if ( reportDBPeptide )	delimitedCell ( os, getDBPeptide () );
		if ( peptideModType == "Mods In Peptide" ) {
			ostringstream ost;
			if ( !getNTerm ().empty () ) ost << getNTerm () << '-';
			ost << getPeptide ();
			if ( !getCTerm ().empty () ) ost << '-' << getCTerm ();
			if ( !getNeutralLoss ().empty () ) ost << '+' << getNeutralLoss ();
			delimitedCell ( os, ost.str () );
		}
		if ( reportNextAA )		delimitedCell ( os, proteinInfo.getNextAA ( getEndAA (), reportNextAA ) );
		if ( peptideModType != "Off" && peptideModType != "Mods In Peptide" ) {
			string cMod = SCModInfo::getConstModsString ( getDBPeptide (), searchIndex );
			string vMod = SCModInfo::getModsString ( searchIndex, specID, getDBPeptide () );
			if ( peptideModType == "Variable Mods Only" )		delimitedCell ( os, vMod );
			else if ( peptideModType == "Constant Mods Only" )	delimitedCell ( os, cMod );
			else if ( peptideModType == "All Mods (1 column)" )	{
				delimitedCell ( os, cMod.empty () ? vMod : ( vMod.empty () ? cMod : cMod + ';' + vMod ) );
			}
			else if ( peptideModType == "All Mods (2 columns)" ) {
				delimitedCell ( os, cMod );
				delimitedCell ( os, vMod );
			}
		}
		if ( reportProteinMods )delimitedCell ( os, SCModInfo::getModsString ( searchIndex, specID, getDBPeptide (), startAA ) );
		if ( reportLength )		delimitedCell ( os, getPeptideLength () );
		if ( reportComposition )delimitedCell ( os, checkComposition () );
		if ( reportMModValue ) {
			double m = getMModValue ();
			if ( m )delimitedCell ( os, m, 0 );
			else	delimitedEmptyCell ( os );
		}
		if ( reportMissedCleavages )delimitedCell ( os, getMissedCleavages () );
	}
	if ( showTimes ) printDelimitedPeak2 ( os, joint );
	if ( reportMSMSInfo ) specID->printDelimitedMSMSInfoCell ( os );
	if ( getPeptide ().empty () || quanRatio == 0 ) {
		PeptidePositionQuan::printQuanBlankDelimited ( os, searchIndex );
	}
	else {
		PeptidePositionQuan::printDelimited ( os, quanRatio, joint );
	}
	if ( getPeptide ().empty () ) {
		if ( reportStartAA )delimitedCell ( os, "---" );
		if ( reportEndAA )	delimitedCell ( os, "---" );
	}
	else {
		if ( reportStartAA )delimitedCell ( os, startAA );
		if ( reportEndAA )	delimitedCell ( os, getEndAA () );
	}
	if ( reportElemComp )	delimitedCell ( os, getElemComp () );
}
void PeptidePosition::printDelimitedPeak1 ( ostream& os ) const
{
	if ( reportMPlusH )		delimitedCell ( os, getMPlusH (), 4 );
	if ( reportMOverZ )		delimitedCell ( os, getMOverZ (), 4 );
	if ( reportCharge )		delimitedCell ( os, getCharge () );
	if ( reportIntensity )	delimitedCellSigFig ( os, getIntensity (), 3 );
}
void PeptidePosition::printDelimitedPeak2 ( ostream& os, bool joint ) const
{
	if ( reportTime ) {
		delimitedCell ( os, fractionNames [searchIndex][specID->getFraction ()-1] );
		bool sPlates = joint ? spottingPlatesFlag : spottingPlates [searchIndex];
		bool sNum = joint ? spectrumNumberFlag : spectrumNumber [searchIndex];
		specID->printDelimitedCell ( os, sPlates, sPlates || sNum );
	}
}
bool PeptidePosition::outputQuanResults ( ostream& os, const string& searchName, int numRepeats, bool area ) const
{
	return PeptidePositionQuan::outputQuanResults ( os, quanRatio, searchName, numRepeats, area );
}
DoubleVector PeptidePosition::getIntensityRatios () const
{
	return PeptidePositionQuan::getIntensityRatios ( quanRatio );
}
DoubleVector PeptidePosition::getAreaRatios () const
{
	return PeptidePositionQuan::getAreaRatios ( quanRatio );
}
void PeptidePosition::initialiseComposition ( const StringVector& compIons, const StringVector& massCompIons, const string& maskType )
{
	compMask = 0;
	unsigned int compM = 0;
	StringVector comps;
	copy ( compIons.begin (), compIons.end (), back_inserter ( comps ) );
	copy ( massCompIons.begin (), massCompIons.end (), back_inserter ( comps ) );
	for ( StringVectorSizeType i = 0 ; i < comps.size () ; i++ ) {
		if ( i == 32 ) {
			ErrorHandler::genError ()->error ( "A maximum of 32 composition options can be selected at one time.\n" );
		}
		compM = ( i == 0 ) ? 1 : compM * 2;
		string cIon = comps [i];
		string::size_type cIonLen = cIon.length ();
		string s;
		if ( cIon [cIonLen-1] == ')' ) {
			int pos = cIon.rfind ( '(' );
			string theAA = cIon.substr ( pos+1, cIonLen-pos-2 );
			string theMod = cIon.substr ( 0, pos-1 );
			if ( theAA == "N-term" ) {
				s = "n" + theMod;
			}
			else if ( theAA == "C-term" ) {
				s = "c" + theMod;
			}
			else if ( theAA == "Neutral loss" ) {
				s = "." + theMod;
			}
			else {
				s = theAA + '(' + theMod + ')';
			}
		}
		else {	// This is an amino acid
			s = cIon;
		}
		compMaskMap [s] = compM;
		compMask |= compM;
	}
	compMaskType = maskType;
}
void PeptidePosition::initialiseParams ( const ParameterList* p )
{
	static string enzyme;
	if ( params.empty () ) {	// First parameter set
		enzyme = p->getStringValue ( "enzyme" );
		if ( enzyme != "No enzyme" ) {
			init_fasta_enzyme_function ( enzyme );
			enzymeInit = true;
		}
	}
	else {						// Subsequent parameter set
		if ( p->getStringValue ( "enzyme" ) != enzyme ) {
			enzymeInit = false;
		}
	}
	params.push_back ( p );
	instrument.push_back ( p->getStringValue ( "instrument_name" ) );
	InstrumentInfo instInfo ( instrument.back () );
	chargeReducedFragmentation.push_back ( instInfo.getChargeReducedFragmentation () );
	searchKey.push_back ( p->getStringValue ( "search_key" ) );
	PeptidePositionQuan::initialiseParams ( p );
}
void PeptidePosition::initialiseDataSetInfo ( double timeWindowStart, double timeWindowEnd )
{
	rtIntervalStart = timeWindowStart;
	rtIntervalEnd = timeWindowEnd;
	PeptidePositionQuan::initialiseDataSetInfo ( timeWindowStart, timeWindowEnd );
}
int PeptidePosition::getNumQuanRatioColumns ()
{
	return PeptidePositionQuan::getNumQuanRatioColumns ();
}
void PeptidePosition::runMSProduct ( int i, double score ) const
{
	if ( !hasHitPeptide () ) return;
	static int idx = 0;
	idx++;
	string scoreStr = gen_ftoa ( score, "%.1f" );
	string outputFilename =  gen_itoa ( idx ) + "_" + stripFilenameChars ( getPeptide () ) + "_" + scoreStr + "_" + gen_itoa ( getCharge () );
	PPTempFile tempFile ( "msprod", ".txt" );
	string tempFileFullPath = tempFile.getFullPath ();

	string command;
#ifdef VIS_C
	command += "\"";
#endif
	command += getSystemCall ( "mssearch.cgi" );
	command += " - ";
	command += getCommandLineNVPair ( "search_name", "msproduct" );
	command += " ";
	command += getCommandLineNVPair ( "output_filename", outputFilename );
	command += " ";
	command += getCommandLineNVPair ( "results_to_file", "1" );
	command += " ";
	command += getCommandLineNVPair ( "output_type", "XML" );
	command += " ";
	command += getCommandLineNVPair ( "report_title", "MS-Product" );
	command += " ";
	command += getCommandLineNVPair ( "version", Version::instance ().getVersion () );
	command += " ";
	command += getCommandLineNVPair ( "data_source", "List of Files" );
	command += " ";
	command += getCommandLineNVPair ( "use_instrument_ion_types", "1" );
	command += " ";
	command += getCommandLineNVPair ( "search_key", searchKey [i] );
	command += " ";
	command += getCommandLineNVPair ( "instrument_name", instrument [i] );
	command += " ";
	command += parentTolerances [i]->getCommandLineNVPair ( "msms_parent_mass" );
	command += fragmentTolerances [i]->getCommandLineNVPair ( "fragment_masses" );
	command += getCommandLineNVPair ( "parent_mass_convert", "monoisotopic" );
	command += " ";
	command += specID->getCommandLineNVPair ();
	command += getCommandLineNVPair ( "max_charge", mmsi->getCharge () );
	command += " ";
	command += hitPeptide->getCommandLineNVPair ( 1 );
	command += MSMSPeakFilterOptions::getCommandLineNVPair ( ProgramLink::getParams () );
	command += "> ";
#ifdef VIS_C
	command += "\"" + tempFileFullPath + "\"";
#else
	command += tempFileFullPath;
#endif
#ifdef VIS_C
	command +=  "\"";
#endif
	genSystem ( command, "", true );
	genUnlink ( tempFileFullPath );
}
void PeptidePosition::initialise ( const SearchResultsPtrVector& searchResults )
{
	SetString errorUnitsSet;
	SetString fractionNamesSet;
	for ( SearchResultsPtrVectorSizeType i = 0 ; i < searchResults.size () ; i++ ) {
		StringVector fn = searchResults [i]->getFractionNames ();
		fractionNames.push_back ( fn );
		for ( StringVectorConstIterator ii = fn.begin () ; ii != fn.end () ; ii++ ) {
			fractionNamesSet.insert ( *ii );
		}

		spottingPlates.push_back ( searchResults [i]->getSpottingPlate () );
		if ( spottingPlates.back () ) spottingPlatesFlag = true;

		rawTypes.push_back ( searchResults [i]->getRawTypes () );
		if ( rawTypes.back () [0] != WIFF ) quanMultiNormalFlag = true;

		parentTolerances.push_back ( searchResults [i]->getParentTolerance () );
		errorUnitsSet.insert ( parentTolerances.back ()->getUnitsString () );

		sysErrorStr.push_back ( searchResults [i]->getSysErrorStr () );
		offsets.push_back ( searchResults [i]->getOffsets () );
		fragmentTolerances.push_back ( searchResults [i]->getFragmentTolerance () );
	}
	multipleErrorUnits = errorUnitsSet.size () > 1;
	multipleFractionNames = fractionNamesSet.size () > 1;
	PeptidePositionQuan::initialise ( searchResults );
}
void PeptidePosition::addSiteScores ( SiteScores& siteScores, int line ) const
{
	siteScores.add ( getDBPeptide (), SCModInfo::getModsString ( searchIndex, specID, getDBPeptide (), startAA ), startAA, line );	// peptide, mods, start, tenLogP
}
void PeptidePosition::initialiseAACalculator ( const MapStringConstModPtr& constMods )
{
	aaCalc = new AACalculator ( true, constMods );
}
bool SearchResultsProteinInfo::reportScore = false;
bool SearchResultsProteinInfo::reportNumUnique = false;
bool SearchResultsProteinInfo::reportPeptideCount = false;
bool SearchResultsProteinInfo::reportBestScore = false;
bool SearchResultsProteinInfo::reportCoverage = false;
bool SearchResultsProteinInfo::reportBestDiscScore = false;
bool SearchResultsProteinInfo::reportBestExpectVal = false;

int SearchResultsProteinInfo::getColspan () const
{
	int colspan = 0;
	if ( reportScore )			colspan++;
	if ( reportNumUnique )		colspan++;
	if ( reportPeptideCount )	colspan++;
	if ( reportBestScore )		colspan++;
	return colspan;
}
void SearchResultsProteinInfo::printHeaderHTML ( ostream& os, int index, const string& styleID ) const
{
	if ( reportScore )			tableHeader ( os, "Protein Score", styleID );
	if ( reportNumUnique )		tableHeader ( os, "Num Unique", styleID );
	if ( reportPeptideCount )	tableHeader ( os, "Peptide Count", styleID );
	if ( reportBestScore )		tableHeader ( os, "Best Score", styleID );
}
void SearchResultsProteinInfo::printHTML ( ostream& os, bool empty, const string& styleID ) const
{
	if ( score == 0.0 || empty ) {
		if ( reportScore )			tableEmptyCell ( os, styleID );
		if ( reportNumUnique )		tableEmptyCell ( os, styleID );
		if ( reportPeptideCount )	tableEmptyCell ( os, styleID );
		if ( reportBestScore )		tableEmptyCell ( os, styleID );
	}
	else {
		if ( reportScore )			tableCell ( os, score, 1, false, styleID );
		if ( reportNumUnique )		tableCell ( os, numUnique, false, styleID );
		if ( reportPeptideCount )	tableCell ( os, peptideCount, false, styleID );
		if ( reportBestScore )		tableCell ( os, bestScore, 1, false, styleID );
	}
}
void SearchResultsProteinInfo::printHeaderDelimited ( ostream& os, int index ) const
{
	if ( reportScore )			delimitedHeader ( os, "Protein Score" );
	if ( reportNumUnique )		delimitedHeader ( os, "Num Unique" );
	if ( reportPeptideCount )	delimitedHeader ( os, "Peptide Count" );
	if ( reportBestScore )		delimitedHeader ( os, "Best Score" );
}
void SearchResultsProteinInfo::printDelimited ( ostream& os ) const
{
	if ( score == 0.0 ) {
		printDelimitedEmpty ( os );
	}
	else {
		if ( reportScore )			delimitedCell ( os, score, 1 );
		if ( reportNumUnique )		delimitedCell ( os, numUnique );
		if ( reportPeptideCount )	delimitedCell ( os, peptideCount );
		if ( reportBestScore )		delimitedCell ( os, bestScore, 1 );
	}
}
void SearchResultsProteinInfo::printDelimitedEmpty ( ostream& os ) const
{
	if ( reportScore )			delimitedEmptyCell ( os );
	if ( reportNumUnique )		delimitedEmptyCell ( os );
	if ( reportPeptideCount )	delimitedEmptyCell ( os );
	if ( reportBestScore )		delimitedEmptyCell ( os );
}
SearchResultsProteinInfo::~SearchResultsProteinInfo () {}
void SearchResultsProteinInfo::setParams ( const ParameterList* params )
{
	reportScore			= params->getBoolValue ( "report_prot_score" );
	reportNumUnique		= params->getBoolValue ( "report_num_unique" );
	reportPeptideCount	= params->getBoolValue ( "report_peptide_count" );
	reportBestScore		= params->getBoolValue ( "report_best_score" );
	reportCoverage		= params->getBoolValue ( "report_coverage" );
	reportBestDiscScore	= params->getBoolValue ( "report_best_disc_score" );
	reportBestExpectVal	= params->getBoolValue ( "report_best_expect" );
}
bool PPPeptideHitInfo::reportUnmatched = false;
bool PPPeptideHitInfo::reportNumPks = false;
bool PPPeptideHitInfo::reportRank = false;
bool PPPeptideHitInfo::reportScore = false;
bool PPPeptideHitInfo::reportScoreDiff = false;
bool PPPeptideHitInfo::reportExpectation = false;
bool PPPeptideHitInfo::reportPValue = false;
bool PPPeptideHitInfo::reportMValue = false;
bool PPPeptideHitInfo::reportNumPrecursor = false;
bool PPPeptideHitInfo::reportGradient = false;
bool PPPeptideHitInfo::reportOffset = false;
bool PPPeptideHitInfo::reportDiscScore = false;
bool PPPeptideHitInfo::reportRepeats = false;
SearchResultsProteinHit::SearchResultsProteinHit ( const SearchResultsProteinInfo* proteinHitInfo, const string& accessionNumber ) :
	proteinHitInfo ( proteinHitInfo ),
	accessionNumber ( accessionNumber )
{
}
SearchResultsProteinHit::~SearchResultsProteinHit () {}

StringVector SearchResults::getAccessionNumbers ( const string& id )
{
	StringVector sv;
	const SearchResultsProteinHitPtrVector& pHits = proteinHits[id];
	for ( SearchResultsProteinHitPtrVectorSizeType i = 0 ; i < pHits.size () ; i++ ) {
		sv.push_back ( pHits[i]->getAccessionNumber () );
	}
	stable_sort ( sv.begin (), sv.end () );
	sv.erase ( unique ( sv.begin (), sv.end () ), sv.end () );
	return sv;
}
const SearchResultsProteinInfo* SearchResults::getSearchResultsProteinInfo ( const string& aNum, const string& id )
{
	if ( aNum.empty () ) return emptyProteinHit;
	MapAccessionNumberToSearchResultsProteinInfoConstIterator cur = proteinMap [id].find ( aNum );
	if ( cur != proteinMap [id].end () ) {
		return (*cur).second;
	}
	return emptyProteinHit;
}
const SearchResultsProteinInfo* SearchResults::getSearchResultsProteinInfo2 ( const string& aNum, const string& id )
{
	MapAccNumSearchResultsProteinHitPtrConstIterator cur = compProtHits [id].find ( aNum );
	if ( cur != compProtHits [id].end () ) {
		return (*cur).second->getProteinHitInfo ();
	}
	return 0;
}
void SearchResults::setProteinMap ( const string& id )
{
	for ( SearchResultsProteinHitPtrVectorSizeType i = 0 ; i < proteinHits[id].size () ; i++ ) {
		proteinMap [id][proteinHits [id][i]->getAccessionNumber ()] = proteinHits [id][i]->getProteinHitInfo ();
	}
}
void SearchResults::setPeptideHitsAndMap ( const string& id )
{
	peptideHits [id] = 0;
	SearchResultsProteinHitPtrVector& pHits = proteinHits [id];
	for ( SearchResultsProteinHitPtrVectorSizeType a = 0 ; a < pHits.size () ; a++ ) {
		SearchResultsProteinHit* mmxph = pHits [a];
		for ( int j = 0 ; j < mmxph->getNumMSMSHits () ; j++ ) {
			SearchResultsPeptideHit* srph = mmxph->getMSMSHit ( j );
			const PeptidePosition* pp = srph->getPeptidePosition ();
			const PPPeptideHitInfo* sri = srph->getPeptideHitInfo ();
			SetSearchResultsPeptideHitConstIterator cur = peptideMap [id].find ( srph );
			if ( cur == peptideMap [id].end () || sri->getScore () > (*cur)->getScore () ) {
				peptideMap [id].insert ( srph );
			}
			peptideHits [id]++;
		}
	}
}
PairSearchResultsPeptideHitPtrVectorConstIterator SearchResults::getSearchPeptidePositionIters ( const string& accessionNumber, const string& id )
{
	MapAccNumSearchResultsPeptideHitPtrVectorConstIterator cur = compPepHits [id].find ( accessionNumber );

	if ( cur == compPepHits [id].end () ) {
		ErrorHandler::genError ()->message ( "Comp hits not found" );
	}
	SearchResultsPeptideHitPtrVectorConstIterator i1 = (*cur).second.begin ();
	SearchResultsPeptideHitPtrVectorConstIterator i2 = (*cur).second.end ();
	return make_pair ( i1, i2 );
}
SearchResultsPeptideHit* SearchResults::getPPPeptideHitInfoPair ( SearchResultsPeptideHit* srph, const string& id )
{
	if ( srph == 0 ) return 0;
	SetSearchResultsPeptideHitConstIterator cur = peptideMap [id].find ( srph );
	if ( cur != peptideMap [id].end () ) {
		return (*cur);
	}
	return 0;
}
StringVector SearchResults::getIDList ()
{
	StringVector sv;
	for ( SetStringConstIterator i = idSet.begin () ; i != idSet.end () ; i++ ) {
		sv.push_back ( *i );
	}
	return sv;
}
const double PeptideSpectralInfo::NO_SCORE_DIFF = std::numeric_limits<double>::infinity();
const double PeptideSpectralInfo::INVALID_MASCOT_SCORE = std::numeric_limits<double>::infinity();
bool PeptideSpectralInfo::linearTailFitExpectation = true;
bool PeptideSpectralInfo::noExpectation = false;
PeptideSpectralInfo::PeptideSpectralInfo ( int unmatched, int numPeaks, int rank, double score, double expectation, const string& scDiff, int numPrecursor, double a, double b ) :
	unmatched ( unmatched ),
	numPeaks ( numPeaks ),
	rank ( rank ),
	score ( score ),
	scoreDiff ( scDiff == "---" ? NO_SCORE_DIFF : atof ( scDiff.c_str () ) ),
	expectation ( expectation ),
	numPrecursor ( numPrecursor ),
	a ( a ),
	b ( b )
{
}
void PeptideSpectralInfo::setExpectationFlag ( const string& expectationCalculationType )
{
	noExpectation = ( expectationCalculationType == "None" );
	linearTailFitExpectation = ( expectationCalculationType == "Linear Tail Fit" );
}
double PeptideSpectralInfo::getExpectationValue ( double score ) const
{
	double e;
	if ( ( a == 0 && b == 0 ) || noExpectation ) {
		e = -1.0;
	}
	else {
		if ( linearTailFitExpectation )
			e = pow ( 10.0, (a * score) + b );
		else
			e = 1.0 - exp ( -exp ( (- 1 / b) * (score - a) ) );	// Method of moments
		e *= numPrecursor;
		if ( isnan ( e ) ) e = -1.0;
	}
	return e;
}
double PeptideSpectralInfo::getPValue ( double score ) const
{
	double e = getExpectationValue ( score );
	if ( e == -1.0 ) return -1.0;		// PValue invalid
	else return e / numPrecursor;
}
double PeptideSpectralInfo::getMascotScore ( double score ) const
{
	double pv = getPValue ( score );
	if ( pv == -1.0 ) return INVALID_MASCOT_SCORE;
	else return - 10.0 * log10 ( pv );
}
const double PPPeptideHitInfo::NO_SCORE_DIFF = std::numeric_limits<double>::infinity();
PPPeptideHitInfo::PPPeptideHitInfo () :
	psi ( new PeptideSpectralInfo ( 0, 0, 0, 0.0, std::numeric_limits<double>::infinity(), "---", 0 ) )
{
}
PPPeptideHitInfo::PPPeptideHitInfo ( const PeptideSpectralInfo* psi ) :
	psi ( psi )
{
}
PPPeptideHitInfo::~PPPeptideHitInfo () {}
int PPPeptideHitInfo::getColspan () const
{
	int colspan = 0;
	if ( reportUnmatched )	colspan++;
	if ( reportNumPks )		colspan++;
	if ( reportRank )		colspan++;
	if ( reportScore )		colspan++;
	if ( reportScoreDiff )	colspan++;
	if ( reportExpectation )colspan++;		// Expectation
	if ( reportPValue )		colspan++;
	if ( reportMValue )		colspan++;
	if ( reportNumPrecursor )colspan++;
	if ( reportGradient )	colspan++;
	if ( reportOffset )		colspan++;
	if ( reportDiscScore )	colspan++;
	if ( reportRepeats )	colspan++;
	return colspan;
}
void PPPeptideHitInfo::printHeaderHTML ( ostream& os, const string& styleID )
{
	if ( reportUnmatched )	tableHeader ( os, "# Unmat", styleID );
	if ( reportNumPks )		tableHeader ( os, "# Pks", styleID );
	if ( reportRank )		tableHeader ( os, "Rank", styleID );
	if ( reportScore )		tableHeader ( os, "Score", styleID );
	if ( reportScoreDiff )	tableHeader ( os, "Score Diff", styleID );
	if ( reportExpectation )tableHeader ( os, "Expect", styleID );
	if ( reportPValue )		tableHeader ( os, "P Value", styleID );
	if ( reportMValue )		tableHeader ( os, "-10logP", styleID );
	if ( reportNumPrecursor )tableHeader ( os, "# Precursor", styleID );
	if ( reportGradient )	tableHeader ( os, "Gradient", styleID );
	if ( reportOffset )		tableHeader ( os, "Offset", styleID );
	if ( reportDiscScore )	tableHeader ( os, "Disc Score", styleID );
	if ( reportRepeats )	tableHeader ( os, "# in DB", styleID );
}
void PPPeptideHitInfo::printHTML ( ostream& os, const string& styleID, int colspan, int rowspan ) const
{
	if ( psi->getScore () != 0.0 ) {
		if ( reportUnmatched )	tableCell ( os, psi->getUnmatched (), false, styleID, colspan, rowspan );
		if ( reportNumPks )		tableCell ( os, psi->getNumPeaks (), false, styleID, colspan, rowspan );
		if ( reportRank )		tableCell ( os, psi->getRank (), false, styleID, colspan, rowspan );
		if ( reportScore )		tableCell ( os, psi->getScore (), 1, false, styleID, colspan, rowspan );
		if ( reportScoreDiff ) {
			if ( psi->getScoreDifference () == NO_SCORE_DIFF )
				tableCell ( os, "---", false, false, styleID, colspan, rowspan );
			else
				tableCell ( os, psi->getScoreDifference (), 1, false, styleID, colspan, rowspan );
		}
		if ( reportExpectation )tableCellSigFig ( os, psi->getExpectation (), 2, false, styleID, colspan, rowspan );
		if ( reportPValue )		tableCellSigFig ( os, psi->getPValue (), 2, false, styleID, colspan, rowspan );
		if ( reportMValue )		tableCellSigFig ( os, psi->getMascotScore (), 0, false, styleID, colspan, rowspan );
		if ( reportNumPrecursor )tableCell ( os, psi->getNumPrecursor (), false, styleID, colspan, rowspan );
		if ( reportGradient )	tableCellSigFig ( os, psi->getA (), 5, false, styleID, colspan, rowspan );
		if ( reportOffset )		tableCellSigFig ( os, psi->getB (), 5, false, styleID, colspan, rowspan );
		if ( reportDiscScore ) {
			if ( discriminantScore == DiscriminantScore::MIN_DISC_SCORE )
				tableCell ( os, "---", false, false, styleID, colspan, rowspan );
			else
				tableCell ( os, discriminantScore, 2, false, styleID, colspan, rowspan );
		}
		if ( reportRepeats )	tableCell ( os, repeats, false, styleID, colspan, rowspan );
	}
	else {
		if ( reportUnmatched )	tableEmptyCell ( os, styleID, colspan, rowspan );
		if ( reportNumPks )		tableEmptyCell ( os, styleID, colspan, rowspan );
		if ( reportRank )		tableEmptyCell ( os, styleID, colspan, rowspan );
		if ( reportScore )		tableEmptyCell ( os, styleID, colspan, rowspan );
		if ( reportScoreDiff )	tableEmptyCell ( os, styleID, colspan, rowspan );
		if ( reportExpectation )tableEmptyCell ( os, styleID, colspan, rowspan );	// Expectation
		if ( reportPValue )		tableEmptyCell ( os, styleID, colspan, rowspan );
		if ( reportMValue )		tableEmptyCell ( os, styleID, colspan, rowspan );
		if ( reportNumPrecursor )tableEmptyCell ( os, styleID, colspan, rowspan );
		if ( reportGradient )	tableEmptyCell ( os, styleID, colspan, rowspan );
		if ( reportOffset )		tableEmptyCell ( os, styleID, colspan, rowspan );
		if ( reportDiscScore )	tableEmptyCell ( os, styleID, colspan, rowspan );
		if ( reportRepeats )	tableEmptyCell ( os, styleID, colspan, rowspan );
	}
}
void PPPeptideHitInfo::printHeaderDelimited ( ostream& os )
{
	if ( reportUnmatched )	delimitedHeader ( os, "# Unmat" );
	if ( reportNumPks )		delimitedHeader ( os, "# Pks" );
	if ( reportRank )		delimitedHeader ( os, "Rank" );
	if ( reportScore )		delimitedHeader ( os, "Score" );
	if ( reportScoreDiff )	delimitedHeader ( os, "Score Diff" );
	if ( reportExpectation )delimitedHeader ( os, "Expect" );
	if ( reportPValue )		delimitedHeader ( os, "P Value" );
	if ( reportMValue )		delimitedHeader ( os, "-10logP" );
	if ( reportNumPrecursor )delimitedHeader ( os, "# Precursor" );
	if ( reportGradient )	delimitedHeader ( os, "Gradient" );
	if ( reportOffset )		delimitedHeader ( os, "Offset" );
	if ( reportDiscScore )	delimitedHeader ( os, "Disc Score" );
	if ( reportRepeats )	delimitedHeader ( os, "# in DB" );
}
void PPPeptideHitInfo::printDelimited ( ostream& os ) const
{
	if ( psi->getScore () != 0.0 ) {
		if ( reportUnmatched )	delimitedCell ( os, psi->getUnmatched () );
		if ( reportNumPks )		delimitedCell ( os, psi->getNumPeaks () );
		if ( reportRank )		delimitedCell ( os, psi->getRank () );
		if ( reportScore )		delimitedCell ( os, psi->getScore (), 1 );
		if ( reportScoreDiff ) {
			if ( psi->getScoreDifference () == NO_SCORE_DIFF )
				delimitedCell ( os, "---" );
			else
				delimitedCell ( os, psi->getScoreDifference (), 1 );
		}
		if ( reportExpectation )delimitedCellSigFig ( os, psi->getExpectation (), 2 );
		if ( reportPValue )		delimitedCellSigFig ( os, psi->getPValue (), 2 );
		if ( reportMValue )		delimitedCellSigFig ( os, psi->getMascotScore (), 0 );
		if ( reportNumPrecursor )delimitedCell ( os, psi->getNumPrecursor () );
		if ( reportGradient )	delimitedCellSigFig ( os, psi->getA (), 5 );
		if ( reportOffset )		delimitedCellSigFig ( os, psi->getB (), 5 );
		if ( reportDiscScore ) {
			if ( discriminantScore == DiscriminantScore::MIN_DISC_SCORE )
				delimitedCell ( os, "---" );
			else
				delimitedCell ( os, discriminantScore, 2 );
		}
		if ( reportRepeats )	delimitedCell ( os, repeats );
	}
	else {
		if ( reportUnmatched )	delimitedEmptyCell ( os );
		if ( reportNumPks )		delimitedEmptyCell ( os );
		if ( reportRank )		delimitedEmptyCell ( os );
		if ( reportScore )		delimitedEmptyCell ( os );
		if ( reportScoreDiff )	delimitedEmptyCell ( os );
		if ( reportExpectation )delimitedEmptyCell ( os );
		if ( reportPValue )		delimitedEmptyCell ( os );
		if ( reportMValue )		delimitedEmptyCell ( os );
		if ( reportNumPrecursor )delimitedEmptyCell ( os );
		if ( reportGradient )	delimitedEmptyCell ( os );
		if ( reportOffset )		delimitedEmptyCell ( os );
		if ( reportDiscScore )	delimitedEmptyCell ( os );
		if ( reportRepeats )	delimitedEmptyCell ( os );
	}
}
SearchResultsPeptideHit::SearchResultsPeptideHit ( const SpecID* spID, const MSMSSpectrumInfo* mmsi, int searchIndex ) :
	peptidePosition ( spID, mmsi, searchIndex )
{
}
SearchResultsPeptideHit::SearchResultsPeptideHit ( const SpecID* specID, const MSMSSpectrumInfo* mmsi, const PeptideSpectralInfo* psi, double error, const HitPeptide* hitPeptide, int startAA, const string& accNum, int searchIndex ) :
	peptideHitInfo ( psi ),
	peptidePosition ( accNum, hitPeptide, error, specID, mmsi, startAA, searchIndex )
{
}
void SearchResultsPeptideHit::printHTML ( ostream& os, const ProteinInfo& proteinInfo, int searchIndex, const string& styleID, const SCMSTagLink& smtl, const string& id, bool joint, bool showTimes ) const
{
	if ( hasSpecID () ) {
		peptidePosition.printHTML ( os, proteinInfo, styleID, smtl, id, joint, showTimes );
		if ( PeptidePosition::getRunMSProductFlag () ) peptidePosition.runMSProduct ( searchIndex, peptideHitInfo.getScore () );
	}
	else
		PeptidePosition::printHTMLEmpty ( os, searchIndex, styleID );
	peptideHitInfo.printHTML ( os, styleID );
}
void SearchResultsPeptideHit::printHTMLPeak1 ( ostream& os, const string& styleID ) const
{
	if ( hasSpecID () )
		peptidePosition.printHTMLPeak1 ( os, styleID );
	else
		PeptidePosition::printHTMLEmptyPeak1 ( os, styleID );
}
void SearchResultsPeptideHit::printHTMLPeak2 ( ostream& os, int searchIndex, const string& styleID, const SCMSTagLink& smtl, bool joint ) const
{
	if ( hasSpecID () )
		peptidePosition.printHTMLPeak2 ( os, styleID, smtl, joint );
	else
		PeptidePosition::printHTMLEmptyPeak2 ( os, searchIndex, styleID );
}
void SearchResultsPeptideHit::printDelimited ( ostream& os, const ProteinInfo& proteinInfo, int searchNumber ) const
{
	if ( hasSpecID () )
		peptidePosition.printDelimited ( os, proteinInfo, searchNumber == -1 );
	else
		PeptidePosition::printDelimitedEmpty ( os, searchNumber );
	peptideHitInfo.printDelimited ( os );
}
void SearchResultsPeptideHit::printDelimitedEmpty ( ostream& os, int searchNumber )
{
	static PPPeptideHitInfo emptyPeptideHit;
	PeptidePosition::printDelimitedEmpty ( os, searchNumber );
	emptyPeptideHit.printDelimited ( os );
}
void SearchResultsPeptideHit::printHTMLEmpty ( ostream& os, int searchNumber, const string& styleID )
{
	static PPPeptideHitInfo emptyPeptideHit;
	PeptidePosition::printHTMLEmpty ( os, searchNumber, styleID );
	emptyPeptideHit.printHTML ( os, styleID );
}
void SearchResultsPeptideHit::addSiteScores ( SiteScores& siteScores, int line ) const
{
	if ( hasSpecID () ) {
		peptidePosition.addSiteScores ( siteScores, line );
	}
}
bool PPProteinHitQuanInfo::reportArea = false;
bool PPProteinHitQuanInfo::reportIntensity = false;
bool PPProteinHitQuanInfo::reportMedian = false;
bool PPProteinHitQuanInfo::reportIQR = false;
bool PPProteinHitQuanInfo::reportMean = false;
bool PPProteinHitQuanInfo::reportStDev = false;
bool PPProteinHitQuanInfo::reportNum = false;

double PPProteinHitQuanInfo::numStdDev = 2.0;
string PPProteinHitQuanInfo::stdevMinusHTMLHeader;
string PPProteinHitQuanInfo::stdevPlusHTMLHeader;
string PPProteinHitQuanInfo::stdevMinusDelimHeader;
string PPProteinHitQuanInfo::stdevPlusDelimHeader;

void PPProteinHitQuanInfo::initQuan ()
{
	if ( QuantitationRatio::getQuanReport () ) {
		numQuanPeaks = PeptidePosition::getNumQuanRatioColumns ();
		DoubleVector dv ( numQuanPeaks, 0.0 );
		IntVector iv ( numQuanPeaks, 0 );
		if ( reportArea ) {
			medianAreaRatio = dv;
			lowQAreaRatio = dv;
			highQAreaRatio = dv;
			meanAreaRatio = dv;
			stDevAreaRatio = dv;
			numAreaRatio = iv;
		}
		if ( reportIntensity ) {
			medianIntensityRatio = dv;
			lowQIntensityRatio = dv;
			highQIntensityRatio = dv;
			meanIntensityRatio = dv;
			stDevIntensityRatio = dv;
			numIntensityRatio = iv;
		}
	}
}
int PPProteinHitQuanInfo::getColspan () const
{
	int colspan = 0;
	if ( reportArea ) {
		if ( reportMedian )	colspan += numQuanPeaks;
		if ( reportIQR )	colspan += 2 * numQuanPeaks;
		if ( reportMean )	colspan += numQuanPeaks;
		if ( reportStDev )	colspan += 2 * numQuanPeaks;
		if ( reportNum )	colspan += numQuanPeaks;
	}
	if ( reportIntensity ) {
		if ( reportMedian )	colspan += numQuanPeaks;
		if ( reportIQR )	colspan += 2 * numQuanPeaks;
		if ( reportMean )	colspan += numQuanPeaks;
		if ( reportStDev )	colspan += 2 * numQuanPeaks;
		if ( reportNum )	colspan += numQuanPeaks;
	}
	return colspan;
}
string PPProteinHitQuanInfo::getRatioStr ( int i ) const
{
	if ( PeptidePositionQuan::getQuanMSMSFlag () )	return QuantitationMulti::getProteinRatioString ( i );
	else											return QuantitationRatio::getRatioString ( i );
}
void PPProteinHitQuanInfo::printHeaderHTML ( ostream& os, const string& styleID ) const
{
	if ( reportArea ) {
		for ( DoubleVectorSizeType i = 0 ; i < numQuanPeaks ; i++ ) {
			if ( reportMedian ) tableHeader ( os, "Med " + getRatioStr ( i ) + " A", styleID );
			if ( reportIQR ) {
				tableHeader ( os, "Q1", styleID );
				tableHeader ( os, "Q3", styleID );
			}
			if ( reportMean )	tableHeader ( os, "Mn " + getRatioStr ( i ) + " A", styleID );
			if ( reportStDev ) {
				tableHeader ( os, stdevMinusHTMLHeader, styleID );
				tableHeader ( os, stdevPlusHTMLHeader, styleID );
			}
			if ( reportNum )	tableHeader ( os, "Num", styleID );
		}
	}
	if ( reportIntensity ) {
		for ( DoubleVectorSizeType i = 0 ; i < numQuanPeaks ; i++ ) {
			if ( reportMedian ) tableHeader ( os, "Med " + getRatioStr ( i ) + " I", styleID );
			if ( reportIQR ) {
				tableHeader ( os, "Q1", styleID );
				tableHeader ( os, "Q3", styleID );
			}
			if ( reportMean )	tableHeader ( os, "Mn " + getRatioStr ( i ) + " I", styleID );
			if ( reportStDev ) {
				tableHeader ( os, stdevMinusHTMLHeader, styleID );
				tableHeader ( os, stdevPlusHTMLHeader, styleID );
			}
			if ( reportNum )	tableHeader ( os, "Num", styleID );
		}
	}
}
void PPProteinHitQuanInfo::printHTMLEmpty ( ostream& os, const string& styleID ) const
{
	if ( reportArea ) {
		if ( reportMedian ) tableEmptyNCells ( os, numQuanPeaks, styleID );
		if ( reportIQR )	tableEmptyNCells ( os, numQuanPeaks * 2, styleID );
		if ( reportMean )	tableEmptyNCells ( os, numQuanPeaks, styleID );
		if ( reportStDev )	tableEmptyNCells ( os, numQuanPeaks * 2, styleID );
		if ( reportNum )	tableEmptyNCells ( os, numQuanPeaks, styleID );
	}
	if ( reportIntensity ) {
		if ( reportMedian )	tableEmptyNCells ( os, numQuanPeaks, styleID );
		if ( reportIQR )	tableEmptyNCells ( os, numQuanPeaks * 2, styleID );
		if ( reportMean )	tableEmptyNCells ( os, numQuanPeaks, styleID );
		if ( reportStDev )	tableEmptyNCells ( os, numQuanPeaks * 2, styleID );
		if ( reportNum )	tableEmptyNCells ( os, numQuanPeaks, styleID );
	}
}
void PPProteinHitQuanInfo::printHTML ( ostream& os, const string& styleID ) const
{
	if ( reportArea ) {
		for ( DoubleVectorSizeType i = 0 ; i < numQuanPeaks ; i++ ) {
			if ( medianAreaRatio [i] != 0.0 ) {
				if ( reportMedian ) tableCellSigFig ( os, medianAreaRatio [i], 3, false, styleID );
				if ( reportIQR ) {
					if ( lowQAreaRatio [i] != 0.0 ) {
						tableCellSigFig ( os, lowQAreaRatio [i], 3, false, styleID );
						tableCellSigFig ( os, highQAreaRatio [i], 3, false, styleID );
					}
					else tableEmptyNCells ( os, 2, styleID );
				}
				double m = meanAreaRatio [i];
				double s = stDevAreaRatio [i] * numStdDev;
				if ( reportMean )	tableCellSigFig ( os, pow ( 10.0, m ), 3, false, styleID );
				if ( reportStDev ) {
					if ( s != 0.0 ) {
						tableCellSigFig ( os, pow ( 10.0, m - s ), 3, false, styleID );
						tableCellSigFig ( os, pow ( 10.0, m + s ), 3, false, styleID );
					}
					else tableEmptyNCells ( os, 2, styleID );
				}
				if ( reportNum )	tableCell ( os, numAreaRatio [i], false, styleID );
			}
			else {
				if ( reportMedian )	tableEmptyCell ( os, styleID );
				if ( reportIQR )	tableEmptyNCells ( os, 2, styleID );
				if ( reportMean )	tableEmptyCell ( os, styleID );
				if ( reportStDev )	tableEmptyNCells ( os, 2, styleID );
				if ( reportNum )	tableEmptyCell ( os, styleID );
			}
		}
	}
	if ( reportIntensity ) {
		for ( DoubleVectorSizeType j = 0 ; j < numQuanPeaks ; j++ ) {
			if ( medianIntensityRatio [j] != 0.0 ) {
				if ( reportMedian ) tableCellSigFig ( os, medianIntensityRatio [j], 3, false, styleID );
				if ( reportIQR ) {
					if ( lowQIntensityRatio [j] != 0.0 ) {
						tableCellSigFig ( os, lowQIntensityRatio [j], 3, false, styleID );
						tableCellSigFig ( os, highQIntensityRatio [j], 3, false, styleID );
					}
					else tableEmptyNCells ( os, 2, styleID );
				}
				double m = meanIntensityRatio [j];
				double s = stDevIntensityRatio [j] * numStdDev;
				if ( reportMean )	tableCellSigFig ( os, pow ( 10.0, m ), 3, false, styleID );
				if ( reportStDev ) {
					if ( s != 0.0 ) {
						tableCellSigFig ( os, pow ( 10.0, m - s ), 3, false, styleID );
						tableCellSigFig ( os, pow ( 10.0, m + s ), 3, false, styleID );
					}
					else tableEmptyNCells ( os, 2, styleID );
				}
				if ( reportNum )	tableCell ( os, numIntensityRatio [j], false, styleID );
			}
			else {
				if ( reportMedian )	tableEmptyCell ( os, styleID );
				if ( reportIQR )	tableEmptyNCells ( os, 2, styleID );
				if ( reportMean )	tableEmptyCell ( os, styleID );
				if ( reportStDev )	tableEmptyNCells ( os, 2, styleID );
				if ( reportNum )	tableEmptyCell ( os, styleID );
			}
		}
	}
}
void PPProteinHitQuanInfo::printHeaderDelimited ( ostream& os ) const
{
	if ( reportArea ) {
		for ( DoubleVectorSizeType i = 0 ; i < numQuanPeaks ; i++ ) {
			if ( reportMedian ) delimitedHeader ( os, "Med " + getRatioStr ( i ) + " A" );
			if ( reportIQR ) {
				delimitedHeader ( os, "Q1" );
				delimitedHeader ( os, "Q3" );
			}
			if ( reportMean )	delimitedHeader ( os, "Mn " + getRatioStr ( i ) + " A" );
			if ( reportStDev ) {
				delimitedHeader ( os, stdevMinusDelimHeader );
				delimitedHeader ( os, stdevPlusDelimHeader );
			}
			if ( reportNum )	delimitedHeader ( os, "Num" );
		}
	}
	if ( reportIntensity ) {
		for ( DoubleVectorSizeType i = 0 ; i < numQuanPeaks ; i++ ) {
			if ( reportMedian ) delimitedHeader ( os, "Med " + getRatioStr ( i ) + " I" );
			if ( reportIQR ) {
				delimitedHeader ( os, "Q1" );
				delimitedHeader ( os, "Q3" );
			}
			if ( reportMean )	delimitedHeader ( os, "Mn " + getRatioStr ( i ) + " I" );
			if ( reportStDev ) {
				delimitedHeader ( os, stdevMinusDelimHeader );
				delimitedHeader ( os, stdevPlusDelimHeader );
			}
			if ( reportNum )	delimitedHeader ( os, "Num" );
		}
	}
}
void PPProteinHitQuanInfo::printDelimitedEmpty ( ostream& os ) const
{
	if ( reportArea ) {
		if ( reportMedian ) delimitedEmptyNCells ( os, numQuanPeaks );
		if ( reportIQR )	delimitedEmptyNCells ( os, numQuanPeaks * 2 );
		if ( reportMean )	delimitedEmptyNCells ( os, numQuanPeaks );
		if ( reportStDev )	delimitedEmptyNCells ( os, numQuanPeaks * 2 );
		if ( reportNum )	delimitedEmptyNCells ( os, numQuanPeaks );
	}
	if ( reportIntensity ) {
		if ( reportMedian )	delimitedEmptyNCells ( os, numQuanPeaks );
		if ( reportIQR )	delimitedEmptyNCells ( os, numQuanPeaks * 2 );
		if ( reportMean )	delimitedEmptyNCells ( os, numQuanPeaks );
		if ( reportStDev )	delimitedEmptyNCells ( os, numQuanPeaks * 2 );
		if ( reportNum )	delimitedEmptyNCells ( os, numQuanPeaks );
	}
}
void PPProteinHitQuanInfo::printDelimited ( ostream& os ) const
{
	if ( reportArea ) {
		for ( DoubleVectorSizeType i = 0 ; i < numQuanPeaks ; i++ ) {
			if ( medianAreaRatio [i] != 0.0 ) {
				if ( reportMedian ) delimitedCellSigFig ( os, medianAreaRatio [i], 3 );
				if ( reportIQR ) {
					if ( lowQAreaRatio [i] != 0.0 ) {
						delimitedCellSigFig ( os, lowQAreaRatio [i], 3 );
						delimitedCellSigFig ( os, highQAreaRatio [i], 3 );
					}
					else delimitedEmptyNCells ( os, 2 );
				}
				double m = meanAreaRatio [i];
				double s = stDevAreaRatio [i] * numStdDev;
				if ( reportMean )	delimitedCellSigFig ( os, pow ( 10.0, m ), 3 );
				if ( reportStDev ) {
					if ( s != 0.0 ) {
						delimitedCellSigFig ( os, pow ( 10.0, m - s ), 3 );
						delimitedCellSigFig ( os, pow ( 10.0, m + s ), 3 );
					}
					else delimitedEmptyNCells ( os, 2 );
				}
				if ( reportNum )	delimitedCell ( os, numAreaRatio [i] );
			}
			else {
				if ( reportMedian )	delimitedEmptyCell ( os );
				if ( reportIQR )	delimitedEmptyNCells ( os, 2 );
				if ( reportMean )	delimitedEmptyCell ( os );
				if ( reportStDev )	delimitedEmptyNCells ( os, 2 );
				if ( reportNum )	delimitedEmptyCell ( os );
			}
		}
	}
	if ( reportIntensity ) {
		for ( DoubleVectorSizeType j = 0 ; j < numQuanPeaks ; j++ ) {
			if ( medianIntensityRatio [j] != 0.0 ) {
				if ( reportMedian ) delimitedCellSigFig ( os, medianIntensityRatio [j], 3 );
				if ( reportIQR ) {
					if ( lowQIntensityRatio [j] != 0.0 ) {
						delimitedCellSigFig ( os, lowQIntensityRatio [j], 3 );
						delimitedCellSigFig ( os, highQIntensityRatio [j], 3 );
					}
					else delimitedEmptyNCells ( os, 2 );
				}
				double m = meanIntensityRatio [j];
				double s = stDevIntensityRatio [j] * numStdDev;
				if ( reportMean )	delimitedCellSigFig ( os, pow ( 10.0, m ), 3 );
				if ( reportStDev ) {
					if ( s != 0.0 ) {
						delimitedCellSigFig ( os, pow ( 10.0, m - s ), 3 );
						delimitedCellSigFig ( os, pow ( 10.0, m + s ), 3 );
					}
					else delimitedEmptyNCells ( os, 2 );
				}
				if ( reportNum )	delimitedCell ( os, numIntensityRatio [j] );
			}
			else {
				if ( reportMedian )	delimitedEmptyCell ( os );
				if ( reportIQR )	delimitedEmptyNCells ( os, 2 );
				if ( reportMean )	delimitedEmptyCell ( os );
				if ( reportStDev )	delimitedEmptyNCells ( os, 2 );
				if ( reportNum )	delimitedEmptyCell ( os );
			}
		}
	}
}
void PPProteinHitQuanInfo::setQuanParams ( const ParameterList* params )
{
	reportArea		= params->getBoolValue ( "rep_a_lh_area" );
	reportIntensity	= params->getBoolValue ( "rep_a_lh_int" );
	reportMedian	= params->getBoolValue ( "rep_q_median" );
	reportIQR		= params->getBoolValue ( "rep_q_iqr" );
	reportMean		= params->getBoolValue ( "rep_q_mean" );
	reportStDev		= params->getBoolValue ( "rep_q_stdev" );
	reportNum		= params->getBoolValue ( "rep_q_num" );
	numStdDev		= params->getDoubleValue ( "rep_q_n_sdv", 2.0 );
	stdevMinusHTMLHeader	= "-" + gen_ftoa ( numStdDev, "%.1f" ) + "&#963;";
	stdevPlusHTMLHeader		= "+" + gen_ftoa ( numStdDev, "%.1f" ) + "&#963;";
	stdevMinusDelimHeader	= "-" + gen_ftoa ( numStdDev, "%.1f" ) + "sd";
	stdevPlusDelimHeader	= "+" + gen_ftoa ( numStdDev, "%.1f" ) + "sd";
}

PPProteinHitInfo::PPProteinHitInfo () :
	SearchResultsProteinInfo ( 0.0, 0.0, 0, 0 ),
	bestDiscriminantScore ( DiscriminantScore::MIN_DISC_SCORE ),
	totalDiscriminantScore ( DiscriminantScore::MIN_DISC_SCORE )
{
	ppphqi.initQuan ();
}
PPProteinHitInfo::PPProteinHitInfo ( const string& accessionNumber ) :
	SearchResultsProteinInfo ( 0.0, 0.0, 0, 0 ),
	bestDiscriminantScore ( DiscriminantScore::MIN_DISC_SCORE ),
	accessionNumber ( accessionNumber ),
	totalDiscriminantScore ( DiscriminantScore::MIN_DISC_SCORE )
{
	ppphqi.initQuan ();
}
PPProteinHitInfo::PPProteinHitInfo ( const StringVector& str ) :
	SearchResultsProteinInfo ( 0.0, 0.0, 0, 0 ),
	bestDiscriminantScore ( DiscriminantScore::MIN_DISC_SCORE ),
	totalDiscriminantScore ( DiscriminantScore::MIN_DISC_SCORE )
{
	if ( !str.empty () ) {
		accessionNumber = XMLParser::getStringValue	( str [0], "accession_number" );
	}
	ppphqi.initQuan ();
}
int PPProteinHitInfo::getColspan () const
{
	int colspan = SearchResultsProteinInfo::getColspan ();
	if ( reportCoverage )		colspan++;
	if ( reportBestDiscScore )	colspan++;
	if ( reportBestExpectVal )	colspan++;
	colspan += ppphqi.getColspan ();
	return colspan;
}
string PPProteinHitInfo::getRatioStr ( int i ) const
{
	return ppphqi.getRatioStr ( i );
}
void PPProteinHitInfo::printHeaderHTML ( ostream& os, int index, const string& styleID ) const
{
	SearchResultsProteinInfo::printHeaderHTML ( os, index, styleID );
	if ( reportCoverage )		tableHeader ( os, "% Cov", styleID );
	if ( reportBestDiscScore )	tableHeader ( os, "Best Disc Score", styleID );
	if ( reportBestExpectVal )	tableHeader ( os, "Best Expect Val", styleID );
	ppphqi.printHeaderHTML ( os, styleID );
}
void PPProteinHitInfo::printHTML ( ostream& os, bool empty, const string& styleID ) const
{
	SearchResultsProteinInfo::printHTML ( os, empty, styleID );
	if ( score == 0.0 || empty ) {
		if ( reportCoverage )		tableEmptyCell ( os, styleID );
		if ( reportBestDiscScore )	tableEmptyCell ( os, styleID );
		if ( reportBestExpectVal )	tableEmptyCell ( os, styleID );
		ppphqi.printHTMLEmpty ( os, styleID );
	}
	else {
		if ( reportCoverage ) tableCell ( os, coverageMap.getPercentCoverage (), 1, false, styleID );
		if ( reportBestDiscScore ) {
			if ( bestDiscriminantScore == DiscriminantScore::MIN_DISC_SCORE ) {
				tableCell ( os, "---", false, false, styleID );
			}
			else {
				tableCell ( os, bestDiscriminantScore, 2, false, styleID );
			}
		}
		if ( reportBestExpectVal ) tableCellSigFig ( os, bestExpectationValue, 2, false, styleID );
		ppphqi.printHTML ( os, styleID );
	}
}
void PPProteinHitInfo::printHeaderDelimited ( ostream& os, int index ) const
{
	SearchResultsProteinInfo::printHeaderDelimited ( os, index );
	if ( reportCoverage )		delimitedHeader ( os, "% Cov" );
	if ( reportBestDiscScore )	delimitedHeader ( os, "Best Disc Score" );
	if ( reportBestExpectVal )	delimitedHeader ( os, "Best Expect Val" );
	ppphqi.printHeaderDelimited ( os );
}
void PPProteinHitInfo::printDelimited ( ostream& os ) const
{
	SearchResultsProteinInfo::printDelimited ( os );
	if ( score == 0.0 ) {
		if ( reportCoverage )		delimitedEmptyCell ( os );
		if ( reportBestDiscScore )	delimitedEmptyCell ( os );
		if ( reportBestExpectVal )	delimitedEmptyCell ( os );
		ppphqi.printDelimitedEmpty ( os );
	}
	else {
		if ( reportCoverage ) delimitedCell ( os, coverageMap.getPercentCoverage (), 1 );
		if ( reportBestDiscScore ) {
			if ( bestDiscriminantScore == DiscriminantScore::MIN_DISC_SCORE )
				delimitedCell ( os, "---" );
			else
				delimitedCell ( os, bestDiscriminantScore, 2 );
		}
		if ( reportBestExpectVal ) delimitedCellSigFig ( os, bestExpectationValue, 2 );
		ppphqi.printDelimited ( os );
	}
}
void PPProteinHitInfo::printDelimitedEmpty ( ostream& os ) const
{
	SearchResultsProteinInfo::printDelimitedEmpty ( os );
	if ( reportCoverage )		delimitedEmptyCell ( os );
	if ( reportBestDiscScore )	delimitedEmptyCell ( os );
	if ( reportBestExpectVal )	delimitedEmptyCell ( os );
	ppphqi.printDelimitedEmpty ( os );
}
PPXMLProteinHit::PPXMLProteinHit ( const string& aNum, vector <SearchResultsPeptideHit*>& th ) :
	SearchResultsProteinHit ( new PPProteinHitInfo ( aNum ), aNum ),
	tagHits ( th )
{
	double protBestPepScore = 0.0;

	for ( StringVectorSizeType i = 0 ; i < th.size () ; i++ ) {
		protBestPepScore = genMax ( th [i]->getScore (), protBestPepScore );
	}
	getPPProteinHitInfo ()->setBestScore ( protBestPepScore );
}
class CheckTagHitsPValue {
	double minPeptideScore;
	double minBestDiscScore;
	double maxPeptideEValue;
	MapPairStringIntToDouble bestPeptideScore;	// map <pair <hitSequence, charge>, discriminantScore>
public:
	CheckTagHitsPValue ( double minPeptideScore, double minBestDiscScore, double maxPeptideEValue, const SearchResultsPeptideHitPtrVector& tagHits ) :
		minPeptideScore ( minPeptideScore ),
		minBestDiscScore ( minBestDiscScore ),
		maxPeptideEValue ( maxPeptideEValue )
	{
		if ( !sresKeepReplicates ) {
			for ( SearchResultsPeptideHitPtrVectorConstIterator i = tagHits.begin () ; i != tagHits.end () ; i++ ) {
				string peptide = (*i)->getHitSequence ();
				int charge = sresKeepCharges ? (*i)->getCharge () : 0;
				PairStringInt psi = make_pair ( peptide, charge );
				double pValue = (*i)->getPValue ();
				MapPairStringIntToDoubleConstIterator cur = bestPeptideScore.find ( psi );
				if ( cur == bestPeptideScore.end () || ( pValue != 1.0 && pValue < (*cur).second ) ) {
					bestPeptideScore [psi] = pValue;	// map <pair <hitSequence, charge>, pValue>
				}
			}
		}
	}
	bool operator () ( const SearchResultsPeptideHit* th )
	{
		if ( th->getExpectationValue () > maxPeptideEValue ) return true;
		if ( th->getScore () < minPeptideScore ) return true;
		if ( th->getDiscriminantScore () < minBestDiscScore ) return true;
		if ( !sresKeepReplicates ) {
			PairStringInt psi = make_pair ( th->getHitSequence (), sresKeepCharges ? th->getCharge () : 0 );
			MapPairStringIntToDoubleIterator cur = bestPeptideScore.find ( psi );
			if ( th->getPValue () > (*cur).second ) {
				return true;
			}
			else {	// If the peptide is retained increment the best score for the peptide slightly so that a second peptide
					// with the same score isn't retained.
				(*cur).second *= 0.9999999;
			}
		}
		return false;
	}
};
class CheckTagHitsDiscScore {
	double minPeptideScore;
	double minBestDiscScore;
	double maxPeptideEValue;
	MapPairStringIntToDouble bestPeptideScore;	// map <pair <hitSequence, charge>, discriminantScore>
public:
	CheckTagHitsDiscScore ( double minPeptideScore, double minBestDiscScore, double maxPeptideEValue, const SearchResultsPeptideHitPtrVector& tagHits ) :
		minPeptideScore ( minPeptideScore ),
		minBestDiscScore ( minBestDiscScore ),
		maxPeptideEValue ( maxPeptideEValue )
	{
		if ( !sresKeepReplicates ) {
			for ( SearchResultsPeptideHitPtrVectorConstIterator i = tagHits.begin () ; i != tagHits.end () ; i++ ) {
				string peptide = (*i)->getHitSequence ();
				int charge = sresKeepCharges ? (*i)->getCharge () : 0;
				PairStringInt psi = make_pair ( peptide, charge );
				double discScore = (*i)->getDiscriminantScore ();
				MapPairStringIntToDoubleConstIterator cur = bestPeptideScore.find ( psi );
				if ( cur == bestPeptideScore.end () || discScore > (*cur).second ) {
					bestPeptideScore [psi] = discScore;	// map <pair <hitSequence, charge>, discriminantScore>
				}
			}
		}
	}
	bool operator () ( const SearchResultsPeptideHit* th )
	{
		if ( th->getExpectationValue () > maxPeptideEValue ) return true;
		if ( th->getScore () < minPeptideScore ) return true;
		if ( th->getDiscriminantScore () < minBestDiscScore ) return true;
		if ( !sresKeepReplicates ) {
			PairStringInt psi = make_pair ( th->getHitSequence (), sresKeepCharges ? th->getCharge () : 0 );
			MapPairStringIntToDoubleIterator cur = bestPeptideScore.find ( psi );
			if ( th->getDiscriminantScore () < (*cur).second ) {
				return true;
			}
			else {	// If the peptide is retained increment the best score for the peptide slightly so that a second peptide
					// with the same score isn't retained.
				(*cur).second += 0.000001;
			}
		}
		return false;
	}
};
class CheckTagHits2 {
	const MapSpecIDBestDiscriminantScore& bestScores;
	SetString setSpecID;
public:
	CheckTagHits2 ( const MapSpecIDBestDiscriminantScore& bestScores ) :
		bestScores ( bestScores )
	{}
	bool operator () ( const SearchResultsPeptideHit* th )
	{
		const PeptidePosition* pp = th->getPeptidePosition ();
		MapSpecIDBestDiscriminantScore::const_iterator cur2 = bestScores.find ( pp->getSpecID () );
		bool mm = matrixMatch ( cur2->second.second, pp->getDBPeptide () );
		if ( mm && th->getScore () < cur2->second.first.second ) return true;
		if ( !mm && th->getDiscriminantScore () < cur2->second.first.first ) return true;
		pair <SetStringIterator, bool> flag = setSpecID.insert ( pp->getSpecID () );
		if ( flag.second == false ) return true;	// Only one entry allowed per specID
		return false;
	}
};

class CheckTagHitsIntersection {
	string id;
	int numMatches;
public:
	CheckTagHitsIntersection ( const string& id, int numMatches ) :
		id ( id ),
		numMatches ( numMatches )
	{
	}
	bool operator () ( SearchResultsPeptideHit* th )
	{
		int matches = 0;
		for ( int i = 0 ; i < SearchResults::numSearches ; i++ ) {
			if ( SearchResults::vsr [i]->getPPPeptideHitInfoPair ( th, id ) != 0 ) matches++;
		}
		return matches != numMatches;
	}
};
void PPXMLProteinHit::calculateDiscriminantScores ( MapSpecIDAndPeptideDiscriminantScore& dsMap, MapSpecIDBestDiscriminantScore& bestScores, bool setRepeats, const DiscriminantScore& discScore, double maxPeptideEValue )
{
	double protBestPepScore = proteinHitInfo->getBestScore ();
	double bestDiscScore = DiscriminantScore::MIN_DISC_SCORE;
	SetIntPtr repeatsPtr;
	for ( SearchResultsPeptideHitPtrVectorSizeType i = 0 ; i < tagHits.size () ; i++ ) {
		SearchResultsPeptideHit* th = tagHits [i];
		double discriminantScore = discScore.calculateDiscriminantScore ( protBestPepScore, th->getScoreDifference (), th->getExpectationValue () );

		if ( th->getPeptidePosition ()->getMmod () ) discriminantScore -= 0.000001;	// Reduce discriminant score for mass modified peptides
		string specID = th->getPeptidePosition ()->getSpecID ();
		const string& peptide = th->getPeptidePosition ()->getDBPeptide ();
		PairStringString sidpep;
		sidpep.first = specID;
		sidpep.second = peptide;
		string& newPeptide = sidpep.second;
		for ( StringSizeType j = 0 ; j < peptide.length () ; j++ ) {
			if ( peptide [j] == 'I' ) newPeptide [j] = 'L';
			if ( peptide [j] == 'Q' ) newPeptide [j] = 'K';
		}
		MapSpecIDAndPeptideDiscriminantScoreIterator cur1 = dsMap.find ( sidpep );
		if ( cur1 == dsMap.end () ) {
			pair <double, int> info;
			info.first = discriminantScore;
			if ( setRepeats ) {
				info.second = 0;
			}
			dsMap [sidpep] = info;
		}
		else {
			(*cur1).second.first = discriminantScore;
		}
		if ( setRepeats ) repeatsPtr.insert ( &(dsMap [sidpep].second) );
		MapSpecIDBestDiscriminantScoreConstIterator cur = bestScores.find ( specID );
		if ( ( cur == bestScores.end () || (*cur).second.first.first < discriminantScore ) && th->getExpectationValue () <= maxPeptideEValue ) {
			pair <pair<double,double>, string> pds;
			pds.first.first = discriminantScore;
			pds.first.second = th->getScore ();
			pds.second = peptide;
			bestScores [specID] = pds;
		}
		th->setDiscriminantScore ( discriminantScore );
		if ( !setRepeats ) th->setRepeats ( dsMap [sidpep].second );
		bestDiscScore = genMax ( discriminantScore, bestDiscScore );
	}
	if ( setRepeats ) {
		for ( SetIntPtrIterator i = repeatsPtr.begin () ; i != repeatsPtr.end () ; i++ ) {
			(**i)++;		// Increments the peptide repeats counter once per protein
		}
	}
	getPPProteinHitInfo ()->setBestDiscriminantScore ( bestDiscScore );
}
void PPXMLProteinHit::filterPeptides ( bool bestDiscrimOnly, bool noExpectation, double minPeptideScore, double minBestDiscScore, double maxPeptideEValue, MapSpecIDBestDiscriminantScore& bestScores )
{
	if ( bestDiscrimOnly ) {	// This removes all but the best discriminant score for each spectrum
		int i = remove_if ( tagHits.begin (), tagHits.end (), CheckTagHits2 ( bestScores ) ) - tagHits.begin ();
		tagHits.erase ( tagHits.begin () + i, tagHits.end () );
	}
	if ( noExpectation ) {			// Select the best peptide based on discriminant scores
		int i = remove_if ( tagHits.begin (), tagHits.end (), CheckTagHitsDiscScore ( minPeptideScore, minBestDiscScore, maxPeptideEValue, tagHits ) ) - tagHits.begin ();
		tagHits.erase ( tagHits.begin () + i, tagHits.end () );
	}
	else {	// Select the best peptide based on p values
		int i = remove_if ( tagHits.begin (), tagHits.end (), CheckTagHitsPValue ( minPeptideScore, minBestDiscScore, maxPeptideEValue, tagHits ) ) - tagHits.begin ();
		tagHits.erase ( tagHits.begin () + i, tagHits.end () );
	}
}
void PPXMLProteinHit::filterPeptidesIntersection ( const string& id, int numMatches )
{
	int i = remove_if ( tagHits.begin (), tagHits.end (), CheckTagHitsIntersection ( id, numMatches ) ) - tagHits.begin ();
	tagHits.erase ( tagHits.begin () + i, tagHits.end () );
}
class RemoveTagHits {
	VectorPairStringString removePeps;
public:
	RemoveTagHits ( const VectorPairStringString& removePeps ) :
		removePeps ( removePeps ) {}
	bool operator () ( const SearchResultsPeptideHit* th )
	{
		const string& specID = th->getPeptidePosition ()->getSpecID ();
		const string& pep = th->getPeptidePosition ()->getPeptide ();
		for ( int i = 0 ; i < removePeps.size () ; i++ ) {
			if ( specID == removePeps [i].first && pep == removePeps [i].second ) return true;
		}
		return false;
	}
};
class SortTagHitAscending {
	public:
		bool operator () ( const SearchResultsPeptideHit* lhs, const SearchResultsPeptideHit* rhs ) const
		{
			return lhs->getPeptidePosition ()->getSequence () < rhs->getPeptidePosition ()->getSequence ();
		}
};
void PPXMLProteinHit::removePeptides ( const VectorPairStringString& removePeps )
{
	int i = remove_if ( tagHits.begin (), tagHits.end (), RemoveTagHits ( removePeps ) ) - tagHits.begin ();
	tagHits.erase ( tagHits.begin () + i, tagHits.end () );
}
void PPXMLProteinHit::calculateStats ()
{
	double protBestPepScore = 0.0;
	double protScore = 0.0;
	double bestPepScore = 0.0;
	double bestDScore = 0.0;
	double bestEScore = 10E+20;
	double totalDiscScore = 0.0;
	double bestExpectationValue = 10E+20;
	double bestDiscScore = DiscriminantScore::MIN_DISC_SCORE;
	const PeptidePosition* oldPP;
	int numUnique = 0;
	stable_sort ( tagHits.begin (), tagHits.end (), SortTagHitAscending () );
	for ( SearchResultsPeptideHitPtrVectorSizeType j = 0 ; j < tagHits.size () ; j++ ) {
		SearchResultsPeptideHit* th = tagHits [j];
		double pepScore = th->getScore ();
		double dScore = th->getDiscriminantScore ();
		double eVal = th->getExpectationValue ();
		protBestPepScore = genMax ( pepScore, protBestPepScore );
		bestDiscScore = genMax ( dScore, bestDiscScore );
		bestExpectationValue = genMin ( eVal, bestExpectationValue );
		const PeptidePosition* pp = th->getPeptidePosition ();
		if ( j != 0 && ( pp->getSequence () == oldPP->getSequence () || pp->getSpecID () == oldPP->getSpecID () ) ) {		// Same as previous sequence
			if ( eVal < bestEScore ) {
				protScore -= bestPepScore;
				totalDiscScore -= bestDScore;
				bestPepScore = pepScore;
				bestDScore = dScore;
				bestEScore = eVal;
				protScore += pepScore;
				totalDiscScore += dScore;
			}
		}
		else {
			oldPP = pp;
			bestPepScore = pepScore;
			bestDScore = dScore;
			bestEScore = eVal;
			protScore += pepScore;
			totalDiscScore += dScore;
			numUnique++;
		}
	}
	getPPProteinHitInfo ()->setScore ( protScore );
	getPPProteinHitInfo ()->setBestScore ( protBestPepScore );
	getPPProteinHitInfo ()->setBestDiscriminantScore ( bestDiscScore );
	getPPProteinHitInfo ()->setTotalDiscriminantScore ( totalDiscScore );
	getPPProteinHitInfo ()->setBestExpectationValue ( bestExpectationValue );
	getPPProteinHitInfo ()->setNumUnique ( numUnique );
	getPPProteinHitInfo ()->setPeptideCount ( tagHits.size () );
}
void PPXMLProteinHit::setProteinStats ()
{
	string accNum = accessionNumber;
	ProteinInfo proteinInfo ( accNum );
	CoverageMap coverageMap ( proteinInfo.getLength () );
	const PeptidePosition* oldPP;
	bool errorReported = false;
	for ( SearchResultsPeptideHitPtrVectorSizeType i = 0 ; i < tagHits.size () ; i++ ) {
		const PeptidePosition* pp = tagHits [i]->getPeptidePosition ();
		bool flag = i != 0 && ( pp->getSequence () == oldPP->getSequence () || pp->getSpecID () == oldPP->getSpecID () );
		if ( !flag ) {		// Same as previous sequence
			try {
				coverageMap.setCoverage ( pp->getStartAA (), pp->getEndAA () );
			}
			catch ( runtime_error e ) {
				if ( !errorReported ) {
					ErrorHandler::genError ()->message ( "Accession Number: " + accNum + ". " + string ( e.what () ) + "\n" );
					errorReported = true;
				}
			}
			oldPP = pp;
		}
	}
	getPPProteinHitInfo ()->setCoverageMap ( coverageMap );
}
class CheckProteinHits {
	double minProteinScore;
	double minBestDiscScore;
	double maxEValue;
public:
	CheckProteinHits ( double minProteinScore, double minBestDiscScore, double maxEValue ) :
		minProteinScore ( minProteinScore ),
		minBestDiscScore ( minBestDiscScore ),
		maxEValue ( maxEValue )
	{
	}
	bool operator () ( const SearchResultsProteinHit* ph )
	{
		if ( ph->getNumMSMSHits () == 0 ) return true;
		if ( ph->getScore () < minProteinScore ) return true;
		if ( ph->getBestDiscriminantScore () < minBestDiscScore ) return true;
		if ( ph->getBestExpectationValue () > maxEValue ) return true;
		return false;
	}
};
DatabaseResults::DatabaseResults ( const StringVector& database, istream& istr ) :
	database ( database )
{
	XMLIStreamList xstr ( istr, "pre_search_results" );
	for ( int i = 0 ; i < database.size () ; i++ ) {
		string preSearch;
		xstr.getNext ( preSearch );
		numDatabaseEntries.push_back ( XMLParser::getIntValue ( preSearch, "num_entries" ) );
		numSearchedEntries.push_back ( XMLParser::getIntValue ( preSearch, "num_final_indicies" ) );
	}
}
void DatabaseResults::printHTML ( ostream& os ) const
{
	for ( int i = 0 ; i < database.size () ; i++ ) {
		os << "Database: <b>" << database [i] << "</b>" << " ";
		if ( numDatabaseEntries [i] )
			os << "(" << numSearchedEntries [i] << "/" << numDatabaseEntries [i] << " entries searched)";
		else
			os << "No database entries searched";

		os << "<br />";
		os << endl;
	}
}
bool DatabaseResults::isConcat () const
{
	for ( StringVectorSizeType i = 0 ; i < database.size () ; i++ ) {
		if ( database [i].find ( "concat" ) != string::npos ) return true;
	}
	return false;
}
VectorSearchResultsPtr SearchResults::vsr;
int SearchResults::numSearches = 0;
bool SearchResults::noMergeExpectation = false;
MapIDAccNumSearchResultsProteinHitPtrConstIterator SearchResults::compProtHits;
MapIDMapAccNumSearchResultsPeptideHitPtrVector SearchResults::compPepHits;

const string SearchResults::defaultID = "id";

SearchResults::SearchResults ( const string& projectName, const string& resultsName, const string& fname, const SearchCompareParams& params, int fileIndex, const string& searchEndTime, const string& searchTime ) :
	projectName ( projectName ),
	resultsName ( resultsName ),
	fname ( fname ),
	discScoreGraph ( params.getDiscScoreGraph () ),
	modificationScoreThreshold ( params.getModificationScoreThreshold () ),
	discFilename ( fname.substr ( 0, fname.length () - 4 ) + ".disc.txt" ),
	discFilenameExists ( genFileExists ( discFilename ) ),
	searchEndTime ( searchEndTime ),
	searchTime ( searchTime )
{
	numSearches++;
	vsr.push_back ( this );
	double minPeptideScore = params.getMinPPPeptideScore ();
	double maxPeptideEValue = params.getMaxPeptideEValue ();
	bool multisample = params.getMultiSample ();
	SetInt idFilterSet = params.getIDFilterSet ();
	MapIDMapAccNoAndVectorSearchResultsPeptideHit mimab;
	string instrument = getMapTagHits ( mimab, fname, multisample, idFilterSet, minPeptideScore, maxPeptideEValue, params.getXLMinLowScore (), params.getXLMinScoreDiff (), params.getXLMaxLowExpectation () );
	double minBestDiscScore = params.getMinBestDiscScore ( instrument );
	DiscriminantScore discScore ( instrument );
	MapSpecIDBestDiscriminantScore bestScores;
	MapSpecIDAndPeptideDiscriminantScore dsMap;
	ujm->writeMessage ( cout, "Creating protein list" );
	int num = 0;
	for ( MapIDMapAccNoAndVectorSearchResultsPeptideHit::iterator i = mimab.begin () ; i != mimab.end () ; i++ ) {	// Iterate through the IDs
		MapAccNoAndVectorSearchResultsPeptideHit& mapTagHits = (*i).second;
		for ( MapAccNoAndVectorSearchResultsPeptideHit::iterator j = mapTagHits.begin () ; j != mapTagHits.end () ; j++ ) {	// Iterate through the Accession Numbers
			num++;
			if ( num % 10000 == 0 ) ujm->writeMessage ( cout, "Creating protein list, " + gen_itoa ( num ) + " proteins" );
			PPXMLProteinHit* ppxph = new PPXMLProteinHit ( (*j).first, (*j).second );
			proteinHits [(*i).first].push_back ( ppxph );
			ppxph->calculateDiscriminantScores ( dsMap, bestScores, true, discScore, maxPeptideEValue );
		}
		processHits ( proteinHits [(*i).first], noExpectation, minBestDiscScore, bestScores, dsMap, discScore, params, fileIndex, (*i).first );
	}
	for ( MapIDSearchResultsProteinHitPtrVector::iterator k = proteinHits.begin () ; k != proteinHits.end () ; k++ ) {
		string id = (*k).first;
		setPeptideHitsAndMap ( id );
		setProteinMap ( id );
	}
	if ( discScoreGraph ) createHistogram ( dsMap );
	emptyProteinHit = new PPProteinHitInfo;
}

void SearchResults::mergeResults ( const SearchCompareParams& params )
{
	string reportHitsType = params.getReportHitsType ();
	for ( SetStringConstIterator ii = vsr [0]->idSet.begin () ; ii != vsr [0]->idSet.end () ; ii++ ) {
		string id = *ii;
		SearchResultsPeptideHitPtrVector reportedHits;
		copy ( vsr [0]->peptideMap [id].begin (), vsr [0]->peptideMap [id].end (), back_inserter ( reportedHits ) );
		for ( int i = 1 ; i < numSearches ; i++ ) {
			SetSearchResultsPeptideHit& pepPosn = vsr [i]->peptideMap [id];
			SearchResultsPeptideHitPtrVector tempHits;
			merge ( reportedHits.begin (), reportedHits.end (), pepPosn.begin (), pepPosn.end (), back_inserter ( tempHits ), SortPeptidePositionAscending2 () );
			reportedHits = tempHits;
		}
		string prevAccNum;
		SearchResultsPeptideHitPtrVector accHits;
		MapSpecIDBestDiscriminantScore bestScores;
		SearchResultsPeptideHitPtrVector rHits;
		for ( SearchResultsPeptideHitPtrVectorSizeType kk = 0 ; kk < reportedHits.size () ; kk++ ) {
			string aNum = reportedHits [kk]->getAccessionNumber ();
			if ( kk == 0 ) prevAccNum = aNum;
			if ( aNum != prevAccNum ) {
				compPepHits [id][prevAccNum] = accHits;
				prevAccNum = aNum;
				accHits.clear ();
			}
			accHits.push_back ( reportedHits [kk] );
			if ( kk == reportedHits.size () - 1 ) {
				compPepHits [id][aNum] = accHits;
			}
		}
		for ( MapAccNumSearchResultsPeptideHitPtrVectorIterator mm = compPepHits [id].begin () ; mm != compPepHits [id].end () ; mm++ ) {
			compProtHits [id][(*mm).first] = new PPXMLProteinHit ( (*mm).first, (*mm).second );
		}
		{
			double minProteinScore = params.getMinPPProteinScore ();
			double maxProteinEValue = params.getMaxProteinEValue ();
			double maxPeptideEValue = params.getMaxPeptideEValue ();
			double minPeptideScore = params.getMinPPPeptideScore ();
			bool bestDiscrimOnly = params.getBestDiscOnly ();
			for ( MapAccNumSearchResultsProteinHitPtrConstIterator i = compProtHits [id].begin () ; i != compProtHits [id].end () ; i++ ) {
				static_cast <PPXMLProteinHit*> ((*i).second)->filterPeptides ( false, noMergeExpectation, minPeptideScore, -100, maxPeptideEValue, bestScores );
			}
			if ( reportHitsType == "Peptide Intersection" || reportHitsType == "Peptide Difference" ) {
				int numMatches = ( reportHitsType == "Peptide Intersection" ) ? numSearches : 1;
				for ( MapAccNumSearchResultsProteinHitPtrConstIterator i = compProtHits [id].begin () ; i != compProtHits [id].end () ; i++ ) {
					static_cast <PPXMLProteinHit*> ((*i).second)->filterPeptidesIntersection ( id, numMatches );
				}
			}
			for ( MapAccNumSearchResultsProteinHitPtrConstIterator j = compProtHits [id].begin () ; j != compProtHits [id].end () ; j++ ) {
				static_cast <PPXMLProteinHit*> ((*j).second)->calculateStats ();
			}
			if ( params.getReportCoverage () ) {
				for ( MapAccNumSearchResultsProteinHitPtrConstIterator b = compProtHits [id].begin () ; b != compProtHits [id].end () ; b++ ) {
					static_cast <PPXMLProteinHit*> ((*b).second)->setProteinStats ();
				}
			}
		}
	}
}
void SearchResults::processHits ( SearchResultsProteinHitPtrVector& pHits, bool noExpectation, double minBestDiscScore, MapSpecIDBestDiscriminantScore& bestScores, MapSpecIDAndPeptideDiscriminantScore& dsMap, DiscriminantScore& discScore, const SearchCompareParams& params, int fileIndex, const string& id )
{
	double minProteinScore = params.getMinPPProteinScore ();
	double maxProteinEValue = params.getMaxProteinEValue ();
	double maxPeptideEValue = params.getMaxPeptideEValue ();
	double minPeptideScore = params.getMinPPPeptideScore ();
	bool bestDiscrimOnly = params.getBestDiscOnly ();
	for ( PPXMLProteinHitPtrVectorSizeType i = 0 ; i < pHits.size () ; i++ ) {
		static_cast <PPXMLProteinHit*> (pHits[i])->filterPeptides ( bestDiscrimOnly, noExpectation, minPeptideScore, minBestDiscScore, maxPeptideEValue, bestScores );
	}
	for ( PPXMLProteinHitPtrVectorSizeType j = 0 ; j < pHits.size () ; j++ ) {
		static_cast <PPXMLProteinHit*> (pHits[j])->calculateStats ();
	}
	int k = remove_if ( pHits.begin (), pHits.end (), CheckProteinHits ( minProteinScore, minBestDiscScore, maxProteinEValue ) ) - pHits.begin ();
	// Delete pointers
	pHits.erase ( pHits.begin () + k, pHits.end () );
	bestScores.clear ();
	for ( PPXMLProteinHitPtrVectorSizeType x = 0 ; x < pHits.size () ; x++ ) {
		static_cast <PPXMLProteinHit*> (pHits[x])->calculateDiscriminantScores ( dsMap, bestScores, false, discScore, maxPeptideEValue );
	}
	for ( PPXMLProteinHitPtrVectorSizeType y = 0 ; y < pHits.size () ; y++ ) {
		static_cast <PPXMLProteinHit*> (pHits[y])->filterPeptides ( bestDiscrimOnly, noExpectation, minPeptideScore, minBestDiscScore, maxPeptideEValue, bestScores );
	}
	if ( params.getRemovePeptides () ) {
		VectorPairStringString removePeptides = params.getRemovePeptides ( fileIndex, id );
		for ( PPXMLProteinHitPtrVectorSizeType yy = 0 ; yy < pHits.size () ; yy++ ) {
			if ( removePeptides.size () ) static_cast <PPXMLProteinHit*> (pHits[yy])->removePeptides ( removePeptides );
		}
	}
	for ( PPXMLProteinHitPtrVectorSizeType z = 0 ; z < pHits.size () ; z++ ) {
		static_cast <PPXMLProteinHit*> (pHits[z])->calculateStats ();
	}
	int a = remove_if ( pHits.begin (), pHits.end (), CheckProteinHits ( minProteinScore, minBestDiscScore, maxProteinEValue ) ) - pHits.begin ();
	// Delete pointers
	pHits.erase ( pHits.begin () + a, pHits.end () );
	if ( params.getReportCoverage () ) {
		for ( PPXMLProteinHitPtrVectorSizeType b = 0 ; b < pHits.size () ; b++ ) {
			static_cast <PPXMLProteinHit*> (pHits[b])->setProteinStats ();
		}
	}
}
void SearchResults::createHistogram ( const MapSpecIDAndPeptideDiscriminantScore& dsMap )
{
	if ( discFilenameExists ) {	// Use cache file if it there
		if ( genFileSize ( discFilename ) != 0 ) {
			GenIFStream ifs ( discFilename );
			ifs >> histogram;
		}
	}
	else {
		for ( MapSpecIDAndPeptideDiscriminantScoreConstIterator i = dsMap.begin () ; i != dsMap.end () ; i++ ) {
			histogram.add ( (*i).second.first );
		}
	}
}

string SearchResults::getMapTagHits ( MapIDMapAccNoAndVectorSearchResultsPeptideHit& mapTagHits, const string& fname, bool multisample, const SetInt& idFilterSet, double minPeptideScore, double maxPeptideEValue, double xlMinLowScore, double xlMinScoreDiff, double xlMaxLowExpectation )
{
	MapIDExpectationValueInfo mievi;
	bool eValueFlag;
	bool linearTailFitExpectation;
	string instrument;
	string expName;
	IntVector iv;
	bool idFilter = true;
	for ( int i = 0 ; ; i++ ) {
		string curFile;
		if ( genFileExists ( fname ) ) curFile = fname;
		else {
			string f = fname + string ( "_" ) + gen_itoa ( i );
			if ( genFileExists ( f ) ) curFile = f;
			else {
				if ( i == 0 )
					ErrorHandler::genError ()->error ( "Results file not found.\n" );
				else
					break;
			}
		}
		GenIFStream ist ( curFile );
		if ( i == 0 ) {
			XMLIStreamList xstr ( ist, "parameters" );
			string pStr;
			xstr.getNext ( pStr );
			ParameterList* pList = new ParameterList ( pStr, false, false, false ); 
			iv = ProteinInfo::initialise ( pList );
			string outputDir = getProjectDir ( pList );
			expName = outputDir + SLASH + pList->getStringValue  ( "expect_coeff_file" ) + ".xml";
			eValueFlag = getEValues ( expName, mievi );
			PeptidePosition::initialiseParams ( pList );
			ExtraUserMods::instance ().addUserMods2 ( pList, ExtraUserModsForm::getNumUserMods () );
			SCModInfo::init ( numSearches, pList->getStringVectorValue ( "const_mod" ), modificationScoreThreshold );
			instrument = pList->getStringValue ( "instrument_name" );
			string expectationCalculationType = pList->getStringValue ( "expect_calc_method", "None" );
			noExpectation = ( expectationCalculationType == "None" );
			noMergeExpectation |= noExpectation;
			linearTailFitExpectation = ( expectationCalculationType == "Linear Tail Fit" );
			PeptideSpectralInfo::setExpectationFlag ( expectationCalculationType );
			spottingPlate = ProjectFile::getSpottingPlate ( pList );
			fractionNames = ProjectFile::getFractionNameList ( pList );
			numMSSpectra = ProjectFile::getNumMSSpectra ( pList );
			if ( sresMods ) SiteScores::init ( pList );	// Initialise from list of used mods
			rawTypes = ProjectFile::getRawTypes ( pList );
			if ( rawTypes.empty () ) rawTypes.resize ( fractionNames.size () );
			parentToleranceInfo = new ToleranceInfo ( "msms_parent_mass", pList );
			sysError = pList->getDoubleValue ( "msms_parent_mass_systematic_error" );
			sysErrorStr = pList->getStringValue ( "msms_parent_mass_systematic_error" );
			offsets = ProjectFile::getOffsets ( pList, sysError );
			fragmentToleranceInfo = new ToleranceInfo ( "fragment_masses", pList );
			linkInfo = new LinkInfo ( pList );

			databaseResults = new DatabaseResults ( pList->getStringVectorValue ( "database" ), ist );
		}
		ujm->writeMessage ( cout, "Reading results for each spectrum" );
		XMLIStreamList xstr ( ist, "d" );
		int num = 0;
		for ( ; ; ) {
			string spec;
			if ( xstr.getNextBlock ( spec ) ) {
				num++;
				if ( num % 10000 == 0 ) ujm->writeMessage ( cout, "Reading results for each spectrum, " + gen_itoa ( num ) + " spectra" );
				readBlock ( spec, mapTagHits, multisample, idFilterSet, idFilter, eValueFlag, mievi, expName, linearTailFitExpectation, minPeptideScore, maxPeptideEValue, xlMinLowScore, xlMinScoreDiff, xlMaxLowExpectation, iv );
			}
			else break;
		}
		if ( curFile == fname ) break;		// Just one results file
	}
	//sortCLinkPeptideLines ( "" );
	return instrument;
}
void SearchResults::readBlock ( string& spec, MapIDMapAccNoAndVectorSearchResultsPeptideHit& mapTagHits, bool multisample, const SetInt& idFilterSet, bool idFilter, bool eValueFlag, MapIDExpectationValueInfo& mievi, const string& expName, bool linearTailFitExpectation, double minPeptideScore, double maxPeptideEValue, double xlMinLowScore, double xlMinScoreDiff, double xlMaxLowExpectation, const IntVector& iv )
{
	bool idFilterSetFlag = !idFilterSet.empty ();
	string::size_type start = 0;
	string::size_type end;
	string specID = genNextString ( spec, "\t", start, end );
	string spotID;
	SpecID sid ( specID );
	if ( multisample ) {
		if ( spottingPlate ) spotID = sid.getSpotID ();
		else spotID = gen_itoa ( sid.getFraction () );
	}
	else {
		spotID = defaultID;
		if ( idFilterSetFlag ) {
			idFilter = idFilterSet.find ( sid.getFraction () ) != idFilterSet.end ();  
		}
	}
	idSet.insert ( spotID );
	int numPeaks = genNextInt ( spec, "\t", start, end );
	double mOverZ = genNextDouble ( spec, "\t", start, end );
	int charge = genNextInt ( spec, "\t", start, end );
	double intensity = genNextDouble ( spec, "\r\n", start, end );
	if ( spec [start] == '\n' ) start++;
	MSMSSpectrumInfo* mmsi = new MSMSSpectrumInfo ( mOverZ, charge, intensity, numPeaks );
	midmsi [spotID][specID] = mmsi;
	double a = 0.0;
	double b = 0.0;
	int numSpectra = 0;
	if ( !noExpectation ) {
		a = genNextDouble ( spec, "\t", start, end );
		b = genNextDouble ( spec, "\t", start, end );
		numSpectra = genNextInt ( spec, "\r\n", start, end );
		if ( spec [start] == '\n' ) start++;
		if ( eValueFlag ) {
			MapIDExpectationValueInfo::const_iterator cur = mievi.find ( specID );
			if ( cur != mievi.end () ) {
				a = (*cur).second.getGradient ();
				b = (*cur).second.getOffset ();
			}
			else {
				ErrorHandler::genError ()->error ( "The expectation value file " + expName + " is incomplete.\n" + "SpecID: " + specID + " not found.\n" );
			}
		}
	}
	int rank;
	int numUnmatched;
	double error;
	string nTerm;
	string peptide;
	string cTerm;
	string neutralLoss;
	int startAA;
	double score;
	double maxScore = -std::numeric_limits<double>::max();
	double expectation;
	string scoreDiff;
	string accNum;
	SpecID* spID = 0;
	HitPeptide* hitPeptide = 0;
	PeptideSpectralInfo* psi;
	SCModInfo scmi;
	bool diskFlag = discScoreGraph && !discFilenameExists;
	for ( ; ; ) {
		string s = genNextString ( spec, "\t", start, end );
		if ( end == string::npos ) break;
		bool mod = false;
		if ( s [0] == '+' ) {
			startAA = atoi ( s.substr ( 1 ).c_str () );
			accNum = genNextString ( spec, "\r\n", start, end );
			if ( spec [start] == '\n' ) start++;
		}
		else {
			if ( s [0] == 'M' ) {
				numUnmatched = genNextInt ( spec, "\t", start, end );
				getPeptideFromResults ( spec, start, nTerm, peptide, cTerm, neutralLoss );
				score = genNextDouble ( spec, "\r\n", start, end );
				mod = true;
			}
			else if ( s [0] == 'X' ) {
				readXLinkBlock ( spec, s, start, end, maxScore, a, b, numSpectra, linearTailFitExpectation, minPeptideScore, maxPeptideEValue, xlMinLowScore, xlMinScoreDiff, xlMaxLowExpectation, idFilter, spID, specID, numPeaks, mmsi );
				break;
			}
			else {
				rank = atoi ( s.c_str () );
				numUnmatched = genNextInt ( spec, "\t", start, end );
				error = genNextDouble ( spec, "\t", start, end );
				getPeptideFromResults ( spec, start, nTerm, peptide, cTerm, neutralLoss );
				startAA = genNextInt ( spec, "\t", start, end );
				score = genNextDouble ( spec, "\t", start, end );
				scoreDiff = genNextString ( spec, "\t", start, end );
				accNum = genNextString ( spec, "\r\n", start, end );
				maxScore = genMax ( maxScore, score );
			}
			expectation = getEval ( score, a, b, numSpectra, linearTailFitExpectation );
			if ( spec [start] == '\n' ) start++;
			if ( !mod ) {
				bool scoreFlag = score >= minPeptideScore && expectation <= maxPeptideEValue;
				if ( diskFlag || scoreFlag ) {
					hitPeptide = new HitPeptide ( nTerm, peptide, cTerm, neutralLoss );
					psi = new PeptideSpectralInfo ( numUnmatched, numPeaks, rank, score, expectation, scoreDiff, numSpectra, a, b );
				}
			}
			scmi.addHit ( nTerm, peptide, cTerm, neutralLoss, score, expectation, numSpectra );
		}
		if ( !mod ) {
			bool scoreFlag = score >= minPeptideScore && expectation <= maxPeptideEValue && idFilter;
			if ( diskFlag || scoreFlag ) {
				if ( spID == 0 ) spID = new SpecID ( specID );
				accNum = getFullAccNum ( iv, accNum );
				mapTagHits [spotID][accNum].push_back ( new SearchResultsPeptideHit ( spID, mmsi, psi, error, hitPeptide, startAA, accNum, numSearches-1 ) );
			}
		}
	}
	if ( spID ) scmi.add ( spID );
}
void SearchResults::readXLinkBlock ( string& spec, string& s, string::size_type& start, string::size_type& end, double maxScore, double a, double b, int numSpectra, bool linearTailFitExpectation, double minPeptideScore, double maxPeptideEValue, double xlMinLowScore, double xlMinScoreDiff, double xlMaxLowExpectation, bool idFilter, SpecID* spID, const string& specID, int numPeaks, MSMSSpectrumInfo* mmsi )
{
	bool first = true;
	for ( ; ; ) {
		if (!first ) {
			s = genNextString ( spec, "\t", start, end );
			if ( end == string::npos ) break;
		}
		char letter = s [0];
		int rank = atoi ( s.substr (1).c_str () );
		int numUnmatched = genNextInt ( spec, "\t", start, end );
		double error = genNextDouble ( spec, "\t", start, end );
		string nTerm;
		string peptide;
		string cTerm;
		string neutralLoss;
		getPeptideFromResults ( spec, start, nTerm, peptide, cTerm, neutralLoss );
		int startAA = genNextInt ( spec, "\t", start, end );
		string accNum = genNextString ( spec, "\t", start, end );
		string nTerm2;
		string peptide2;
		string cTerm2;
		string neutralLoss2;
		getPeptideFromResults ( spec, start, nTerm2, peptide2, cTerm2, neutralLoss2 );
		int startAA2 = genNextInt ( spec, "\t", start, end );
		string accNum2 = genNextString ( spec, "\t", start, end );
		double score;
		double xScore1;
		int xRank1 = 0;
		double xScore2;
		int xRank2 = 0;
		int dummy;										// This field used to hold the number of crosslinking ions. On searches after v5.14.1 is set to zero.
		double xFirstScore;
		XLinkDecoyInfo xldi;
		score = genNextDouble ( spec, "\t", start, end );
		xScore1 = genNextDouble ( spec, "\t", start, end );
		xRank1 = genNextInt ( spec, "\t", start, end );
		xScore2 = genNextDouble ( spec, "\t", start, end );
		xRank2 = genNextInt ( spec, "\t", start, end );
		xFirstScore = genNextDouble ( spec, "\t", start, end );
		dummy = genNextInt ( spec, "\r\n", start, end );
		if ( spec [start] == '\n' ) start++;
		xldi.add ( accNum, startAA, accNum2, startAA2 );
		for ( ; ; ) {
			string::size_type saveStart = start;
			string::size_type saveEnd = end;
			string t = genNextString ( spec, "\t", start, end );
			if ( end == string::npos ) break;
			if ( t [0] == '+' ) {
				startAA = atoi ( t.substr ( 1 ).c_str () );
				accNum = genNextString ( spec, "\t", start, end );
				startAA2 = genNextInt ( spec, "\t", start, end );
				accNum2 = genNextString ( spec, "\r\n", start, end );
				if ( spec [start] == '\n' ) start++;
				xldi.add ( accNum, startAA, accNum2, startAA2 );
			}
			else {
				start = saveStart;
				end = saveEnd;
				break;
			}
		}
		string scoreDiff = gen_ftoa ( score - maxScore, "%.1f" );
		double expectation = getEval ( score, a, b, numSpectra, linearTailFitExpectation );
		if ( ( rank == 1 || sresKeepReplicates ) && score >= minPeptideScore && expectation <= maxPeptideEValue && score-xFirstScore >= xlMinLowScore && ( scoreDiff == "---" || atof ( scoreDiff.c_str () ) >= xlMinScoreDiff ) && idFilter ) {
			if ( spID == 0 ) {
				spID = new SpecID ( specID );
			}
			HitPeptide* hitPeptide = new HitPeptide ( nTerm, peptide, cTerm, neutralLoss );
			HitPeptide* hitPeptide2 = new HitPeptide ( nTerm2, peptide2, cTerm2, neutralLoss2 );
			PeptideSpectralInfo* psi = new PeptideSpectralInfo ( numUnmatched, numPeaks, rank, score, expectation, scoreDiff, numSpectra, a, b );
			if ( psi->getExpectationValue ( psi->getScore () - xFirstScore ) <= xlMaxLowExpectation ) {
				xldi.calc ( 0 );
				cLinkPeptideHits.push_back ( new SearchResultsCrosslinkPeptideHit ( spID, mmsi, psi, error, hitPeptide, hitPeptide2, xldi, numSearches-1, xScore1, xRank1, xScore2, xRank2, xFirstScore ) );
			}
		}
		first = false;
	}
}
void SearchResults::sortCLinkPeptideLines ( const string& sortType )
{
	if ( !cLinkPeptideHits.empty () ) {
		sortCLinkPeptideLines ( sortType, cLinkPeptideHits.begin (), cLinkPeptideHits.end () );
		for ( int i = 0 ; i < cLinkPeptideHits.size () ; i++ ) {
			if ( i == 0 || cLinkPeptideHits [i]->getAccessionNumbers () != cLinkPeptideHits [i-1]->getAccessionNumbers () ) {
				cLinkProteinHits.push_back ( new SearchResultsCrosslinkProteinHit () );
			}
			cLinkProteinHits.back ()->push_back ( cLinkPeptideHits [i] );
		}
		stable_sort ( cLinkProteinHits.begin (), cLinkProteinHits.end (), sortCrosslinkProteinHitsByBestExpectationValue () );
		for ( int j = 0 ; j < cLinkProteinHits.size () ; j++ ) {
			if ( QuantitationRatio::getAreaRatioReport () )	cLinkProteinHits [j]->quanProteinStats ( true );
			if ( QuantitationRatio::getIntRatioReport () )	cLinkProteinHits [j]->quanProteinStats ( false );
		}
	}
}
void SearchResults::sortCLinkPeptideLines ( const string& sortType, const SearchResultsCrosslinkPeptideHitPtrVectorIterator& begin, const SearchResultsCrosslinkPeptideHitPtrVectorIterator& end )
{
	if		( sortType == "m/z" )				stable_sort ( begin, end, sortCrosslinkHitsByMOverZ () );
	else if	( sortType == "M+H" )				stable_sort ( begin, end, sortCrosslinkHitsByMPlusH () );
	else if	( sortType == "Error" )				stable_sort ( begin, end, sortCrosslinkHitsByError () );
	else if	( sortType == "Intensity" )			stable_sort ( begin, end, sortCrosslinkHitsByIntensity () );
	else if	( sortType == "Fraction/RT" )		stable_sort ( begin, end, sortCrosslinkHitsBySpot () );
	else if	( sortType == "RT" )				stable_sort ( begin, end, sortCrosslinkHitsByRT () );
	else if	( sortType == "Start Residue" )		stable_sort ( begin, end, sortCrosslinkHitsByStartResidue () );
	else if	( sortType == "End Residue" )		stable_sort ( begin, end, sortCrosslinkHitsByEndResidue () );
	else if	( sortType == "Crosslink AA" )		stable_sort ( begin, end, sortCrosslinkHitsByXLinkPosn () );
	else if	( sortType == "Peptide Score" )		stable_sort ( begin, end, sortCrosslinkHitsByScore () );
	else if	( sortType == "Expectation Value" )	stable_sort ( begin, end, sortCrosslinkHitsByExpectationValue () );
	else ErrorHandler::genError ()->error ( "Invalid sort type for crosslink report.\n" );
}
string SearchResults::getFullAccNum ( const IntVector& iv, const string& a )
{
	int n = gen_strcharcount ( a, '$' );
	if ( n == 0 || n == 2 )							// no database index
		return gen_itoa ( iv [0] ) + '$' + a;
	else {											// database index
		int pos = a.find ( '$' );
		int n = atoi ( a.substr ( 0, pos ).c_str () ) - 1;
		return gen_itoa ( iv [n] ) + '$' + a.substr ( pos+1 );
	}
}
double SearchResults::getEval ( double score, double a, double b, int numSpectra, bool linearTailFitExpectation ) const
{
	double e;
	if ( ( a == 0 && b == 0 ) || noExpectation ) {
		e = -1.0;
	}
	else {
		if ( linearTailFitExpectation )
			e = pow ( 10.0, (a * score) + b );
		else
			e = 1.0 - exp ( -exp ( (- 1 / b) * (score - a) ) );	// Method of moments
		e *= numSpectra;
		if ( isnan ( e ) ) e = -1.0;
	}
	return e;
}
void SearchResults::getPeptideFromResults ( string& spec, string::size_type& start, string& nTerm, string& peptide, string& cTerm, string& neutralLoss )
{
	string::size_type previousStart = start;
	string::size_type end;
	string s1 = genNextString ( spec, "\t", start, end );
	nTerm = "";
	if ( s1 [0] == '-' ) nTerm = s1.substr ( 1 );
	else start = previousStart;							// No N Terminus

	peptide = "";
	peptide = genNextString ( spec, "\t", start, end );	// Should always be there

	previousStart = start;
	string s2 = genNextString ( spec, "\t", start, end );
	cTerm = "";
	if ( s2.find_first_of ( "\r\n" ) == string::npos && s2 [0] == '-' ) cTerm = s2.substr ( 1 );	// This is required for M entries where a negative score can be misinterpreted as a c-term
	else start = previousStart;							// No C Terminus

	previousStart = start;
	string s3 = genNextString ( spec, "\t", start, end );
	neutralLoss = "";
	if ( s3 [0] == '+' ) neutralLoss = s3.substr ( 1 );
	else start = previousStart;							// No Neutral Loss
}
bool SearchResults::getEValues ( const string& fname, MapIDExpectationValueInfo& mievi )
{
	string fileName = fname;
	for ( int i = 0 ; ; i++ ) {
		string curFile;
		if ( genFileExists ( fname ) ) curFile = fname;	// All results in one file
		else {
			string f = fname + string ( "_" ) + gen_itoa ( i );
			if ( genFileExists ( f ) ) curFile = f;		// Results spread across several files
			else {
				if ( i == 0 )
					return false;						// No file(s) in the first place
				else									// Keep going until you run out of files
					break;
			}
		}
		if ( i == 0 ) ujm->writeMessage ( cout, "Getting expectation values" );
		GenIFStream ist ( curFile );						// Open the current file
		XMLIStreamList xslDataSet ( ist, "d" );
		for ( ; ; ) {
			string spec;
			if ( xslDataSet.getNextBlock ( spec ) ) {
				string::size_type start = 0;
				string::size_type end;
				string specID = genNextString ( spec, "\t", start, end );
				string dummy = genNextString ( spec, "\r\n", start, end );
				if ( spec [start] == '\n' ) start++;
				double gradient = genNextDouble ( spec, "\t", start, end );
				double offset = genNextDouble ( spec, "\t", start, end );
				int numSpectra = genNextInt ( spec, "\r\n", start, end );
				if ( spec [start] == '\n' ) start++;
				mievi [specID] = ExpectationValueInfo ( gradient, offset, numSpectra );
			}
			else break;
		}
		if ( curFile == fname ) break;						// Just one results file
	}
	return true;
}
void SearchResults::drawHistogram ( ostream& os ) const
{
	if ( discScoreGraph ) { 
		histogram.drawGraph ( os );
		if ( !genFileExists ( discFilename ) ) {	// Create cache file if it doesn't exist
			if ( histogram.size () ) {
				GenOFStream ofs ( discFilename );
				ofs << histogram;
			}
		}
	}
}
StringVector SearchResults::getMergedAccNumbers ( const string& id )
{
	StringVector sv;
	for ( MapAccNumSearchResultsProteinHitPtrConstIterator i = compProtHits [id].begin () ; i != compProtHits [id].end () ; i++ ) {
		if ( (*i).second->getNumMSMSHits () ) {
			sv.push_back ( (*i).first );
		}
	}
	stable_sort ( sv.begin (), sv.end () );
	sv.erase ( unique ( sv.begin (), sv.end () ), sv.end () );
	return sv;
}
void SearchResults::printCLinkProteinHitsLinks ( ostream& os, int idx ) const
{
	string xLinkType;
	string printedXLinkType;
	for ( int i = 0 ; i < cLinkProteinHits.size () ; i++ ) {
		xLinkType = cLinkProteinHits [i]->getXLinkType ();
		if ( xLinkType != printedXLinkType ) {
			os << "<a href=\"#" << gen_strstrip ( xLinkType ) << idx+1 << "\">";
			os << xLinkType << " Crosslinked Peptide Hits " << idx+1;
			os << "</a><br />" << endl;
			printedXLinkType = xLinkType;
		}
	}
}
void SearchResults::printCLinkProteinHitsHTML ( ostream& os, const SResLink& sresLink, int idx ) const
{
	string xLinkType;
	string printedXLinkType;
	for ( int i = 0 ; i < cLinkProteinHits.size () ; i++ ) {
		xLinkType = cLinkProteinHits [i]->getXLinkType ();
		if ( xLinkType != printedXLinkType ) {
			os << "<div class=\"results_section_header\">";
			os << "<a name=\"" << gen_strstrip ( xLinkType ) << idx+1 << "\">";
			os << xLinkType << " Crosslinked Peptide Hits " << idx+1;
			os << "</a>";
			os << "</div>" << endl;
			printedXLinkType = xLinkType;
		}
		cLinkProteinHits [i]->printHTML ( os, sresLink, idx, linkInfo );
		if ( i != cLinkProteinHits.size () - 1 ) os << "<hr />" << endl;
	}
}
void SearchResults::printCLinkProteinHitsDelimited ( ostream& os, int idx ) const
{
	if ( !cLinkProteinHits.empty () ) {
		cLinkProteinHits [0]->printHeaderDelimited ( os, idx );
		for ( int i = 0 ; i < cLinkProteinHits.size () ; i++ ) {
			cLinkProteinHits [i]->printDelimited ( os, linkInfo );
		}
	}
}
void SearchResults::setQuanCLink ( int idx )
{
	if ( !cLinkPeptideHits.empty () ) {
		sortCLinkPeptideLines ( "Fraction/RT", cLinkPeptideHits.begin (), cLinkPeptideHits.end () );
		int currentFraction = -1;
		int last = cLinkPeptideHits.size () - 1;
		int num;
		GenElapsedTime et;
		SetPairStringInt spsi;
		bool quanMSMSFlag = PeptidePositionQuan::getQuanMSMSFlag ();
		for ( std::vector <const SearchResultsCrosslinkPeptideHit*>::size_type k = 0 ; k < cLinkPeptideHits.size () ; k++ ) {
			int fraction = cLinkPeptideHits [k]->getFraction ();
			if ( fraction != currentFraction ) {
				ujm->writeMessage ( cout, "Processing quantitation results for fraction " + gen_itoa ( fraction ) );
				if ( currentFraction != -1 ) PeptidePositionQuan::deleteQuan ( idx );
				PeptidePositionQuan::initialiseQuan ( idx, fraction );							// Init new fraction
				currentFraction = fraction;
				num = 0;
				spsi.clear ();
			}
			num++;
			if ( et.getElapsedTime () > 10 ) {
				ujm->writeMessage ( cout, "Processing quantitation results for fraction " + gen_itoa ( fraction ) + ", " + gen_itoa ( num ) + " peptides" );
				et.reset ();
			}
			PairSetPairStringIntIteratorBool pspsiib = spsi.insert ( make_pair ( cLinkPeptideHits [k]->getAccessionNumbers () + QuantitationRatio::getUnmodifiedPeptide ( cLinkPeptideHits [k]->getSequence1 () ) + '\t' + QuantitationRatio::getUnmodifiedPeptide ( cLinkPeptideHits [k]->getSequence2 () ), cLinkPeptideHits [k]->getCharge () ) );
			if ( quanMSMSFlag || pspsiib.second ) {
				cLinkPeptideHits [k]->setQuanResults ( linkInfo );
			}
			if ( k == last ) PeptidePositionQuan::deleteQuan ( idx );
		}
	}
}
MSMSSpectrumInfo::MSMSSpectrumInfo ( double mOverZ, int charge, double intensity, int numPeaks ) :
	mOverZ ( mOverZ ),
	intensity ( intensity ),
	charge ( charge ),
	numPeaks ( numPeaks )
{
}
