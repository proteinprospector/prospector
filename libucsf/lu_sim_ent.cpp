/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_sim_ent.cpp                                                *
*                                                                             *
*  Created    : May 17th 2001                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <algorithm>
#include <lg_stdio.h>
#include <lg_string.h>
#include <lgen_file.h>
#include <lu_aa_calc.h>
#include <lu_aa_info.h>
#include <lu_acc_link.h>
#include <lu_acc_num.h>
#include <lu_sim_ent.h>
#include <lu_html_form.h>
#include <lu_html.h>
#include <lu_mass.h>
#include <lu_check_db.h>
#include <lu_pq_vector.h>
#include <lu_param_list.h>
#include <lu_cgi_val.h>
#include <lu_db_entry.h>
#include <lu_fasta.h>
#include <lu_getfil.h>
using std::vector;
using std::string;
using std::ostream;
using std::find;
using std::endl;
using std::istringstream;
using std::getline;
using std::stringstream;
using std::random_shuffle;
using std::reverse;
using std::make_pair;

class SingleDatabaseEntry : public SingleEntry {
	string database;
	StringVector names;
	StringVector species;
	StringVector accessionNumbers;
	StringVector uniprotIDs;
	StringVector organismNames;
	StringVector geneNames;
	StringVector proteinExistences;
	StringVector sequenceVersion;
public:
	SingleDatabaseEntry ( FastaServer& fs, int indexNumber, int dnaReadingFrame, int openReadingFrame, const string& protein );
	int getIndex () const { return indexNum; }
	void printHTML ( ostream& os ) const;
	void printXMLBody ( ostream& os ) const;
};

SingleEntryParameters::SingleEntryParameters ( const ParameterList* params, const string& tempDatabase ) :
	tempDBFullPath ( tempDatabase )
{
	database		= params->getStringValue	( "database", "" );
	accessMethod	= params->getStringValue	( "access_method", "Index Number" );
	maxNTermAA		= params->getIntValue		( "n_term_aa_limit", 0 );

	if ( tempDatabase.empty () ) {
		if ( database == "User Protein" )	parseForUserProtein ( params );
		else								parseForEntryData ( params );
	}

	if ( maxNTermAA ) {
		for ( StringVectorSizeType i = 0 ; i < userProteinSequence.size () ; i++ ) {
			if ( userProteinSequence [i].size () > maxNTermAA ) {
				userProteinSequence [i].resize ( maxNTermAA );
			}
		}
	}
}
void SingleEntryParameters::parseForUserProtein ( const ParameterList* params )
{
	StringVector sv = params->getPQStringVectorValue ( "user_protein_sequence" );
	bool dbSearch = params->getStringValue	( "access_method" ) == "";				// This is a database search program
	string prot;
	string n;
	string name;
	for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
		bool cLine = false;
		if ( sv [i][0] == '>' ) {
			n = sv [i].substr (1);
			cLine = true;
		}
		else {
			prot += gen_strstrip ( sv [i] );	// Strip out all spaces
			name = n;
		}
		if ( !prot.empty () && cLine ) {
			if ( dbSearch )	addDBUserProteinSequence ( prot );
			else			addUserProteinSequence ( prot );
			names.push_back ( name );
			prot = "";
		}
	}
	if ( !prot.empty () ) {
		if ( dbSearch )	addDBUserProteinSequence ( prot );
		else			addUserProteinSequence ( prot );
		names.push_back ( name );
	}
	if ( userProteinSequence.empty () ) {
		ErrorHandler::genError ()->error ( "No user protein sequence entered." );
	}
}
void SingleEntryParameters::parseForEntryData ( const ParameterList* params )
{
	const char* value;
	if ( params->getValue ( "entry_data", value ) ) {
		if ( is_dna_database ( database ) ) {
			if ( accessMethod == "Index Number" ) {
				bool flag = getPostQuery3Vectors ( value, indexNum, dnaReadingFrame, openReadingFrame );
				if ( !flag ) {
					ErrorHandler::genError ()->error ( "Incorrectly formatted List of Entries for DNA database.\n Correct format is:\n(Index Number) (DNA Reading Frame) (Open Reading Frame)\n" );
				}
			}
			else {
				bool flag = getPostQuery3Vectors ( value, accessionNum, dnaReadingFrame, openReadingFrame );
				if ( !flag ) {
					ErrorHandler::genError ()->error ( "Incorrectly formatted List of Entries for DNA database.\n Correct format is:\n(Accession Number) (DNA Reading Frame) (Open Reading Frame)\n" );
				}
			}
		}
		else {
			if ( accessMethod == "Index Number" ) {
				getPostQueryVector ( value, indexNum, '\n' );
				if ( indexNum.empty () ) {
					StringVector temp;
					getPostQueryVector ( value, temp, '\n' );
					if ( temp.empty () )
						ErrorHandler::genError ()->error ( "No index numbers entered.\n" );
					else
						ErrorHandler::genError ()->error ( "You have not entered index numbers which must be integer numbers. If you entered accession numbers change Retrieve Entry by to Accession Number.\n" );
				}
			}
			else {
				getPostQueryVector ( value, accessionNum, '\n' );
				if ( accessionNum.empty () )
					ErrorHandler::genError ()->error ( "No accession numbers entered.\n" );
			}
		}
		BoolDeque bv = params->getBoolDequeValue ( "e", getNumEntries () );
		BoolDequeConstIterator bvi = find ( bv.begin (), bv.end (), true );		// If all are false then no e parameters so don't delete anything.
		if ( bvi != bv.end () ) {
			deleteEntries ( bv, accessionNum );
			deleteEntries ( bv, indexNum );
			deleteEntries ( bv, dnaReadingFrame );
			deleteEntries ( bv, openReadingFrame );
		}
	}
}
int SingleEntryParameters::getNumEntries () const
{
	if ( !userProteinSequence.empty () ) {
		return userProteinSequence.size ();
	}
	else {
		return genMax ( accessionNum.size (), indexNum.size () );
	}
}
void SingleEntryParameters::addDBUserProteinSequence ( const string& value )
{
	string ups = gen_strstrip ( value );	// Strip out all spaces
	ups = gen_strstripnum ( ups );			// Strip out all numbers
	ups = gen_strstripchar ( ups, '*' );	// Strip out all asterisks
	ups = gen_strstripchar ( ups, '/' );	// Strip out all / characters
	ups = gen_strstripchar ( ups, '^' );	// Strip out all ^ characters
	if ( ups.empty () ) return;
	if ( ups [ups.length () - 1] == '.' ) ups = ups.substr ( 0, ups.length () - 1 );	// Delete a trailing dot
	if ( ups [ups.length () - 1] == '-' ) ups = ups.substr ( 0, ups.length () - 1 );	// Delete a trailing dash
	if ( ups.empty () ) return;
	if ( !genIsUpper ( ups ) ) {
		if ( genIsLower ( ups ) ) {
			ups = genToUpper ( ups );	// Try converting to upper case if all lower case
		}
		else {
			ups = convertDBSequence ( ups );	// Check if user entered 3 letter code
			if ( ups.empty () ) {
				ErrorHandler::genError ()->error ( "Invalid user protein sequence entered." );
			}
		}
	}
	userProteinSequence.push_back ( ups );
}
void SingleEntryParameters::addUserProteinSequence ( const string& value )
{
	static bool prot = false;
	static bool dna = false;
	static int numAdded = 0;
	string ups = gen_strstrip ( value );	// Strip out all spaces
	ups = gen_strstripnum ( ups );			// Strip out all numbers
	ups = gen_strstripchar ( ups, '*' );	// Strip out all asterisks
	ups = gen_strstripchar ( ups, '/' );	// Strip out all / characters
	ups = gen_strstripchar ( ups, '^' );	// Strip out all ^ characters
	if ( ups.empty () ) return;
	if ( ups [ups.length () - 1] == '.' ) ups = ups.substr ( 0, ups.length () - 1 );	// Delete a trailing dot
	if ( ups [ups.length () - 1] == '-' ) ups = ups.substr ( 0, ups.length () - 1 );	// Delete a trailing dash
	if ( ups.empty () ) return;
	string invalidAA = checkSequence ( ups );	// This is the invalid AA to report
	bool flag = invalidAA.empty ();
	if ( !flag ) {
		if ( genIsLower ( ups ) ) ups = genToUpper ( ups );	// Try converting to upper case if all lower case
		string invAA2 = checkSequence ( ups );
		flag = invAA2.empty ();
		if ( !flag ) {
			ups = convertSequence ( ups );	// Check if user entered 3 letter code
			flag = !ups.empty ();
		}
	}
	if ( !flag ) {
		ErrorHandler::genError ()->error (
			"The input sequence contains the invalid amino acid: '" + invalidAA + "'" );
	}
	int len = ups.length ();
	if ( len > 10 && static_cast<double> ( genNumChars ( ups, "ACGT" ) ) / static_cast<double> ( len ) > 0.85 ) {
		if ( prot ) {
			ErrorHandler::genError ()->error (
				"DNA Sequences cannot be mixed with protein sequences." );
		}
		int numUnknowns;
		for ( int i = 1 ; i <= 6 ; i++ ) {
			char* prot = new char [len]; 
			readProteinFromDNA ( i, maxNTermAA, const_cast<char*>(ups.c_str ()), len, prot, numUnknowns );
			string sProt ( prot );
			int st = 0;
			int orf = 0;
			for ( ; ; ) {
				int end = sProt.find ( '.', st );
				string s = ( end == string::npos ) ? sProt.substr ( st ) : sProt.substr ( st, end-st );
				if ( !s.empty () ) {
					userProteinSequence.push_back ( s );
					orf++;
					numAdded++;
					indexNum.push_back ( numAdded );
					dnaReadingFrame.push_back ( i );
					openReadingFrame.push_back ( orf );
				}
				if ( end == string::npos ) break;
				st = end + 1;
			}
			delete [] prot;
			
		}
		dna = true;
	}
	else {
		if ( dna ) {
			ErrorHandler::genError ()->error (
				"DNA Sequences cannot be mixed with protein sequences." );
		}
		userProteinSequence.push_back ( ups );
		prot = true;
	}
}
string SingleEntryParameters::checkSequence ( const string& s ) const
{
	for ( StringSizeType i = 0 ; i < s.length () ; i++ ) {
		char c = s [i];
		if ( !aa_composition_mask [c] && c != 'B' && c != 'X' && c != 'Z' ) {
			return string ( 1, c );
		}
	}
	return "";
}
string SingleEntryParameters::convertDBSequence ( const string& s ) const
{
	string ret;
	for ( StringSizeType i = 0 ; i < s.length () ; i += 3 ) {
		char c = AAInfo::getInfo ().getAA ( s.substr ( i, 3 ) );
		if ( c == 0 ) return "";
		ret += c;
	}
	return ret;
}
string SingleEntryParameters::convertSequence ( const string& s ) const
{
	string ret;
	for ( StringSizeType i = 0 ; i < s.length () ; i += 3 ) {
		char c = AAInfo::getInfo ().getAA ( s.substr ( i, 3 ) );
		if ( c == 0 || ( !aa_composition_mask [c] && c != 'B' && c != 'X' && c != 'Z' ) ) {
			return "";
		}
		ret += c;
	}
	return ret;
}
void SingleEntryParameters::printHTML ( ostream& os ) const
{
	ParameterList::printHTML ( os, "Database", database );
	if ( userProteinSequence.empty () ) {
		for ( int i = 0 ; i < getNumEntries () ; i++ ) {
			if ( accessMethod == "Index Number" )
				ParameterList::printHTML ( os, "Index Number", indexNum [i] );
			else
				ParameterList::printHTML ( os, "Accession Number", accessionNum [i] );
			if ( is_dna_database ( database ) ) {
				ParameterList::printHTML ( os, "DNA Reading Frame", dnaReadingFrame [i] );
				ParameterList::printHTML ( os, "Open Reading Frame", openReadingFrame [i] );
			}
		}
	}
}
void SingleEntryParameters::putCGI ( ostream& os ) const
{
	printCGIString ( os, "database", database );
}
void SingleEntryParameters::putIndexCGI ( ostream& os, bool dnaFlag, int indexNum, int dnaReadingFrame, int openReadingFrame )
{
	printCGIString ( os, "access_method", "Index Number" );
	stringstream s;
	s << indexNum;
	if ( dnaFlag )
		s << " " << dnaReadingFrame << " " << openReadingFrame + 1;

	printCGIString ( os, "entry_data", s.str () );
}
void SingleEntryParameters::putIndexHiddenFormEntry ( ostream& os, bool dnaFlag, int indexNum, int dnaReadingFrame, int openReadingFrame )
{
	printHTMLFORMHidden ( os, "access_method", "Index Number" );
	stringstream s;
	s << indexNum;
	if ( dnaFlag )
		s << " " << dnaReadingFrame << " " << openReadingFrame + 1;

	printHTMLFORMHidden ( os, "entry_data", s.str () );
}
void SingleEntryParameters::putHiddenFormEntry ( ostream& os ) const
{
	printHTMLFORMHidden ( os, "database", database );
}
void SingleEntryParameters::copyToHiddenFormEntry ( ostream& os, const ParameterList* params )
{
	params->copyToHiddenFormEntry ( os, "database" );
}
PairStringBool SingleEntryParameters::createTemporaryDatabase ( bool randomFlag, bool reverseFlag, const string& tempDir ) const
{
	static string key = genToLower ( genRandomString ( 10 ) );
	PPTempFile pptf ( "", "" );					// Move the directory to the temporary directory
	tempDBFullPath = pptf.getFullPathDir () + SLASH;

	if ( tempDir.empty () )	{
		tempDBFullPath += key;
	}
	else
		tempDBFullPath += tempDir;

	bool newDB;
	if ( genCreateNewDirectory ( tempDBFullPath ) ) {
		tempDBFullPath += SLASH;
		tempDBFullPath += "UserProtein.fasta";
		FILE* fp = gen_fopen_binary ( tempDBFullPath, "w", "Creating temporary database" );
		for ( int i = 0 ; i < getNumEntries () ; i++ ) {
			writeProteinEntry ( fp, names [i], userProteinSequence [i] );
		}
		gen_fclose ( fp, "Creating temporary database" );
		if ( randomFlag || reverseFlag ) {
			string suffix;
			if ( randomFlag ) suffix = ".random.fasta";
			if ( reverseFlag ) suffix = ".reverse.fasta";
			FILE* fp = gen_fopen_binary ( tempDBFullPath.substr ( 0, tempDBFullPath.length () - 6 ) + suffix, "w", "Creating temporary database" );
			for ( int i = 0 ; i < getNumEntries () ; i++ ) {
				string p = userProteinSequence [i];
				if ( randomFlag )	random_shuffle ( p.begin (), p.end () );
				else				reverse ( p.begin (), p.end () );
				writeProteinEntry ( fp, names [i], p );
			}
			gen_fclose ( fp, "Creating temporary database" );
		}
		newDB = true;
	}
	else {
		newDB = false;
		tempDBFullPath += SLASH;
		tempDBFullPath += "UserProtein.fasta";
	}
	return make_pair ( tempDBFullPath, newDB );
}
void SingleEntryParameters::writeProteinEntry ( FILE* fp, const string& n, const string& protein )
{
	fputc ( '>', fp );
	if ( !n.empty () ) gen_fwrite ( const_cast <char*>(n.c_str ()), n.length () * sizeof (char), 1, fp, "Creating temporary database" );
	fputc ( '\n', fp );
	gen_fwrite ( const_cast <char*>(protein.c_str ()), protein.length () * sizeof (char), 1, fp, "Creating temporary database" );
	fputc ( '\n', fp );
}
vector <SingleEntry*> getSingleEntry ( const SingleEntryParameters& p )
{
	AccessionNumberMap* am;
	if ( p.getDatabase () != "User Protein" && p.getAccessMethod () == "Accession Number" ) {
		am = getAccessionNumberMap ( p.getDatabase () );
	}
	vector <SingleEntry*> singleEntries;
	for ( int i = 0 ; i < p.getNumEntries () ; i++ ) {
		if ( p.getDatabase () == "User Protein" ) {
			singleEntries.push_back ( new SingleEntry ( p.getUserProteinSequence ( i ), p.getName ( i ), p.getDNADatabase (), i+1, p.getDNAReadingFrame ( i ), p.getOpenReadingFrame ( i ) ) );
		}
		else {
			int indexNum;
			if ( p.getAccessMethod () == "Accession Number" ) {
				indexNum = am->getIndexNumber ( p.getAccessionNum ( i ).c_str () );
				if ( indexNum == -1 ) {
					ErrorHandler::genError ()->error ( "The accession number \"" + p.getAccessionNum ( i ) + "\" does not exist in the database.\n" );
				}
			}
			else
				indexNum = p.getIndexNum ( i );
			FastaServer fs ( p.getDatabase () );
			fs.setMaxNTermAA ( p.getMaxNTermAA () );
			DatabaseEntry de ( indexNum, p.getDNAReadingFrame ( i ), p.getOpenReadingFrame ( i ) - 1 );
			singleEntries.push_back ( new SingleDatabaseEntry ( fs, indexNum, p.getDNAReadingFrame (i), p.getOpenReadingFrame (i), fs.getProtein ( de ) ) );
		}
	}
	if ( p.getDatabase () != "User Protein" && p.getAccessMethod () == "Accession Number" ) {
		delete am;
	}
	return singleEntries;
}
SingleDatabaseEntry::SingleDatabaseEntry ( FastaServer& fs, int indexNumber, int dnaReadingFrame, int openReadingFrame, const string& protein ) :
	SingleEntry ( protein, "", fs.getDNADatabase (), indexNumber, dnaReadingFrame, openReadingFrame ),
	database ( fs.getFileName () )
{
	for ( fs.firstLine ( indexNum ) ; fs.isDoneLine () ; fs.nextLine () ) {
		names.push_back				( fs.getLineName () );
		species.push_back			( fs.getLineSpecies () );
		accessionNumbers.push_back	( fs.getLineAccessionNumber () );
		uniprotIDs.push_back		( fs.getLineUniprotID () );
		organismNames.push_back		( fs.getLineOrganismName () );
		geneNames.push_back			( fs.getLineGeneName () );
		proteinExistences.push_back	( fs.getLineProteinExistence () );
		sequenceVersion.push_back	( fs.getLineSequenceVersion () );
	}
}
void SingleDatabaseEntry::printHTML ( ostream& os ) const
{
	startJavascript ( os );
	AccessionNumberLinkInfo anli;
	anli.printHTML ( os, database );
	endJavascript ( os );
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		os << "<b>Acc. #: </b>";
		anli.write2 ( os, accessionNumbers [i], true );
		os << " ";
		if ( !uniprotIDs [i].empty () ) {
			os << "<b>Uniprot ID: </b>";
			os << uniprotIDs [i];
			os << " ";
		}
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
	SingleEntry::printHTML ( os );
}
void SingleDatabaseEntry::printXMLBody ( ostream& os ) const
{
	for ( StringVectorSizeType i = 0 ; i < names.size () ; i++ ) {
		os << "<comment_line>" << endl;
		ParameterList::printXML ( os, "protein_name", names [i] );
		ParameterList::printXML ( os, "species", species [i] );
		ParameterList::printXML ( os, "accession_number", accessionNumbers [i] );
		os << "</comment_line>" << endl;
	}
	SingleEntry::printXMLBody ( os );
}
SingleEntry::SingleEntry ( const string& protein, const string& singleEntryName, bool dnaDatabase, int indexNumber, int dnaReadingFrame, int openReadingFrame ) :
	prot ( protein ),
	singleEntryName ( singleEntryName ),
	ppi ( prot.c_str () ),
	pmw ( prot.c_str () ),
	aaComposition ( AACalculator::calculateAAComposition ( prot ).getFormula () ),
	dnaDatabase ( dnaDatabase ),
	indexNum ( indexNumber ),
	drf ( dnaReadingFrame ),
	orf ( openReadingFrame )
{
}
void SingleEntry::printHTML ( ostream& os ) const
{
	ParameterList::printHTML ( os, "Index Number", indexNum );
	if ( dnaDatabase ) {
		ParameterList::printHTML ( os, "DNA Reading Frame", drf );
		ParameterList::printHTML ( os, "Open Reading Frame Number", orf );
	}
	ParameterList::printHTML ( os, "pI of Protein", ppi );
	ParameterList::printHTML ( os, "Protein MW", pmw );
	ParameterList::printHTML ( os, "Amino Acid Composition", aaComposition );
	if ( !singleEntryName.empty () ) ParameterList::printHTML ( os, "Name", singleEntryName );
	os << "<br />" << endl;
}
void SingleEntry::printXML ( ostream& os ) const
{
	os << "<database_entry>" << endl;
	printXMLBody ( os );
	os << "</database_entry>" << endl;
}
void SingleEntry::printXMLBody ( ostream& os ) const
{
	ParameterList::printXML ( os, "index_number", indexNum );
	if ( dnaDatabase ) {
		ParameterList::printXML ( os, "dna_reading_frame", drf );
		ParameterList::printXML ( os, "open_reading_frame", orf );
	}
	ParameterList::printXML ( os, "protein_pi", ppi );
	ParameterList::printXML ( os, "protein_mw", pmw );
	ParameterList::printXML ( os, "protein_aa_composition", aaComposition );
}
