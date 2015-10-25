/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pre_srch.cpp                                               *
*                                                                             *
*  Created    : July 18th 2000                                                *
*                                                                             *
*  Purpose    : Functions to do all pre-searches.                             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2000-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <algorithm>
#include <lgen_file.h>
#include <lgen_reg_exp.h>
#include <lu_fas_ind.h>
#include <lu_fas_sp.h>
#include <lu_fasta.h>
#include <lu_acc_num.h>
#include <lu_pq_vector.h>
#include <lu_pre_srch.h>
#include <lu_html.h>
#include <lu_html_form.h>
#include <lu_cgi_val.h>
#include <lu_param_list.h>
#include <lu_species.h>
#include <lu_table.h>
#include <lu_check_db.h>

#ifdef VIS_C
#pragma warning( disable : 4503 )	// Disable warnings about decorated name length being too long
#endif

using std::string;
using std::ostream;
using std::endl;
using std::vector;
using std::set_intersection;
using std::set_difference;
using std::sort;
using std::stable_sort;
using std::ios;
using std::cout;
using std::find;
using std::back_inserter;
using std::make_pair;
using std::unique;

struct TaxonomyHit {
	string scientificName;
	int node;
	int numEntries;
	TaxonomyHit ( const string& scientificName, int node, int numEntries ) :
		scientificName ( scientificName ),
		node ( node ),
		numEntries ( numEntries ) {}
};

class SortTaxonomyAscending {
	public:
		int operator () ( const TaxonomyHit& a, const TaxonomyHit& b ) const
		{
			return genStrcasecmp ( a.scientificName, b.scientificName ) < 0;
		}
};

class PreSearch {
protected:
	string database;
	bool concat;
	int numEntries;
	IntVector finalIndicies;
	string getMessage () const;
	static string FULL_MW_RANGE;
	static string FULL_PI_RANGE;
	static string SPECIES;
	static void getIndiciesFromAccessionNumbers ( const string& filename, const StringVector& accessionNumbers, IntVector& indicies );
public:
	PreSearch ( const string& database, bool concat = false );
	virtual ~PreSearch () = 0;
	virtual void doSearch ( FastaServer* fs ) = 0;
	virtual void putCGI ( ostream& os ) const = 0;
	virtual void putHiddenFormEntryJavascript ( ostream& os ) const = 0;
	virtual void printHTML ( ostream& os ) const = 0;
	virtual void printXMLResults ( ostream& os ) const;
	const IntVector& getIndicies () const { return finalIndicies; };
	int getNumIndicies () const { return finalIndicies.size (); };
	static StringVector addConcatAccessionNumbers ( const StringVector& accessionNumbers );
};

class IndexPreSearch : public PreSearch {
	string inputFilename;
	string inputProgramName;
	string inputDirectory;
	string databaseDate;
public:
	IndexPreSearch ( const ParameterList* params, const string& database, int dbIdx );
	void doSearch ( FastaServer* fs );
	void printHTML ( ostream& os ) const;
	void putCGI ( ostream& os ) const;
	void putHiddenFormEntryJavascript ( ostream& os ) const;
};

class AccessionPreSearch : public PreSearch {
	StringVector accessionNumbers;
	static string ACCESSION_NUMS;
public:
	AccessionPreSearch ( const string& database, const StringVector& accessionNumbers );
	void doSearch ( FastaServer* fs );
	void printHTML ( ostream& os ) const;
	void putCGI ( ostream& os ) const;
	void putHiddenFormEntryJavascript ( ostream& os ) const;
};

class UserProteinPreSearch : public PreSearch {
	string prefix;
public:
	UserProteinPreSearch ( const string& prefix );
	void doSearch ( FastaServer* fs );
	void printHTML ( ostream& os ) const;
	void putCGI ( ostream& os ) const {}
	void putHiddenFormEntryJavascript ( ostream& os ) const;
};

class StandardPreSearch : public PreSearch {
	string prefix;

	bool fullMWRange;
	double lowMass;
	string lowMassStr;
	double highMass;
	string highMassStr;
	bool fullPIRange;
	double lowPI;
	string lowPIStr;
	double highPI;
	string highPIStr;
	StringVector names;
	StringVector species;
	bool speciesRemove;
	StringVector speciesNames;
	StringVector addAccNums;

	StringVector fullTaxonomyList;

	int numMWIndicies;
	int numPIIndicies;
	int numTaxIndicies;
	int numPreNameIndicies;

	vector <TaxonomyHit> taxonomyHits;

	static string PROT_LOW_MASS;
	static string PROT_HIGH_MASS;
	static string LOW_PI;
	static string HIGH_PI;
	static string SPECIES_REMOVE;
	static string SPECIES_NAMES;
	static string ADD_ACC_NUMS;
	static string NAMES;

	static bool reportTaxonomy;

	void floatingPointListSearch ( const string& filename, double low, double high, IntVector& serialNumber );
	void name_regexp_pre_search ( FastaServer* fs, IntVector& indicies, const StringVector& names );
	void outputHTMLSpeciesSearchResults ( ostream& os ) const;
public:
	StandardPreSearch ( const ParameterList* params, const string& database, const StringVector& addAccNums, const string& prefix = "" );
	void doSearch ( FastaServer* fs );
	void printHTML ( ostream& os ) const;
	void printXMLResults ( ostream& os ) const;
	void putCGI ( ostream& os ) const;
	void putHiddenFormEntryJavascript ( ostream& os ) const;
	static void setReportTaxonomy ( bool flag ) { reportTaxonomy = flag; }  
};

string PreSearchInfo::RESULTS_FROM_FILE = "results_from_file";
bool PreSearchInfo::RESULTS_FROM_FILE_DEFAULT = false;

PreSearchInfo::PreSearchInfo ( const ParameterList* params, const string& prefix ) :
	fromFile ( params->getBoolValue ( RESULTS_FROM_FILE, RESULTS_FROM_FILE_DEFAULT ) )
{
	StringVector db = params->getStringVectorValue ( "database" );
	MapStringToStringVector accessionNumbers;
	MapStringToStringVector addAccessionNumbers;
	setAccessionNumbers ( accessionNumbers, params->getPQStringVectorValue ( "accession_nums" ) );
	setAccessionNumbers ( addAccessionNumbers, params->getPQStringVectorValue ( "add_accession_numbers" ) );
	DBSearchFlags dbFlags ( db );
	for ( int i = 0 ; i < db.size () ; i++ ) {
		string menuDatabase = db [i];
		preSearchIndex.push_back ( i );
		if ( menuDatabase == "User Protein" ) {
			preSearch.push_back ( new UserProteinPreSearch ( prefix ) );
			if ( dbFlags.getConcatFlag () ) preSearchIndex.push_back ( i );
		}
		else {
			PairStringString pss;
			int newSize = preSearchIndex.size();
			string searchDatabase = menuDatabase;
			if ( getConcatDBPair ( menuDatabase, pss ) ) {
				searchDatabase = pss.first;
				preSearchIndex.push_back ( i );
			}
			if ( fromFile ) {
				preSearch.push_back ( new IndexPreSearch ( params, searchDatabase, newSize ) );
			}
			else {
				PairBoolStringVector pbsv = getAccNumSearch ( menuDatabase, accessionNumbers );
				if ( pbsv.first )
					preSearch.push_back ( new AccessionPreSearch ( searchDatabase, pbsv.second ) );
				else
					preSearch.push_back ( new StandardPreSearch ( params, searchDatabase, getAddAccessionNumbers ( menuDatabase, addAccessionNumbers ), prefix ) );
			}
		}
	}
}
PreSearchInfo::~PreSearchInfo ()
{
	for ( int i = 0 ; i < preSearch.size () ; i++ ) {
		delete preSearch [i];
	}
}
PairBoolStringVector PreSearchInfo::getAccNumSearch ( const string& database, const MapStringToStringVector& accessionNumbers )
{
	bool flag = !accessionNumbers.empty ();
	MapStringToStringVectorConstIterator cur;
	if ( flag ) {
		cur = accessionNumbers.find ( "" );	// database is blank
		if ( cur != accessionNumbers.end () ) flag = true;
		else {
			cur = accessionNumbers.find ( database );
			flag = ( cur != accessionNumbers.end () );
		}
	}
	if ( flag )	return make_pair ( true, (*cur).second );
	else		return make_pair ( false, StringVector () );
}
StringVector PreSearchInfo::getAddAccessionNumbers ( const string& database, const MapStringToStringVector& addAccessionNumbers )
{
	if ( !addAccessionNumbers.empty () ) {
		MapStringToStringVectorConstIterator cur = addAccessionNumbers.find ( "" );
		if ( cur != addAccessionNumbers.end () ) return (*cur).second;
		else {
			cur = addAccessionNumbers.find ( database );
			if ( cur != addAccessionNumbers.end () ) return (*cur).second;
		}
	}
	return StringVector ();
}
void PreSearchInfo::setAccessionNumbers ( MapStringToStringVector& accessionNumbers, const StringVector& aNum )
{
	StringVector sv;
	string d;
	for ( int i = 0 ; i < aNum.size () ; i++ ) {
		string acc = aNum [i];
		if ( acc [0] == '>' ) {
			if ( !sv.empty () ) {
				accessionNumbers [d] = sv;
				sv.clear ();
			}
			d = gen_strtrim ( acc.substr ( 1 ) );
		}
		else {
			if ( isPrefix ( acc, "[-" ) ) {
				acc = acc.substr ( 1, acc.length () - 2 );
			}
			sv.push_back ( acc );
		}
	}
	if ( !sv.empty () ) accessionNumbers [d] = sv;
	for ( MapStringToStringVectorIterator j = accessionNumbers.begin () ; j != accessionNumbers.end () ; j++ ) {
		stable_sort ( (*j).second.begin (), (*j).second.end () );
		(*j).second.erase ( unique ( (*j).second.begin (), (*j).second.end () ), (*j).second.end () );
	}
}
void PreSearchInfo::doSearch ( const vector <FastaServer*>& fs ) const
{
	int prev = -1;
	for ( int i = 0 ; i < fs.size () ; i++ ) {
		int cur = preSearchIndex [i];
		if ( cur != prev ) {
			preSearch [cur]->doSearch ( fs [i] );
			prev = cur;
		}
	}
}
const IntVector& PreSearchInfo::getIndicies ( int num ) const
{
	return preSearch [preSearchIndex[num]]->getIndicies ();
}
int PreSearchInfo::getNumIndicies ( int num ) const
{
	return preSearch [preSearchIndex[num]]->getNumIndicies ();
}
void PreSearchInfo::printHTML ( ostream& os ) const
{
	for ( int i = 0 ; i < preSearch.size () ; i++ ) {
		preSearch [i]->printHTML ( os );
	}
}
void PreSearchInfo::printBodyXML ( ostream& os ) const
{
	for ( int i = 0 ; i < preSearch.size () ; i++ ) {
		preSearch [i]->printXMLResults ( os );
	}
}
void PreSearchInfo::putCGI ( ostream& os ) const
{
	printCGI ( os, RESULTS_FROM_FILE, fromFile );
	preSearch [0]->putCGI ( os );
}
void PreSearchInfo::putHiddenFormEntryJavascript ( ostream& os ) const
{
	printHTMLFORMJavascriptHidden ( os, RESULTS_FROM_FILE, fromFile );
	preSearch [0]->putHiddenFormEntryJavascript ( os );
}
void PreSearchInfo::setReportTaxonomy ()
{
	StandardPreSearch::setReportTaxonomy ( true );
}
string PreSearch::FULL_MW_RANGE = "full_mw_range";
string PreSearch::FULL_PI_RANGE = "full_pi_range";
string PreSearch::SPECIES = "species";
PreSearch::PreSearch ( const string& database, bool concat ) :
	database ( database ),
	concat ( concat )
{
}
PreSearch::~PreSearch () {}
void PreSearch::getIndiciesFromAccessionNumbers ( const string& filename, const StringVector& accessionNumbers, IntVector& indicies )
{
	AccessionNumberMap* am = getAccessionNumberMap ( filename );
	for ( StringVectorSizeType i = 0 ; i < accessionNumbers.size () ; i++ ) {
		const string& acc = accessionNumbers [i];
		int indexNum = am->getIndexNumber ( acc.c_str () );

		if ( indexNum == -1 ) {
			if ( acc.length () > 100 )
				ErrorHandler::genError ()->warning ( "Accession number too long. Each accession number needs to be on a separate line.\n" );
			else
				ErrorHandler::genError ()->warning ( "The accession number " + acc + " does not exist in the database.\n" );
		}
		else
			indicies.push_back ( indexNum );
	}
	delete am;
}
StringVector PreSearch::addConcatAccessionNumbers ( const StringVector& accessionNumbers )
{
	StringVector an = accessionNumbers;
	for ( StringVectorSizeType i = 0 ; i < accessionNumbers.size () ; i++ ) {
		an.push_back ( "-" + accessionNumbers [i] );
	}
	return an;
}
void PreSearch::printXMLResults ( ostream& os ) const
{
	os << "<pre_search_results>" << endl;
	ParameterList::printXML ( os, "num_entries", numEntries );
	ParameterList::printXML ( os, "num_final_indicies", finalIndicies.size () );
	os << "</pre_search_results>" << endl;
}
string PreSearch::getMessage () const
{
	string message = "Pre Search Results";
	message += " (";
	message += database;
	message += ")";
	return message;
}

IndexPreSearch::IndexPreSearch ( const ParameterList* params, const string& database, int dbIdx ) :
	PreSearch ( database ),
	inputFilename	( params->getStringValue	( "input_filename" ) ),
	inputProgramName( params->getStringValue	( "input_program_name", "msfit" ) ),
	inputDirectory	( params->getStringValue	( "results_input_dir", string ( "results" ) + SLASH + string ( inputProgramName ) ) )
{
	ParameterList fParams ( adjustPPOutputPath ( inputDirectory ) + inputFilename + ".res" );
	string dbIdxStr = ( dbIdx == 1 ) ? "" : gen_itoa ( dbIdx );
	string db = fParams.getStringValue ( "database" + dbIdxStr );
	if ( db != database ) {
		ErrorHandler::genError ()->error ( "When searching previously saved results you need to select exactly the same database(s).\n" );
	}
	databaseDate = fParams.getStringValue ( "database_date" + dbIdxStr );
	const char* value;
	if ( fParams.getValue ( "indicies" + dbIdxStr, value ) )	getPostQueryVector ( value, finalIndicies, ' ' );
}
void IndexPreSearch::doSearch ( FastaServer* fs )
{
	numEntries = fs->getNumEntries ();
}
void IndexPreSearch::printHTML ( ostream& os ) const
{
	os << "<p>" << endl;
	ExpandableJavascriptBlock ejb ( getMessage () );
	ejb.printHeader ( os );

	os << "<a href =\"" << ( genIsFullPath ( inputDirectory ) ? "file://" : "../" );
	os << genTranslateSlash ( inputDirectory ) << "/" << inputFilename << ".htm\">Saved Results</a>";
	os << "<br />";

	if ( databaseDate != "" ) {
		if ( databaseDate != getIndexFileLastModifiedTime ( database ) ) {
			os << "<b>WARNING!!! The database index file has been regenerated since the index numbers were stored.<br />";
			os << " i.e. The program FA-Index was re-run so the saved Hits may not be properly retrieved.</b><br />";
		}
	}
	ejb.printFooter ( os );
	os << "</p>" << endl;
}
void IndexPreSearch::putCGI ( ostream& os ) const
{
	printCGIString ( os, "input_filename", inputFilename );
	printCGIString ( os, "input_program_name", inputProgramName );
	printCGIString ( os, "results_input_dir", inputDirectory );
}
void IndexPreSearch::putHiddenFormEntryJavascript ( ostream& os ) const
{
	printHTMLFORMJavascriptHidden ( os, "input_filename", inputFilename );
	printHTMLFORMJavascriptHidden ( os, "input_program_name", inputProgramName );
	printHTMLFORMJavascriptHidden ( os, "results_input_dir", inputDirectory );
}
string AccessionPreSearch::ACCESSION_NUMS = "accession_nums";
AccessionPreSearch::AccessionPreSearch ( const string& database, const StringVector& accessionNumbers ) :
	PreSearch ( database, database.find ( "concat" ) != string::npos ),
	accessionNumbers ( accessionNumbers )
{
}
void AccessionPreSearch::doSearch ( FastaServer* fs )
{
	numEntries = fs->getNumEntries ();
	getIndiciesFromAccessionNumbers ( fs->getFileName (), concat ? addConcatAccessionNumbers ( accessionNumbers ) : accessionNumbers, finalIndicies );
}
void AccessionPreSearch::printHTML ( ostream& os ) const
{
	os << "<p>" << endl;
	ExpandableJavascriptBlock ejb ( getMessage () );
	ejb.printHeader ( os );
	os << "Accession number search selects ";
	os << "<b>" << finalIndicies.size () << "</b>";
	os << ( finalIndicies.size () == 1 ? " entry." : " entries." );
	ejb.printFooter ( os );
	os << "</p>" << endl;
}
void AccessionPreSearch::putCGI ( ostream& os ) const
{
	printCGIContainer ( os, ACCESSION_NUMS, accessionNumbers );
}
void AccessionPreSearch::putHiddenFormEntryJavascript ( ostream& os ) const
{
	printHTMLFORMJavascriptHiddenContainer ( os, ACCESSION_NUMS, accessionNumbers );
}
UserProteinPreSearch::UserProteinPreSearch ( const string& prefix ) :
	PreSearch ( "User Protein" ),
	prefix ( prefix )
{
}
void UserProteinPreSearch::printHTML ( ostream& os ) const
{
	os << "<p>" << endl;
	ExpandableJavascriptBlock ejb ( getMessage () );
	ejb.printHeader ( os );
	os << "User Protein entry contains ";
	os << "<b>" << finalIndicies.size () << "</b>";
	os << ( finalIndicies.size () == 1 ? " entry." : " entries." );
	ejb.printFooter ( os );
	os << "</p>" << endl;
}
void UserProteinPreSearch::doSearch ( FastaServer* fs )
{
	numEntries = fs->getNumEntries ();
	finalIndicies.resize ( numEntries );
	for ( int i = 0 ; i < numEntries ; i++ ) {
		finalIndicies [i] = i + 1;
	}
}
void UserProteinPreSearch::putHiddenFormEntryJavascript ( ostream& os ) const
{
	StringVector sv;
	sv.push_back ( "All" );
	printHTMLFORMJavascriptHiddenContainer2 ( os, SPECIES, sv );
	printHTMLFORMJavascriptHidden ( os, prefix + FULL_MW_RANGE, true );
	printHTMLFORMJavascriptHidden ( os, FULL_PI_RANGE, true );
}
string StandardPreSearch::PROT_LOW_MASS = "prot_low_mass";
string StandardPreSearch::PROT_HIGH_MASS = "prot_high_mass";
string StandardPreSearch::LOW_PI = "low_pi";
string StandardPreSearch::HIGH_PI = "high_pi";
string StandardPreSearch::SPECIES_REMOVE = "species_remove";
string StandardPreSearch::SPECIES_NAMES = "species_names";
string StandardPreSearch::ADD_ACC_NUMS = "add_accession_numbers";
string StandardPreSearch::NAMES = "names";

bool StandardPreSearch::reportTaxonomy = false;
StandardPreSearch::StandardPreSearch ( const ParameterList* params, const string& database, const StringVector& addAccNums, const string& prefix ) :
	PreSearch ( database, database.find ( "concat" ) != string::npos ),
	prefix ( prefix ),
	fullMWRange	( params->getBoolValue	( prefix + FULL_MW_RANGE, false ) ),
	lowMass		( params->getDoubleValue( prefix + PROT_LOW_MASS, 1000.0 ) ),
	lowMassStr	( params->getStringValue( prefix + PROT_LOW_MASS, "1000.0" ) ),
	highMass	( params->getDoubleValue( prefix + PROT_HIGH_MASS, 100000.0 ) ),
	highMassStr	( params->getStringValue( prefix + PROT_HIGH_MASS, "100000.0" ) ),
	fullPIRange	( params->getBoolValue	( FULL_PI_RANGE, false ) ),
	lowPI		( params->getDoubleValue( LOW_PI, 3.0 ) ),
	lowPIStr	( params->getStringValue( LOW_PI, "3.0" ) ),
	highPI		( params->getDoubleValue( HIGH_PI, 10.0 ) ),
	highPIStr	( params->getStringValue( HIGH_PI, "10.0" ) ),
	species		( params->getStringVectorValue ( SPECIES ) ),
	speciesRemove( params->getBoolValue( SPECIES_REMOVE, false ) ),
	names		( params->getPQStringVectorValue ( NAMES ) ),
	speciesNames( params->getPQStringVectorValue ( SPECIES_NAMES ) ),
	addAccNums	( addAccNums )
{
	if ( !speciesNames.empty () )
		fullTaxonomyList = speciesNames;
	else
		fullTaxonomyList = species;

	if ( fullTaxonomyList.empty () ) {
		ErrorHandler::genError ()->error ( "No taxonomy options have been selected.<br />" );
	}
	StringVectorIterator iter = find ( fullTaxonomyList.begin (), fullTaxonomyList.end (), "All" );
	if ( iter != fullTaxonomyList.end () ) {
		if ( fullTaxonomyList.size () == 1 ) {
			fullTaxonomyList.clear ();
		}
		else {
			ErrorHandler::genError ()->error ( "If taxonomy is set to All then it must be the only taxonomy selection.<br />" );
		}
	}
}
void StandardPreSearch::putCGI ( ostream& os ) const
{
	printCGIContainer ( os, SPECIES, species );
	if ( speciesRemove )
		printCGI ( os, SPECIES_REMOVE, speciesRemove );
	if ( fullMWRange )
		printCGI ( os, prefix + FULL_MW_RANGE, fullMWRange );
	else {
		printCGI ( os, prefix + PROT_LOW_MASS, lowMassStr );
		printCGI ( os, prefix + PROT_HIGH_MASS, highMassStr );
	}
	if ( fullPIRange )
		printCGI ( os, FULL_PI_RANGE, fullPIRange );
	else {
		printCGI ( os, LOW_PI, lowPIStr );
		printCGI ( os, HIGH_PI, highPIStr );
	}
	printCGIContainer ( os, NAMES, names );
	printCGIContainer ( os, SPECIES_NAMES, speciesNames );
	printCGIContainer ( os, ADD_ACC_NUMS, addAccNums );
}
void StandardPreSearch::putHiddenFormEntryJavascript ( ostream& os ) const
{
	printHTMLFORMJavascriptHiddenContainer2 ( os, SPECIES, species );
	if ( speciesRemove )
		printHTMLFORMJavascriptHidden ( os, SPECIES_REMOVE, speciesRemove );
	if ( fullMWRange )
		printHTMLFORMJavascriptHidden ( os, prefix + FULL_MW_RANGE, fullMWRange );
	else {
		printHTMLFORMJavascriptHidden ( os, prefix + PROT_LOW_MASS, lowMassStr );
		printHTMLFORMJavascriptHidden ( os, prefix + PROT_HIGH_MASS, highMassStr );
	}
	if ( fullPIRange )
		printHTMLFORMJavascriptHidden ( os, FULL_PI_RANGE, fullPIRange );
	else {
		printHTMLFORMJavascriptHidden ( os, LOW_PI, lowPIStr );
		printHTMLFORMJavascriptHidden ( os, HIGH_PI, highPIStr );
	}
	printHTMLFORMJavascriptHiddenContainer ( os, NAMES, names );
	printHTMLFORMJavascriptHiddenContainer ( os, SPECIES_NAMES, speciesNames );
	printHTMLFORMJavascriptHiddenContainer ( os, ADD_ACC_NUMS, addAccNums );
}
void StandardPreSearch::doSearch ( FastaServer* fs )
{
	string filename = fs->getFileName ();
	numEntries = fs->getNumEntries ();

	UpdatingJavascriptMessage ujm;
	ujm.writeMessage ( cout, "Doing Pre Search." );
	numMWIndicies = numEntries;
	numPIIndicies = numEntries;
	IntVector fullIndicies;
	IntVector mwIndicies;
	IntVector pIIndicies;
	IntVector mwPiIndicies;
	IntVector& indiciesRef = fullIndicies;
	if ( fullPIRange && fullMWRange ) {
		fullIndicies.resize ( numEntries );
		for ( int i = 0 ; i < numEntries ; i++ ) {
			fullIndicies [i] = i + 1;
		}
		indiciesRef = fullIndicies;
	}
	else {
		if ( !fullMWRange ) {
			ujm.writeMessage ( cout, "Doing Protein MW search." );
			if ( fs->getDNADatabase () ) ErrorHandler::genError ()->error ( "The molecular weight pre-search option is not available for DNA databases.<br />" );
			if ( highMass <= lowMass ) {
				ErrorHandler::genError ()->error ( "The protein MW pre-search high mass is less than or equal to the low mass.<br />" );
			}
			floatingPointListSearch ( SeqdbDir::instance ().getSeqdbDir () + filename + ".mw", lowMass, highMass, mwIndicies );
			numMWIndicies = mwIndicies.size ();
			if ( fullPIRange ) indiciesRef = mwIndicies;
		}
		if ( !fullPIRange ) {
			ujm.writeMessage ( cout, "Doing Protein pI search." );
			if ( fs->getDNADatabase () ) ErrorHandler::genError ()->error ( "The pI pre-search option is not available for DNA databases.<br />" );
			if ( highPI <= lowPI ) {
				ErrorHandler::genError ()->error ( "The protein pI pre-search high pI is less than or equal to the low pI.<br />" );
			}
			floatingPointListSearch ( SeqdbDir::instance ().getSeqdbDir () + filename + ".pi", lowPI, highPI, pIIndicies );
			numPIIndicies = pIIndicies.size ();
			if ( fullMWRange ) indiciesRef = pIIndicies;
		}
		if ( !fullMWRange && !fullPIRange ) {
			set_intersection ( mwIndicies.begin (), mwIndicies.end (), pIIndicies.begin (), pIIndicies.end (), back_inserter ( mwPiIndicies ) );
			indiciesRef = mwPiIndicies;
		}
	}
	if ( !fullTaxonomyList.empty () ) {
		ujm.writeMessage ( cout, "Doing Taxonomy search." );
		TaxonomySearch taxSearch ( fs->getFileName (), fs->getNumEntries (), fullTaxonomyList, speciesNames.empty () );
		if ( reportTaxonomy ) {
			IntVector taxNumEntries = taxSearch.getTaxListSize ();
			IntVector taxList = taxSearch.getTaxList ();
			NCBISpeciesMap* nsm = 0;
			for ( IntVectorSizeType i = 0 ; i < taxList.size () ; i++ ) {
				string scientificName;
				int t = taxList [i];
				if ( t == -1 ) scientificName = "Unreadable";
				else if ( t == 0 ) scientificName = "Unknown";
				else {
					if ( nsm == 0 ) nsm = new NCBISpeciesMap ( true );
					scientificName = nsm->getScientificName ( t );
				}
				taxonomyHits.push_back ( TaxonomyHit ( scientificName, t, taxNumEntries [i] ) );
			}
			delete nsm;
			stable_sort ( taxonomyHits.begin (), taxonomyHits.end (), SortTaxonomyAscending () );
		}
		const IntVector& taxonomyIndicies = taxSearch.getIndicies ();
		numTaxIndicies = taxonomyIndicies.size ();
		if ( speciesRemove )
			set_difference ( indiciesRef.begin (), indiciesRef.end (), taxonomyIndicies.begin (), taxonomyIndicies.end (), back_inserter ( finalIndicies ) );
		else
			set_intersection ( indiciesRef.begin (), indiciesRef.end (), taxonomyIndicies.begin (), taxonomyIndicies.end (), back_inserter ( finalIndicies ) );
	}
	else
		finalIndicies = indiciesRef;

	numPreNameIndicies = finalIndicies.size ();
	if ( names.size () ) {
		ujm.writeMessage ( cout, "Doing Name search." );
		name_regexp_pre_search ( fs, finalIndicies, names );
	}

	if ( !addAccNums.empty () ) {
		ujm.writeMessage ( cout, "Inserting additional accession numbers." );
		IntVector addAccNumberIndicies;
		getIndiciesFromAccessionNumbers ( filename, concat ? addConcatAccessionNumbers ( addAccNums ) : addAccNums, addAccNumberIndicies );

		finalIndicies.insert ( finalIndicies.end (), addAccNumberIndicies.begin (), addAccNumberIndicies.end () );
		sort ( finalIndicies.begin (), finalIndicies.end () );
	}
	ujm.deletePreviousMessage ( cout );
}
void StandardPreSearch::floatingPointListSearch ( const string& filename, double low, double high, IntVector& serialNumber )
{
	GenIFStream inFile ( filename, ios::binary );
	for ( ; ; ) {
		IndexedDouble id;
		inFile.read ( (char*) &id, sizeof (IndexedDouble) );
		if ( inFile.eof () ) break;
		double number = id.number;
		if ( number >= low ) {
			if ( number > high ) break;
			serialNumber.push_back ( id.index );
		}
	}
	sort ( serialNumber.begin (), serialNumber.end () );
}
void StandardPreSearch::name_regexp_pre_search ( FastaServer* fs, IntVector& indicies, const StringVector& names )
{
	int numNames = names.size ();
	vector <RegularExpression*> regExp;
	int i, j, k;

	for ( i = 0 ; i < numNames ; i++ ) {
		regExp.push_back ( new RegularExpression ( names [i].c_str (), true ) );
	}
	int numIndicies = indicies.size ();
	for ( i = 0, k = 0 ; i < numIndicies ; i++ ) {
		for ( fs->firstLine ( indicies [i] ) ; fs->isDoneLine () ; fs->nextLine () ) {
			string str = fs->getLineName ();

			for ( j = 0 ; j < numNames ; j++ ) {
				if ( regExp [j]->isPresent ( str.c_str () ) ) {
					indicies [k++] = indicies [i];
					goto label;
				}
			}
		}
		label:;
	}
	indicies.resize ( k );
	for ( i = 0 ; i < numNames ; i++ ) {
		delete regExp [i];
	}
}
void StandardPreSearch::printXMLResults ( ostream& os ) const
{
	os << "<pre_search_results>" << endl;
	ParameterList::printXML ( os, "num_entries", numEntries );
	if ( !fullMWRange ) ParameterList::printXML ( os, "num_mw_indicies", numMWIndicies );
	if ( !fullPIRange ) ParameterList::printXML ( os, "num_pi_indicies", numPIIndicies );
	if ( !fullTaxonomyList.empty () ) ParameterList::printXML ( os, "num_taxonomy_indicies", numTaxIndicies );

	if ( !fullMWRange || !fullPIRange ) {
		if ( !fullTaxonomyList.empty () )
			ParameterList::printXML ( os, "num_mw_pi_taxonomy_indicies", numPreNameIndicies );
		else
			ParameterList::printXML ( os, "num_mw_pi_indicies", numPreNameIndicies );
	}
	ParameterList::printXML ( os, "num_final_indicies", finalIndicies.size () );
	os << "</pre_search_results>" << endl;
}
void StandardPreSearch::printHTML ( ostream& os ) const
{
	os << "<p>" << endl;
	ExpandableJavascriptBlock ejb ( getMessage (), finalIndicies.empty () );
	ejb.printHeader ( os );
	os << "Number of entries in the database: <b>" << numEntries << "</b><br />";

	if ( fullMWRange )	os << "Full Molecular Weight range: ";
	else				os << "Molecular weight search (<b>" << (int)lowMass << " - " << (int)highMass << " Da</b>) selects ";
	os << "<b>" << numMWIndicies << "</b> " << ( numMWIndicies == 1 ? "entry." : "entries." ) << "<br />" << endl;

	if ( fullPIRange ) os << "Full pI range: <b>" << numPIIndicies << "</b>";
	else {
		os << "pI search (<b>";
		genPrint ( os, lowPI, 2 );
		os << " - ";
		os << highPI;
		os << "</b>) selects <b>";
		os << numPIIndicies;
		os << "</b>";
	}
	os << ( numPIIndicies == 1 ? " entry." : " entries." ) << "<br />" << endl;
	if ( !fullTaxonomyList.empty () ) {
		os << "Taxonomy search <b> ";
		for ( StringVectorSizeType i = 0 ; i < fullTaxonomyList.size () ; i++ ) {
			os << fullTaxonomyList [i];
			if ( i != fullTaxonomyList.size () - 1 ) os << ", ";
		}
		os << " </b>";
		os << " selects <b>" << numTaxIndicies << "</b>";
		os << ( numTaxIndicies == 1 ? " entry." : " entries." );
		if ( speciesRemove ) os << " Removed.";
		os << "<br />";
		os << endl;
	}
	if ( !fullMWRange || !fullPIRange ) {
		os << "Combined molecular weight";
		if ( fullTaxonomyList.empty () )
			os << " and pI";
		else
			os << ", pI and taxonomy";
		os << " searches select <b>" << numPreNameIndicies << "</b>";
		os << ( numPreNameIndicies == 1 ? " entry." : " entries." ) << "<br />" << endl;
	}
	if ( finalIndicies.empty () ) os << "<font color=\"#FF0000\" size=\"+2\">" << endl;
	os << "Pre searches select <b>" << finalIndicies.size () << "</b>";
	os << ( finalIndicies.size () == 1 ? " entry." : " entries." ) << "<br />" << endl;
	if ( finalIndicies.empty () ) os << "</font>" << endl;
	if ( !addAccNums.empty () ) {
		int numExtra = concat ? addAccNums.size () * 2 : addAccNums.size ();
		os << "Includes <b>" << numExtra << "</b> additional";
		os << ( numExtra == 1 ? " entry." : " entries." ) << "<br />" << endl;
	}
	if ( reportTaxonomy ) outputHTMLSpeciesSearchResults ( os );
	os << endl;
	ejb.printFooter ( os );
	os << "</p>" << endl;
}
void StandardPreSearch::outputHTMLSpeciesSearchResults ( ostream& os ) const
{
	if ( taxonomyHits.size () ) {
		tableStart ( os, true );
			tableRowStart ( os );
				tableHeader ( os, "Species String" );
				tableHeader ( os, "Taxonomy Node" );
				tableHeader ( os, "Number of Entries" );
			tableRowEnd ( os );
			for ( StringVectorSizeType i = 0 ; i < taxonomyHits.size () ; i++ ) {
				tableRowStart ( os );
					tableDataStart ( os, "", "left", true );
						os << taxonomyHits [i].scientificName << endl;
					tableDataEnd ( os );
					tableDataStart ( os, "", "right", true );
						os << taxonomyHits [i].node << endl;
					tableDataEnd ( os );
					tableDataStart ( os, "", "right", true );
						os << taxonomyHits [i].numEntries << endl;
					tableDataEnd ( os );
				tableRowEnd ( os );
			}
		tableEnd ( os );
		os << "<br />" << endl;
	}
/*
	for ( StringVectorSizeType i = 0 ; i < taxonomyHits.size () ; i++ ) {
		os << taxonomyHits [i].scientificName;
		//os << taxonomyHits [i].scientificName << '\t';
		//os << taxonomyHits [i].node << '\t';
		//os << taxonomyHits [i].numEntries;
		os << "<br />" << endl;
	}
*/
}
