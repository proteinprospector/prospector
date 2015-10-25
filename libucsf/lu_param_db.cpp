/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_param_db.cpp                                               *
*                                                                             *
*  Created    : Dec 18th 2014                                                 *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2014-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <string>
#include <sqlite3.h>
#include <lg_string.h>
#include <lgen_file.h>
#include <lu_aa_info.h>
#include <lu_dig_par.h>
#include <lu_disc_sc.h>
#include <lu_fas_enz.h>
#include <lu_getfil.h>
#include <lu_mass_elem.h>
#include <lu_mut_mtrx.h>
#include <lu_param_db.h>
#include <lu_param_list.h>
#ifdef RAW_DATA
#include <lu_quan_multi.h>
#endif
#include <lu_repos_info.h>
#include <lu_usermod.h>
#include <lu_viewer_form.h>

using std::string;
using std::istringstream;
using std::getline;
using std::make_pair;
using std::back_inserter;
using std::copy;
using std::vector;

class SQLiteTable {
	static ParamDBWrite* pdbw;
protected:
	StringVector n1;
	StringVector v1;
	StringVector n2;
	StringVector v2;
	StringVector n3;
	StringVector v3;
public:
	SQLiteTable () {}
	virtual void create () = 0;
	virtual void insert () = 0;
	virtual void read ();
	void createTables ( const StringVector& sql );
	void insertTable ( const string& tableName, const StringVector& names, const StringVector& values );
	static void setDatabaseConnection ( ParamDBWrite* p );
};
ParamDBWrite* SQLiteTable::pdbw = 0;
void SQLiteTable::read ()
{
	// default maybe select * from table
}
void SQLiteTable::setDatabaseConnection ( ParamDBWrite* p )
{
	pdbw = p;
}
void SQLiteTable::createTables ( const StringVector& sql )
{
	for ( StringVectorSizeType i = 0 ; i < sql.size () ; i++ ) {
		pdbw->executeSQL ( sql [i] );
	}
}
void SQLiteTable::insertTable ( const string& tableName, const StringVector& names, const StringVector& values )
{
	pdbw->insertQuery ( tableName, names, values );
}
class SQLiteList : public SQLiteTable {
	string fileName;
	string tableName;
public:
	SQLiteList ( const string& tableName, const string& fileName );
	void create ();
	virtual void insert ();
};
SQLiteList::SQLiteList ( const string& tableName, const string& fileName ) :
	tableName ( tableName ),
	fileName ( fileName )
{
	n1.push_back ( "name" );
	n1.push_back ( "distribution" );
	n1.push_back ( "deactivate" );

	v1.resize ( n1.size () );
}
void SQLiteList::create ()
{
	StringVector sql;
	sql.push_back ( "CREATE TABLE " + tableName + " ( "\
			+ n1 [0] + " TEXT PRIMARY KEY NOT NULL, "\
			+ n1 [1] + " INTEGER CHECK ( " + n1 [1] + " = 0 OR " + n1 [1] + " = 1 ), "\
			+ n1 [2] + " INTEGER CHECK ( " + n1 [2] + " = 0 OR " + n1 [2] + " = 1 ) "\
		+ ")" );
	createTables ( sql );
}
void SQLiteList::insert ()
{
	v1 [1] = "1";
	v1 [2] = "0";
	string path = MsparamsDir::instance ().getParamPath ( fileName );
	if ( genFileExists ( path ) ) {
		GenCommentedIFStream ifs ( path );
		string line;
		while ( ifs.getUncommentedLine ( line ) ) {
			v1 [0] = genEscapeQuote ( line );
			insertTable ( tableName, n1, v1 );
		}
	}
}

class SQLiteNameValueList : public SQLiteTable {
protected:
	string tableName;
	string fileName;
public:
	SQLiteNameValueList ( const string& tableName, const string& fileName );
	virtual void create ();
	virtual void insert ();
};
SQLiteNameValueList::SQLiteNameValueList ( const string& tableName, const string& fileName ) :
	tableName ( tableName ),
	fileName ( fileName )
{
	n1.push_back ( "name" );
	n1.push_back ( "value" );
	n1.push_back ( "defaultValue" );
	n1.push_back ( "distribution" );
	n1.push_back ( "deactivate" );

	v1.resize ( n1.size () );
}
void SQLiteNameValueList::create ()
{
	StringVector sql;
	sql.push_back ( "CREATE TABLE " + tableName + " ( "\
			+ n1 [0] + " TEXT PRIMARY KEY NOT NULL, "\
			+ n1 [1] + " TEXT NOT NULL, "\
			+ n1 [2] + " TEXT, "\
			+ n1 [3] + " INTEGER CHECK ( " + n1 [3] + " = 0 OR " + n1 [3] + " = 1 ), "\
			+ n1 [4] + " INTEGER CHECK ( " + n1 [4] + " = 0 OR " + n1 [4] + " = 1 ) "\
		+ ")" );
	createTables ( sql );
}
void SQLiteNameValueList::insert ()
{
	v1 [2] = "";
	v1 [3] = "1";
	v1 [4] = "0";
	string path = MsparamsDir::instance ().getParamPath ( fileName );
	if ( genFileExists ( path ) ) {
		GenNameValueStream nvs ( path );
		StringVector sv = nvs.getNameList ();
		for ( int i = 0 ; i < sv.size () ; i++ ) {
			v1 [0] = sv [i];
			v1 [1] = nvs.getStringValue ( sv [i] );
			insertTable ( tableName, n1, v1 );
		}
	}
}
class SQLiteNonUniqueNameValueList : public SQLiteNameValueList {
public:
	SQLiteNonUniqueNameValueList ( const string& tableName, const string& fileName );
	virtual void create ();
	virtual void insert ();
};
SQLiteNonUniqueNameValueList::SQLiteNonUniqueNameValueList ( const string& tableName, const string& fileName ) :
	SQLiteNameValueList ( tableName, fileName )
{
}
void SQLiteNonUniqueNameValueList::create ()
{
	StringVector sql;
	sql.push_back ( "CREATE TABLE " + tableName + " ( "\
			+ n1 [0] + " TEXT NOT NULL, "\
			+ n1 [1] + " TEXT NOT NULL, "\
			+ n1 [2] + " TEXT, "\
			+ n1 [3] + " INTEGER CHECK ( " + n1 [3] + " = 0 OR " + n1 [3] + " = 1 ), "\
			+ n1 [4] + " INTEGER CHECK ( " + n1 [4] + " = 0 OR " + n1 [4] + " = 1 ), "\
			+ "PRIMARY KEY ( "\
				+ n1 [0] + ", "\
				+ n1 [1]\
			+ " ) "\
		+ ")" );
	createTables ( sql );
}
void SQLiteNonUniqueNameValueList::insert ()
{
	v1 [2] = "";
	v1 [3] = "1";
	v1 [4] = "0";
	string path = MsparamsDir::instance ().getParamPath ( fileName );
	if ( genFileExists ( path ) ) {
		GenNameValueStream nvs ( path );
		StringVector sv = nvs.getNameList ();
		for ( StringVectorSizeType i = 0 ; i < sv.size () ; i++ ) {
			v1 [0] = sv [i];
			StringVector sv2 = nvs.getStringVectorValue ( sv [i] );
			for ( StringVectorSizeType j = 0 ; j < sv2.size () ; j++ ) {
				v1 [1] = sv2 [j];
				insertTable ( tableName, n1, v1 );
			}
		}
	}
}
class SQLiteIDNameValueList : public SQLiteTable {
protected:
	string tableName;
	string fileName;
public:
	SQLiteIDNameValueList ( const string& tableName, const string& fileName );
	void create ();
	virtual void insert ();
};
SQLiteIDNameValueList::SQLiteIDNameValueList ( const string& tableName, const string& fileName ) :
	tableName ( tableName ),
	fileName ( fileName )
{
	n1.push_back ( "id" );
	n1.push_back ( "name" );
	n1.push_back ( "value" );
	n1.push_back ( "defaultValue" );
	n1.push_back ( "distribution" );
	n1.push_back ( "deactivate" );

	v1.resize ( n1.size () );
}
void SQLiteIDNameValueList::create ()
{
	StringVector sql;
	sql.push_back ( "CREATE TABLE " + tableName + " ( "\
			+ n1 [0] + " TEXT NOT NULL, "\
			+ n1 [1] + " TEXT NOT NULL, "\
			+ n1 [2] + " TEXT, "\
			+ n1 [3] + " TEXT, "\
			+ n1 [4] + " INTEGER CHECK ( " + n1 [4] + " = 0 OR " + n1 [4] + " = 1 ), "\
			+ n1 [5] + " INTEGER CHECK ( " + n1 [5] + " = 0 OR " + n1 [5] + " = 1 ), "\
			+ "PRIMARY KEY ( "\
				+ n1 [0] + ", "\
				+ n1 [1] + ", "\
				+ n1 [2]\
			+ " ) "\
		+ ")" );
	createTables ( sql );
}
void SQLiteIDNameValueList::insert ()
{
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	string path = MsparamsDir::instance ().getParamPath ( fileName );
	if ( genFileExists ( path ) ) {
		GenCommentedIFStream ifs ( path );
		while ( ifs.getUncommentedLine ( v1 [0] ) ) {
			ParameterList p ( ifs );
			StringVector sv = p.getNameList ();
			int numValues = p.size ();
			for ( int i = 0 ; i < numValues ; i++ ) {
				StringVector sv2 = p.getStringVectorValue ( sv [i] );
				for ( int j = 0 ; j < sv2.size () ; j++ ) {
					v1 [1] = sv [i];
					v1 [2] = sv2 [j];
					insertTable ( tableName, n1, v1 );
				}
			}
		}
	}
}
class SQLiteIDUnsortedNameValueList : public SQLiteIDNameValueList {
public:
	SQLiteIDUnsortedNameValueList ( const string& tableName, const string& fileName );
	void insert ();
};
SQLiteIDUnsortedNameValueList::SQLiteIDUnsortedNameValueList ( const string& tableName, const string& fileName ) :
	SQLiteIDNameValueList ( tableName, fileName )
{
}
void SQLiteIDUnsortedNameValueList::insert ()
{
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( fileName ) );
	while ( ifs.getUncommentedLine ( v1 [0] ) ) {
		string line;
		while ( ifs.getUncommentedLine ( line ) ) {
			if ( line.length () != 0 ) {
				if ( line [0] == '>' ) break;
				istringstream ist ( line );
				ist >> v1 [1];
				string temp;
				getline ( ist, temp );
				v1 [2] = temp.substr ( temp.find_first_not_of ( " \t" ) );	// strip leading space
				insertTable ( tableName, n1, v1 );
			}
		}
	}
}
class SQLiteIDSubIDNameValueList : public SQLiteTable {
protected:
	string tableName;
	string fileName;
public:
	SQLiteIDSubIDNameValueList ( const string& tableName, const string& fileName );
	void create ();
	virtual void insert () = 0;
};
SQLiteIDSubIDNameValueList::SQLiteIDSubIDNameValueList ( const string& tableName, const string& fileName ) :
	tableName ( tableName ),
	fileName ( fileName )
{
	n1.push_back ( "id" );
	n1.push_back ( "sub_id" );
	n1.push_back ( "name" );
	n1.push_back ( "value" );
	n1.push_back ( "defaultValue" );
	n1.push_back ( "distribution" );
	n1.push_back ( "deactivate" );

	v1.resize ( n1.size () );
}
void SQLiteIDSubIDNameValueList::create ()
{
	StringVector sql;
	sql.push_back ( "CREATE TABLE " + tableName + " ( "\
			+ n1 [0] + " TEXT NOT NULL, "\
			+ n1 [1] + " TEXT NOT NULL, "\
			+ n1 [2] + " TEXT NOT NULL, "\
			+ n1 [3] + " TEXT, "\
			+ n1 [4] + " TEXT, "\
			+ n1 [5] + " INTEGER CHECK ( " + n1 [5] + " = 0 OR " + n1 [5] + " = 1 ), "\
			+ n1 [6] + " INTEGER CHECK ( " + n1 [6] + " = 0 OR " + n1 [6] + " = 1 ), "\
			+ "PRIMARY KEY ( "\
				+ n1 [0] + ", "\
				+ n1 [1] + ", "\
				+ n1 [2] + ", "\
				+ n1 [3]\
			+ " ) "\
		+ ")" );
	createTables ( sql );
}
class SQLiteIDSubIDSubIDNameValueList : public SQLiteTable {
protected:
	string tableName;
	string fileName;
public:
	SQLiteIDSubIDSubIDNameValueList ( const string& tableName, const string& fileName );
	void create ();
	virtual void insert () = 0;
};
SQLiteIDSubIDSubIDNameValueList::SQLiteIDSubIDSubIDNameValueList ( const string& tableName, const string& fileName ) :
	tableName ( tableName ),
	fileName ( fileName )
{
	n1.push_back ( "id" );
	n1.push_back ( "sub_id" );
	n1.push_back ( "sub_id_2" );
	n1.push_back ( "name" );
	n1.push_back ( "value" );
	n1.push_back ( "defaultValue" );
	n1.push_back ( "distribution" );
	n1.push_back ( "deactivate" );

	v1.resize ( n1.size () );
}
void SQLiteIDSubIDSubIDNameValueList::create ()
{
	StringVector sql;
	sql.push_back ( "CREATE TABLE " + tableName + " ( "\
			+ n1 [0] + " TEXT NOT NULL, "\
			+ n1 [1] + " TEXT NOT NULL, "\
			+ n1 [2] + " TEXT NOT NULL, "\
			+ n1 [3] + " TEXT NOT NULL, "\
			+ n1 [4] + " TEXT, "\
			+ n1 [5] + " TEXT, "\
			+ n1 [6] + " INTEGER CHECK ( " + n1 [6] + " = 0 OR " + n1 [6] + " = 1 ), "\
			+ n1 [7] + " INTEGER CHECK ( " + n1 [7] + " = 0 OR " + n1 [7] + " = 1 ), "\
			+ "PRIMARY KEY ( "\
				+ n1 [0] + ", "\
				+ n1 [1] + ", "\
				+ n1 [2] + ", "\
				+ n1 [3] + ", "\
				+ n1 [4]\
			+ " ) "\
		+ ")" );
	createTables ( sql );
}
class SQLiteParameterLists : public SQLiteIDSubIDNameValueList {
public:
	SQLiteParameterLists ();
	SQLiteParameterLists ( const string& fileName );
	virtual void insert ();
};
SQLiteParameterLists::SQLiteParameterLists () :
	SQLiteIDSubIDNameValueList ( "ParameterLists", "" )
{
}
SQLiteParameterLists::SQLiteParameterLists ( const string& fileName ) :
	SQLiteIDSubIDNameValueList ( "ParameterLists", fileName )
{
}
void SQLiteParameterLists::insert ()
{
	StringVector fileList;
	fileList.push_back ( "batchtag" );
	fileList.push_back ( "dbstat" );
	fileList.push_back ( "msbridge" );
	fileList.push_back ( "mscomp" );
	fileList.push_back ( "msdigest" );
	fileList.push_back ( "msfilter" );
	fileList.push_back ( "msfit" );
	fileList.push_back ( "msfitupload" );
	fileList.push_back ( "mshomology" );
	fileList.push_back ( "msisotope" );
	fileList.push_back ( "msnonspecific" );
	fileList.push_back ( "mspattern" );
	fileList.push_back ( "msproduct" );
	fileList.push_back ( "msseq" );
	fileList.push_back ( "mstag" );
	fileList.push_back ( "msviewer" );
	fileList.push_back ( "searchCompare" );
	v1 [4] = "";
	v1 [5] = "1";
	v1 [6] = "0";
	for ( StringVectorSizeType i = 0 ; i < fileList.size () ; i++ ) {
		ParameterList pList ( MsparamsDir::instance ().getParamPath ( fileList [i] + SLASH + "default.xml" ), false, false, false, false );
		StringVector pNames;
		StringVectorVector pValues;
		pList.getParams ( pNames, pValues );
		v1 [0] = fileList [i];
		v1 [1] = "default";
		for ( StringVectorSizeType j = 0 ; j < pNames.size () ; j++ ) {
			v1 [2] = pNames [j];
			for ( StringVectorSizeType k = 0 ; k < pValues [j].size () ; k++ ) {
				v1 [3] = pValues [j] [k];
				insertTable ( tableName, n1, v1 );
			}
		}
	}
}

class SQLiteAATable : public SQLiteIDNameValueList {
public:
	SQLiteAATable ();
	void insert ();
};
SQLiteAATable::SQLiteAATable () :
	SQLiteIDNameValueList ( "AA", "" )
{
}
void SQLiteAATable::insert ()
{
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	CharVector cv = AAInfo::getInfo ().getAminoAcids ();
	for ( StringVectorSizeType i = 0 ; i < cv.size () ; i++ ) {
		char c = cv [i];
		v1 [0] = AAInfo::getInfo ().getLongName ( c );

		v1 [1] = "one_letter_code";
		v1 [2] = string ( 1, c );
		insertTable ( tableName, n1, v1 );

		v1 [1] = "formula";
		v1 [2] = AAInfo::getInfo ().getElementalString ( c );
		insertTable ( tableName, n1, v1 );

		v1 [1] = "dw_substituent_a";
		string subA = AAInfo::getInfo ().getSubstituentA ( c );
		if ( subA != "0" ) {
			v1 [2] = subA;
			insertTable ( tableName, n1, v1 );
		}

		v1 [1] = "dw_substituent_b";
		string subB = AAInfo::getInfo ().getSubstituentB ( c );
		if ( subB != "0" ) {
			v1 [2] = subB;
			insertTable ( tableName, n1, v1 );
		}

		v1 [1] = "pk_c_term";
		string ctpk = AAInfo::convert_pk ( AAInfo::getInfo ().getCTermPK ( c ) );
		if ( ctpk != "n/a" ) {
			v1 [2] = ctpk;
			insertTable ( tableName, n1, v1 );
		}

		v1 [1] = "pk_n_term";
		string ntpk = AAInfo::convert_pk ( AAInfo::getInfo ().getNTermPK ( c ) );
		if ( ntpk != "n/a" ) {
			v1 [2] = ntpk;
			insertTable ( tableName, n1, v1 );
		}

		v1 [1] = "pk_acidic_sc";
		string ascpk = AAInfo::convert_pk ( AAInfo::getInfo ().getAcidicSCPK ( c ) );
		if ( ascpk != "n/a" ) {
			v1 [2] = ascpk;
			insertTable ( tableName, n1, v1 );
		}

		v1 [1] = "pk_basic_sc";
		string bscpk = AAInfo::convert_pk ( AAInfo::getInfo ().getBasicSCPK ( c ) );
		if ( bscpk != "n/a" ) {
			v1 [2] =  bscpk;
			insertTable ( tableName, n1, v1 );
		}
	}
}
class SQLiteApplets : public SQLiteIDNameValueList {
public:
	SQLiteApplets ();
	void insert ();
};
SQLiteApplets::SQLiteApplets () :
	SQLiteIDNameValueList ( "Applets", "" )
{
}
void SQLiteApplets::insert ()
{
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	StringVector fileList;
	fileList.push_back ( "dbstat_hist" );
	fileList.push_back ( "error_hist" );
	fileList.push_back ( "fit_graph" );
	fileList.push_back ( "hist" );
	fileList.push_back ( "mmod_hist" );
	fileList.push_back ( "pr_graph" );
	fileList.push_back ( "sp_graph" );
	for ( StringVectorSizeType i = 0 ; i < fileList.size () ; i++ ) {
		string path = MsparamsDir::instance ().getParamPath ( fileList [i] + ".par.txt" );
		GenNameValueStream nvs ( path );
		if ( genFileExists ( path ) ) {
			v1 [0] = fileList [i];
			GenCommentedIFStream ifs ( path );
			ParameterList p ( ifs );
			StringVector sv = p.getNameList ();
			int numValues = p.size ();
			for ( int j = 0 ; j < numValues ; j++ ) {
				StringVector sv2 = p.getStringVectorValue ( sv [j] );
				for ( int k = 0 ; k < sv2.size () ; k++ ) {
					v1 [1] = sv [j];
					v1 [2] = sv2 [k];
					insertTable ( tableName, n1, v1 );
				}
			}
		}
	}
}
class SQLiteDBHosts : public SQLiteIDNameValueList {
public:
	SQLiteDBHosts ();
	void insert ();
};
SQLiteDBHosts::SQLiteDBHosts () :
	SQLiteIDNameValueList ( "DBHosts", "dbhosts.txt" )
{
}
void SQLiteDBHosts::insert ()
{
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( fileName ) );
	string line;
	int phase = 0;
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	while ( ifs.getUncommentedLine ( line ) ) {
		if ( phase == 0 ) v1 [0] = line;
		if ( phase == 1 ) v1 [1] = "ftp_url";
		if ( phase == 1 ) {
			if ( isPrefix ( line, "ftp" ) ) {
				v1 [2] = line;
				insertTable ( tableName, n1, v1 );
			}
			else
				phase++;							// This section finished used the line below
		}
		if ( phase == 2 ) v1 [1] = "username";
		if ( phase == 3 ) v1 [1] = "password";
		if ( phase == 4 ) v1 [1] = "compression_ratio";
		if ( phase == 5 ) v1 [1] = "decoy_instruction";
		if ( phase != 1 ) {
			if ( phase != 0 && phase != 6 ) {		// If phase is 6 ignore the line as unused
				v1 [2] = line;
				insertTable ( tableName, n1, v1 );
			}
			phase++;
		}
		if ( phase == 7 ) phase = 0;
	}
}

class SQLiteElements : public SQLiteIDNameValueList {
public:
	SQLiteElements ();
	void insert ();
};
SQLiteElements::SQLiteElements () :
	SQLiteIDNameValueList ( "Elements", "elements.txt" )
{
}
void SQLiteElements::insert ()
{
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( fileName ) );
	string line;
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	while ( ifs.getUncommentedLine ( line ) ) {
		istringstream ist ( line );
		ist >> v1 [0];

		v1 [1] = "valency";
		ist >> v1 [2];
		insertTable ( tableName, n1, v1 );

		int numIsotopes;
		ist >> numIsotopes;
		for ( int i = 1 ; i <= numIsotopes ; i++ ) {
			v1 [1] = "isotope_mass_" + gen_itoa ( i );
			ist >> v1 [2];
			insertTable ( tableName, n1, v1 );
			v1 [1] = "isotope_abundance_" + gen_itoa ( i );
			ist >> v1 [2];
			insertTable ( tableName, n1, v1 );
		}
	}
}
class SQLiteEnzymes : public SQLiteIDNameValueList {
public:
	SQLiteEnzymes ();
	void insert ();
};
SQLiteEnzymes::SQLiteEnzymes () :
	SQLiteIDNameValueList ( "Enzymes", "" )
{
}
void SQLiteEnzymes::insert ()
{
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	StringVector enzymeList = DigestTable::instance ().getNames ();
	for ( string::size_type i = 0 ; i < enzymeList.size () ; i++ ) {
		string eList = enzymeList [i];
		if ( eList.find ( '/' ) == string::npos ) {
			v1 [0] = eList;

			v1 [1] =  "specificity";
			v1 [2] = DigestTable::instance ().getBreakMask ( eList );
			insertTable ( tableName, n1, v1 );

			v1 [1] =  "exceptions";
			v1 [2] = DigestTable::instance ().getExcludeMask ( eList );
			insertTable ( tableName, n1, v1 );

			v1 [1] =  "terminus";
			v1 [2] = DigestTable::instance ().getSpecificity ( eList ) == 'C' ? "C" : "N";
			insertTable ( tableName, n1, v1 );
		}
	}
}
class SQLiteEnzymeCombinations : public SQLiteList {
public:
	SQLiteEnzymeCombinations ();
	void insert ();
};
SQLiteEnzymeCombinations::SQLiteEnzymeCombinations () :
	SQLiteList ( "EnzymeCombinations", "" )
{
}
void SQLiteEnzymeCombinations::insert ()
{
	//################################### There was a problem with the distribution file that the entry CNBr/V8 E appears twice
	v1 [1] = "1";
	v1 [2] = "0";
	StringVector enzymeList = DigestTable::instance ().getNames ();
	for ( string::size_type i = 0 ; i < enzymeList.size () ; i++ ) {
		string eList = enzymeList [i];
		if ( eList.find ( '/' ) != string::npos ) {
			v1 [0] = eList;
			insertTable ( "EnzymeCombinations", n1, v1 );
		}
	}
}
class SQLiteExpectation : public SQLiteParameterLists {
public:
	SQLiteExpectation ();
	void create () {}			// Don't create a new table
	void insert ();
};
SQLiteExpectation::SQLiteExpectation () :
	SQLiteParameterLists ( "expectation.xml" )
{
}
void SQLiteExpectation::insert ()
{
	char* info = getFileAsCharPtr ( MsparamsDir::instance ().getParamPath ( fileName ) );

	string xparams = XMLParser::getStringValue ( info, "parameters" );
	ParameterList expParamList ( xparams, false, false, false );
	StringVector pNames;
	StringVectorVector pValues;
	expParamList.getParams ( pNames, pValues );
	v1 [4] = "";
	v1 [5] = "1";
	v1 [6] = "0";
	v1 [0] = "expectation";
	v1 [1] = "default";
	for ( StringVectorSizeType i = 0 ; i < pNames.size () ; i++ ) {
		v1 [2] = pNames [i];
		for ( StringVectorSizeType j = 0 ; j < pValues [i].size () ; j++ ) {
			v1 [3] = pValues [i] [j];
			insertTable ( tableName, n1, v1 );
		}
	}
	StringVector cparams = XMLParser::getStringVectorValue ( info, "copy_parameter" );
	v1 [0] = "expectation";
	v1 [1] = "copy";
	v1 [3] = "";
	for ( StringVectorSizeType k = 0 ; k < pNames.size () ; k++ ) {
		v1 [2] = pNames [k];
		insertTable ( tableName, n1, v1 );
	}

	delete [] info;
}
class SQLiteIdxLinks : public SQLiteNameValueList {
public:
	SQLiteIdxLinks ();
	void insert ();
};
SQLiteIdxLinks::SQLiteIdxLinks () :
	SQLiteNameValueList ( "IdxLinks", "" )
{
}
void SQLiteIdxLinks::insert ()
{
	v1 [2] = "";
	v1 [3] = "1";
	v1 [4] = "0";
	MSDigestLinkNameValueStream nvs ( false );
	StringVector sv = nvs.getNameList ();
	for ( int i = 0 ; i < sv.size () ; i++ ) {
		v1 [0] = sv [i];
		v1 [1] = nvs.getStringValue ( sv [i] );
		insertTable ( tableName, n1, v1 );
	}
}
class SQLiteImmonium : public SQLiteIDNameValueList {
public:
	SQLiteImmonium ();
	void insert ();
};
SQLiteImmonium::SQLiteImmonium () :
	SQLiteIDNameValueList ( "ImmoniumParams", "imm.txt" )
{
}
void SQLiteImmonium::insert ()
{
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( fileName ) );
	string line;
	StringVector sv;
	while ( ifs.getUncommentedLine ( line ) ) {
		sv.push_back ( line );
	}
	v1 [0] = "params";
	v1 [1] = "immonium_tolerance";
	v1 [2] = sv [0];
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	insertTable ( tableName, n1, v1 );

	v1 [1] = "min_fragment_ion_mass";
	v1 [2] = sv [1];
	insertTable ( tableName, n1, v1 );

	for ( int i = 2 ; i < sv.size () ; i++ ) {
		string str = sv [i];
		string::size_type start = 0;
		string::size_type end = 0;

		string formula = genNextString ( str, "|", start, end );
		string aminoAcids = genNextString ( str, "|", start, end );
		v1 [0] = aminoAcids + gen_ftoa ( formula_to_nominal_mass ( formula.c_str () ), "%.0f" );

		v1 [1] = "formula";
		v1 [2] = formula;
		insertTable ( tableName, n1, v1 );

		v1 [1] = "amino_acids";
		v1 [2] = aminoAcids;
		insertTable ( tableName, n1, v1 );

		v1 [1] = "major_flag";
		v1 [2] = genNextString ( str, "|", start, end );
		insertTable ( tableName, n1, v1 );

		v1 [1] = "immonium_flag";
		v1 [2] = genNextString ( str, "|", start, end );
		insertTable ( tableName, n1, v1 );

		v1 [1] = "exclude_amino_acids";
		v1 [2] = str.substr ( start );


		insertTable ( tableName, n1, v1 );
	}
}
class SQLiteIndicies : public SQLiteIDNameValueList {
public:
	SQLiteIndicies ();
	void insert ();
};
SQLiteIndicies::SQLiteIndicies () :
	SQLiteIDNameValueList ( "Indicies", "indicies.txt" )
{
}
void SQLiteIndicies::insert ()
{
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( fileName ) );
	while ( ifs.getUncommentedLine ( v1 [0] ) ) {
		string line;
		while ( ifs.getUncommentedLine ( line ) ) {
			if ( line.length () != 0 ) {
				if ( line [0] == '>' ) break;
				istringstream ist ( line );
				ist >> v1 [1];
				string temp;
				getline ( ist, temp );
				v1 [2] = temp.substr ( temp.find_first_not_of ( " \t" ) );	// strip leading space
				insertTable ( tableName, n1, v1 );
			}
		}
	}
}
class SQLiteInfo : public SQLiteNameValueList {
public:
	SQLiteInfo ();
	void insert ();
};
SQLiteInfo::SQLiteInfo () :
	SQLiteNameValueList ( "Info", "" )
{
}
void SQLiteInfo::insert ()
{
	VectorPairStringString vpss;
	vpss.push_back ( make_pair ( string("seqdb"),							string("seqdb")		) );
	vpss.push_back ( make_pair ( string("database_suffix"),					string(".fasta")	) );
	vpss.push_back ( make_pair ( string("upload_temp"),						string("temp")		) );
	vpss.push_back ( make_pair ( string("max_upload_content_length"),		string("0")			) );
	vpss.push_back ( make_pair ( string("r_command"),						string("/usr/bin/R")) );
	//vpss.push_back ( make_pair ( string("r_command_win"),					string("C:\\Program Files\\R\\R-2.14.1\\bin\\R")) );
	//vpss.push_back ( make_pair ( string("r_command_unix"),				string("/usr/bin/R")) );
	vpss.push_back ( make_pair ( string("ucsf_banner"),						string("true")		) );
	vpss.push_back ( make_pair ( string("mssearch_logging"),				string("false")		) );
	vpss.push_back ( make_pair ( string("mssearch_parameter_logging"),		string("false")		) );
	vpss.push_back ( make_pair ( string("msform_logging"),					string("false")		) );
	vpss.push_back ( make_pair ( string("msform_parameter_logging"),		string("false")		) );
	vpss.push_back ( make_pair ( string("msbatch_logging"),					string("false")		) );
	vpss.push_back ( make_pair ( string("msbatch_parameter_logging"),		string("false")		) );
	vpss.push_back ( make_pair ( string("msdisplay_logging"),				string("false")		) );
	vpss.push_back ( make_pair ( string("msdisplay_parameter_logging"),		string("false")		) );
	vpss.push_back ( make_pair ( string("searchCompare_logging"),			string("false")		) );
	vpss.push_back ( make_pair ( string("searchCompare_parameter_logging"),	string("false")		) );
	vpss.push_back ( make_pair ( string("delete_log_days"),					string("0")			) );
	vpss.push_back ( make_pair ( string("timeout"),							string("0")			) );
	vpss.push_back ( make_pair ( string("max_msprod_sequences"),			string("2")			) );
	vpss.push_back ( make_pair ( string("max_msfit_peaks"),					string("1000")		) );
	vpss.push_back ( make_pair ( string("msfit_max_reported_hits_limit"),	string("500")		) );
	vpss.push_back ( make_pair ( string("faindex_parallel"),				string("false")		) );
	//vpss.push_back ( make_pair ( string("viewer_repository"),				string("")			) );
	//vpss.push_back ( make_pair ( string("centroid_dir"),					string("")			) );
	//vpss.push_back ( make_pair ( string("centroid_dir_win"),					string("")			) );
	//vpss.push_back ( make_pair ( string("centroid_dir_unix"),					string("")			) );
	//vpss.push_back ( make_pair ( string("raw_dir"),						string("")			) );
	//vpss.push_back ( make_pair ( string("raw_dir_win"),						string("")			) );
	//vpss.push_back ( make_pair ( string("raw_dir_unix"),						string("")			) );
	//vpss.push_back ( make_pair ( string("user_repository"),					string("")			) );
	//vpss.push_back ( make_pair ( string("user_repository_win"),					string("")			) );
	//vpss.push_back ( make_pair ( string("user_repository_unix"),					string("")			) );
	vpss.push_back ( make_pair ( string("multi_process"),					string("false")		) );
	vpss.push_back ( make_pair ( string("msms_max_spectra"),				string("500")		) );
	vpss.push_back ( make_pair ( string("duplicate_scans"),					string("false")		) );
	//vpss.push_back ( make_pair ( string("mpi_run"),							string("")			) );
	//vpss.push_back ( make_pair ( string("mpi_args"),						string("")			) );
	vpss.push_back ( make_pair ( string("min_password_length"),				string("6")			) );
	vpss.push_back ( make_pair ( string("db_host"),							string("localhost")	) );
	vpss.push_back ( make_pair ( string("db_port"),							string("0")			) );
	vpss.push_back ( make_pair ( string("db_name"),							string("ppsd")		) );
	vpss.push_back ( make_pair ( string("db_user"),							string("prospector")) );
	vpss.push_back ( make_pair ( string("db_password"),						string("pp")		) );
	//vpss.push_back ( make_pair ( string("btag_daemon_name"),				string("")			) );
	vpss.push_back ( make_pair ( string("btag_daemon_remote"),				string("false")		) );
	vpss.push_back ( make_pair ( string("max_btag_searches"),				string("1")			) );
	vpss.push_back ( make_pair ( string("email"),							string("false")		) );
	vpss.push_back ( make_pair ( string("server_name"),						string("localhost")	) );
	vpss.push_back ( make_pair ( string("server_port"),						string("80")		) );
	//vpss.push_back ( make_pair ( string("virtual_dir"),						string("")			) );
	vpss.push_back ( make_pair ( string("job_status_refresh_time"),			string("5")			) );
	vpss.push_back ( make_pair ( string("daemon_loop_time"),				string("5")			) );
	vpss.push_back ( make_pair ( string("single_server"),					string("false")		) );
	vpss.push_back ( make_pair ( string("aborted_jobs_delete_days"),		string("0")			) );
	vpss.push_back ( make_pair ( string("session_delete_days"),				string("0")			) );
	//vpss.push_back ( make_pair ( string("preload_database"),				string("")			) );
	vpss.push_back ( make_pair ( string("join_results_files"),				string("true")		) );
	vpss.push_back ( make_pair ( string("expectation_search_first"),		string("false")		) );
	vpss.push_back ( make_pair ( string("raw_data_forwarding"),				string("false")		) );
	//vpss.push_back ( make_pair ( string("virtual_dir_proxy"),				string("")			) );
	v1 [2] = "";
	v1 [3] = "1";
	v1 [4] = "0";
	for ( VectorPairStringStringSizeType i = 0 ; i < vpss.size () ; i++ ) {	// This is the default set
		string n = vpss [i].first;
		v1 [0] = n;
		v1 [1] = vpss [i].second;
		insertTable ( "Info", n1, v1 );
	}
}
/*	Could be useful for creating the actual set
	for ( VectorPairStringStringSizeType i = 0 ; i < vpss.size () ; i++ ) {
		string n = vpss [i].first;
		v1 [0] = n;
		v1 [1] = InfoParams::instance ().getStringValue ( n );
		v1 [2] = vpss [i].second;
		insertTable ( "Info", n1, v1 );
	}
*/
class SQLiteMGF : public SQLiteIDNameValueList {
public:
	SQLiteMGF ();
	void insert ();
};
SQLiteMGF::SQLiteMGF () :
	SQLiteIDNameValueList ( "MGF", "mgf.xml" )
{
}
void SQLiteMGF::insert ()
{
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	char* info = getFileAsCharPtr ( MsparamsDir::instance ().getParamPath ( fileName ) );

	StringVector xparams = XMLParser::getStringVectorValue ( info, "mgf_type" );
	for ( StringVectorSizeType i = 0 ; i < xparams.size () ; i++ ) {
		ParameterList pList ( xparams [i], false, false, false );
		StringVector pNames;
		StringVectorVector pValues;
		pList.getParams ( pNames, pValues );
		v1 [0] = pList.getStringValue ( "name" );
		for ( StringVectorSizeType j = 0 ; j < pNames.size () ; j++ ) {
			if ( pNames [j] != "name" ) {
				v1 [1] = pNames [j];
				for ( StringVectorSizeType k = 0 ; k < pValues [j].size () ; k++ ) {
					v1 [2] = pValues [j] [k];
					insertTable ( tableName, n1, v1 );
				}
			}
		}
	}

	delete [] info;
}

class SQLiteModifications : public SQLiteIDNameValueList {
public:
	SQLiteModifications ();
	void insert ();
};
SQLiteModifications::SQLiteModifications () :
	SQLiteIDNameValueList ( "Mods", "" )
{
}
void SQLiteModifications::insert ()
{
	StringVector namesList;
	StringVector specificityList;
	StringVector elementalFormulaStrList;
	StringVector typeList;
	Usermod::getAllUsermodInfo ( namesList, specificityList, elementalFormulaStrList, typeList );
	MapStringToInt modTypes;
	int modTypeIdx = 0;
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	SetString ss;
	for ( int i = 0 ; i < namesList.size () ; i++ ) {
		PairSetStringIteratorBool flag = ss.insert ( namesList [i] );
		if ( flag.second ) {
			v1 [0] = namesList [i];
			v1 [1] = "formula";
			v1 [2] = elementalFormulaStrList [i];

			insertTable ( tableName, n1, v1 );
		}
	}
}

class SQLiteModificationSpecificity : public SQLiteIDSubIDNameValueList {
public:
	SQLiteModificationSpecificity ();
	void insert ();
};
SQLiteModificationSpecificity::SQLiteModificationSpecificity () :
	SQLiteIDSubIDNameValueList ( "ModSpecificities", "" )
{
}
void SQLiteModificationSpecificity::insert ()
{
	StringVector namesList;
	StringVector specificityList;
	StringVector elementalFormulaStrList;
	StringVector typeList;
	Usermod::getAllUsermodInfo ( namesList, specificityList, elementalFormulaStrList, typeList );
	MapStringToInt modTypes;
	int modTypeIdx = 0;
	v1 [4] = "";
	v1 [5] = "1";
	v1 [6] = "0";
	for ( int i = 0 ; i < namesList.size () ; i++ ) {
		v1 [0] = namesList [i];
		v1 [1] = specificityList [i];

		v1 [2] = "type";
		v1 [3] = typeList [i];
		insertTable ( tableName, n1, v1 );
	}
}
#ifdef RAW_DATA
class SQLiteQuan : public SQLiteNonUniqueNameValueList {
public:
	SQLiteQuan ();
	void create ();
	void insert ();
};
SQLiteQuan::SQLiteQuan () :
	SQLiteNonUniqueNameValueList ( "Quan", "quan.txt" )
{
}
void SQLiteQuan::create ()
{
	StringVector sql;
	sql.push_back ( "CREATE TABLE " + tableName + " ( "\
			+ n1 [0] + " TEXT NOT NULL, "\
			+ n1 [1] + " TEXT, "\
			+ n1 [2] + " TEXT, "\
			+ n1 [3] + " INTEGER CHECK ( " + n1 [3] + " = 0 OR " + n1 [3] + " = 1 ), "\
			+ n1 [4] + " INTEGER CHECK ( " + n1 [4] + " = 0 OR " + n1 [4] + " = 1 ), "\
			+ "PRIMARY KEY ( "\
				+ n1 [0] + ", "\
				+ n1 [1]\
			+ " ) "\
		+ ")" );
	createTables ( sql );
}
void SQLiteQuan::insert ()
{
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( fileName ) );
	v1 [2] = "";
	v1 [3] = "1";
	v1 [4] = "0";
	while ( ifs.getUncommentedLine ( v1 [0] ) ) {
		string line;
		int num = 0;
		while ( ifs.getUncommentedLine ( line ) ) {
			if ( line [0] == '>' ) {
				if ( num == 0 ) {
					v1 [1] = "";						// Covers the N15 case
					insertTable ( tableName, n1, v1 );
				}
				break;
			}
			v1 [1] = genEscapeQuote ( line );
			insertTable ( tableName, n1, v1 );
			num++;
		}
	}
}

class SQLiteQuanMSMS : public SQLiteIDNameValueList {
public:
	SQLiteQuanMSMS ();
	void insert ();
};
SQLiteQuanMSMS::SQLiteQuanMSMS () :
	SQLiteIDNameValueList ( "QuanMSMS", "" )
{
}
void SQLiteQuanMSMS::insert ()
{
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	StringVector quanTypes = QuanMSMSXMLData::instance ().getQuanMSMSNames ();
	for ( StringVectorSizeType i = 0 ; i < quanTypes.size () ; i++ ) {
		v1 [0] = quanTypes [i];
		QuanMSMSInfo qmi = QuanMSMSXMLData::instance ().getQuanMSMSInfo ( v1 [0] );
		StringVectorSizeType siz = qmi.formulae.size ();
		for ( StringVectorSizeType j = 0 ; j < siz ; j++ ) {
			v1 [1] = qmi.quanPeaks [j] ? "reporter_ion" : "contaminant_ion";
			v1 [2] = qmi.formulae [j];
			insertTable ( tableName, n1, v1 );
		}
	}
}
class SQLiteQuanMSMSPurity : public SQLiteIDSubIDSubIDNameValueList {
public:
	SQLiteQuanMSMSPurity ();
	void insert ();
};
SQLiteQuanMSMSPurity::SQLiteQuanMSMSPurity () :
	SQLiteIDSubIDSubIDNameValueList ( "QuanMSMSPurity", "" )
{
}
void SQLiteQuanMSMSPurity::insert ()
{
	v1 [5] = "";
	v1 [6] = "1";
	v1 [7] = "0";
	StringVector quanTypes = QuanMSMSXMLData::instance ().getQuanMSMSNames ();

	for ( StringVectorSizeType i = 0 ; i < quanTypes.size () ; i++ ) {
		v1 [0] = quanTypes [i];
		string purityPath = MsparamsDir::instance ().getParamPath (  v1 [0] + ".txt" );
		if ( genFileExists ( purityPath ) ) {
			GenCommentedIFStream ifs ( purityPath );
			while ( ifs.getUncommentedLine ( v1 [1] ) ) {
				ParameterList p ( ifs );
				StringVector sv = p.getNameList ();
				int numQuanPeaks = p.size ();
				for ( int j = 0 ; j < numQuanPeaks ; j++ ) {
					v1 [2] = sv [j];
					istringstream ist ( p.getStringValue ( sv [j] ) );
					for ( int k = 1 ; k <= 4 ; k++ ) {
						v1 [3] = "coefficient" + gen_itoa ( k );
						ist >> v1 [4];
						insertTable ( tableName, n1, v1 );
					}
				}
			}
		}
	}
}
#endif

class SQLiteRepository : public SQLiteIDSubIDNameValueList {
public:
	SQLiteRepository ();
	void insert ();
};
SQLiteRepository::SQLiteRepository () :
	SQLiteIDSubIDNameValueList ( "Repository", "" )
{
}
void SQLiteRepository::insert ()
{
	v1 [4] = "";
	v1 [5] = "1";
	v1 [6] = "0";
	if ( RepositoryInfo::exists () ) {
		StringVector sv = RepositoryInfo::instance ().getNames ();
		int paramIdx = 1;
		for ( int i = 0 ; i < sv.size () ; i++ ) {
			v1 [0] = sv [i];
			MapPairStringIntInstrumentType mpsi = RepositoryInfo::instance ().getInstrumentInfo ( v1 [0] );
			for ( MapPairStringIntInstrumentTypeConstIterator j = mpsi.begin () ; j != mpsi.end () ; j++ ) {
				string suffix = (*j).second.getFileSuffix ();
				v1 [1] = suffix.empty () ? "default" : suffix;

				v1 [2] = "instrument";
				v1 [3] = (*j).first.first;
				insertTable ( tableName, n1, v1 );

				MapStringToStringVector mssv = (*j).second.getParams ();
				if ( !mssv.empty () ) {
					v1 [2] = "parameters";
					v1 [3] = gen_itoa ( paramIdx++ );
					insertTable ( tableName, n1, v1 );
				}
			}
		}
	}
	else {
		string instDirPath = MsparamsDir::instance ().getParamPath ( "inst_dir.txt" );
		if ( genFileExists ( instDirPath ) ) {
			GenNameValueStream nvs ( MsparamsDir::instance ().getParamPath ( instDirPath ) );
			StringVector sv = nvs.getNameList ();
			for ( int i = 0 ; i < sv.size () ; i++ ) {
				v1 [0] = sv [i];
				v1 [1] = "default";
				v1 [2] = "instrument";
				v1 [3] = nvs.getStringValue ( sv [i] );
				insertTable ( tableName, n1, v1 );
			}
		}
	}
}
class SQLiteRepositoryParams : public SQLiteParameterLists {
public:
	SQLiteRepositoryParams ();
	void create () {}			// Don't create a new table
	void insert ();
};
SQLiteRepositoryParams::SQLiteRepositoryParams () :
	SQLiteParameterLists ( "repository.xml" )
{
}
void SQLiteRepositoryParams::insert ()
{
	v1 [4] = "";
	v1 [5] = "1";
	v1 [6] = "0";
	if ( RepositoryInfo::exists () ) {
		StringVector sv = RepositoryInfo::instance ().getNames ();
		int paramIdx = 1;
		for ( int i = 0 ; i < sv.size () ; i++ ) {
			MapPairStringIntInstrumentType mpsi = RepositoryInfo::instance ().getInstrumentInfo ( sv [i] );
			for ( MapPairStringIntInstrumentTypeConstIterator j = mpsi.begin () ; j != mpsi.end () ; j++ ) {
				string suffix = (*j).second.getFileSuffix ();
				MapStringToStringVector mssv = (*j).second.getParams ();
				if ( !mssv.empty () ) {
					v1 [0] = "repository";
					v1 [1] = gen_itoa ( paramIdx++ );
					for ( MapStringToStringVectorConstIterator k = mssv.begin () ; k != mssv.end () ; k++ ) {
						v1 [2] = (*k).first;
						StringVector sv2 = (*k).second;
						for ( StringVectorSizeType l = 0 ; l < sv2.size () ; l++ ) {
							v1 [3] = sv2 [l];
							insertTable ( "ParameterLists", n1, v1 );
						}
					}
				}
			}
		}
	}
}
class SQLiteTaxonomyGroups : public SQLiteNonUniqueNameValueList {
public:
	SQLiteTaxonomyGroups ();
	void insert ();
};
SQLiteTaxonomyGroups::SQLiteTaxonomyGroups () :
	SQLiteNonUniqueNameValueList ( "TaxonomyGroups", "taxonomy_groups.txt" )
{
}
void SQLiteTaxonomyGroups::insert ()
{
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( fileName ) );
	v1 [2] = "";
	v1 [3] = "1";
	v1 [4] = "0";
	while ( ifs.getUncommentedLine ( v1 [0] ) ) {
		string line;
		while ( ifs.getUncommentedLine ( line ) ) {
			if ( line [0] == '>' ) break;
			v1 [1] = genEscapeQuote ( line );
			insertTable ( tableName, n1, v1 );
		}
	}
}
class SQLiteViewerConv : public SQLiteIDNameValueList {
public:
	SQLiteViewerConv ();
	void insert ();
};
SQLiteViewerConv::SQLiteViewerConv () :
	SQLiteIDNameValueList ( "ViewerConv", "viewer_conv.txt" )
{
}
void SQLiteViewerConv::insert ()
{
	v1 [3] = "";
	v1 [4] = "1";
	v1 [5] = "0";
	GenCommentedIFStream ifs ( MsparamsDir::instance ().getParamPath ( fileName ) );
	string line;
	int phase = 0;
	while ( ifs.getUncommentedLine ( line ) ) {
		if ( phase == 0 ) v1 [0] = line;
		if ( phase == 1 ) v1 [1] = "conversion_script";
		if ( phase == 2 ) v1 [1] = "num_title_lines";
		if ( phase == 3 ) v1 [1] = "num_header_lines";
		if ( phase == 4 ) v1 [1] = "column_separator";
		if ( phase == 5 ) v1 [1] = "spectrum_identifier";
		if ( phase == 6 ) v1 [1] = "fraction_header";
		if ( phase == 7 ) v1 [1] = "scan_id_header";
		if ( phase == 8 ) v1 [1] = "peptide_header";
		if ( phase == 9 ) v1 [1] = "charge_header";
		if ( phase == 10 ) v1 [1] = "modifications_header";
		if ( phase != 0 ) {
			v1 [2] = line;
			insertTable ( tableName, n1, v1 );
		}
		phase = ( phase == 10 ) ? 0 : ++phase;
	}
}


ParamDB::ParamDB () :
	PPSQLite ()
{
}
ParamDB::~ParamDB ()
{
}

ParamDBWrite::ParamDBWrite ( const string& name, bool append ) :
	ParamDB ()
{
	if ( !append ) {		// Get rid of the old database if creating a new one
		genUnlink ( name );
	}
	rc = sqlite3_open_v2 ( name.c_str (), &db, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL );
}
ParamDBWrite::~ParamDBWrite ()
{
}
void ParamDBWrite::create ()
{
	SQLiteTable::setDatabaseConnection ( this );
	vector <SQLiteTable*> sqt;
	sqt.push_back ( new SQLiteAATable );
	sqt.push_back ( new SQLiteNameValueList ( "AccLinks", "acclinks.txt" ) );
	sqt.push_back ( new SQLiteApplets );
	sqt.push_back ( new SQLiteList ( "B1", "b1.txt" ) );
	sqt.push_back ( new SQLiteNonUniqueNameValueList ( "ChemScore", "chem_score.txt" ) );
	sqt.push_back ( new SQLiteList ( "DbestSpl", "dbEST.spl.txt" ) );
	sqt.push_back ( new SQLiteDBHosts );
	sqt.push_back ( new SQLiteIDNameValueList ( "DiscScore", "disc_score.txt" ) );
	sqt.push_back ( new SQLiteIDNameValueList ( "DiscScore2", "disc_score2.txt" ) );
	sqt.push_back ( new SQLiteIDUnsortedNameValueList ( "Distribution", "distribution.txt" ) );
	sqt.push_back ( new SQLiteElements );
	sqt.push_back ( new SQLiteEnzymes );
	sqt.push_back ( new SQLiteEnzymeCombinations );
	sqt.push_back ( new SQLiteExpectation );
	sqt.push_back ( new SQLiteNameValueList ( "Expectation", "expectation.txt" ) );
	sqt.push_back ( new SQLiteIDNameValueList ( "Fragmentation", "fragmentation.txt" ) );
	sqt.push_back ( new SQLiteIDNameValueList ( "Homology", "homology.txt" ) );
	sqt.push_back ( new SQLiteIdxLinks );
	sqt.push_back ( new SQLiteImmonium );
	sqt.push_back ( new SQLiteIndicies );
	sqt.push_back ( new SQLiteInfo );
	sqt.push_back ( new SQLiteIDNameValueList ( "Instrument", "instrument.txt" ) );
	sqt.push_back ( new SQLiteList ( "LinkAA", "link_aa.txt" ) );
	sqt.push_back ( new SQLiteIDNameValueList ( "Links", "links.txt" ) );
	sqt.push_back ( new SQLiteIDUnsortedNameValueList ( "MatScore", "mat_score.txt" ) );
	sqt.push_back ( new SQLiteMGF () );
	sqt.push_back ( new SQLiteModifications );
	sqt.push_back ( new SQLiteModificationSpecificity );
	sqt.push_back ( new SQLiteParameterLists );
#ifdef RAW_DATA
	sqt.push_back ( new SQLiteQuan );
	sqt.push_back ( new SQLiteQuanMSMS );
	sqt.push_back ( new SQLiteQuanMSMSPurity );
#endif
	sqt.push_back ( new SQLiteRepository );
	sqt.push_back ( new SQLiteRepositoryParams );
	sqt.push_back ( new SQLiteList ( "Taxonomy", "taxonomy.txt" ) );
	sqt.push_back ( new SQLiteTaxonomyGroups );
	sqt.push_back ( new SQLiteNameValueList ( "Unimod", "unimod.txt" ) );
	sqt.push_back ( new SQLiteNameValueList ( "UpIdLinks", "upidlinks.txt" ) );
	sqt.push_back ( new SQLiteViewerConv );
	beginTransaction ();
	for ( int i = 0 ; i < sqt.size () ; i++ ) {								// create the tables
		sqt [i]->create ();
	}
	endTransaction ();
	beginTransaction ();
	for ( int j = 0 ; j < sqt.size () ; j++ ) {								// fill the tables
		sqt [j]->insert ();
	}
	endTransaction ();
}

ParamDBRead::ParamDBRead ( const string& name ) :
	ParamDB ()
{
	rc = sqlite3_open_v2 ( name.c_str (), &db, SQLITE_OPEN_READONLY, NULL );
}
ParamDBRead::~ParamDBRead ()
{
}
