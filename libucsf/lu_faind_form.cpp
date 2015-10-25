/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_faind_form.cpp                                             *
*                                                                             *
*  Created    : January 11th 2005                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_html.h>
#include <lu_table.h>
#include <lu_faind_form.h>
using std::string;
using std::vector;
using std::ostream;
using std::endl;

class FormItemInputFAIndexDatabase : public FormItemText {
public:
	FormItemInputFAIndexDatabase () :
		FormItemText ( "", "", getName (), 50, 100, "Enter filename here" ) {}
	static string getName () { return "database"; }
};

class FormItemDNAToProtein : public FormItemCheckbox {
public:
	FormItemDNAToProtein () :
		FormItemCheckbox ( "DNA to Protein", "", getName (), false ) {}
	static string getName () { return "dna_to_protein"; }
};

class FormItemDeleteDNADatabase : public FormItemCheckbox {
public:
	FormItemDeleteDNADatabase () :
		FormItemCheckbox ( "Delete DNA Database", "", getName (), false ) {}
	static string getName () { return "delete_dna_database"; }
};

class FormItemRandomDatabase : public FormItemCheckbox {
public:
	FormItemRandomDatabase () :
		FormItemCheckbox ( "Random Database", "", getName (), false ) {}
	static string getName () { return "random_database"; }
};

class FormItemReverseDatabase : public FormItemCheckbox {
public:
	FormItemReverseDatabase () :
		FormItemCheckbox ( "Reverse Database", "", getName (), false ) {}
	static string getName () { return "reverse_database"; }
};

class FormItemConcatDatabase : public FormItemCheckbox {
public:
	FormItemConcatDatabase () :
		FormItemCheckbox ( "Concat Database", "", getName (), false ) {}
	static string getName () { return "concat_database"; }
};

class FormItemCreateDatabaseIndicies : public FormItemCheckbox {
public:
	FormItemCreateDatabaseIndicies ( bool val ) :
		FormItemCheckbox ( "Create Database Indicies", "", getName (), val ) {}
	static string getName () { return "create_database_indicies"; }
};

class FormItemCreateSubDatabase : public FormItemCheckbox {
public:
	FormItemCreateSubDatabase () :
		FormItemCheckbox ( "Create Sub Database", "", getName (), true ) {}
	static string getName () { return "create_sub_database"; }
};

class FormItemSubDatabaseID : public FormItemText {
public:
	FormItemSubDatabaseID () :
		FormItemText ( "Suffix for Subset Database", "", getName (), 15, 100, "sub" ) {}
	static string getName () { return "sub_database_id"; }
};

class FormItemFAIndexUserDB : public FormItemText {
public:
	FormItemFAIndexUserDB () :
		FormItemText ( "Database Name", "", getName (), 20, 100, "NCBInr.user" ) {}
	static string getName () { return "database"; }
};

class FormItemFAIndexUserName : public FormItemText {
public:
	FormItemFAIndexUserName () :
		FormItemText ( "Name", "", getName (), 110, 240, "Better than sliced bread growth factor" ) {}
	static string getName () { return "name_field"; }
};

class FormItemFAIndexUserSpecies : public FormItemText {
public:
	FormItemFAIndexUserSpecies () :
		FormItemText ( "Species", "", getName (), 30, 100, "HUMAN" ) {}
	static string getName () { return "species"; }
};

class FormItemFAIndexUserAccessionNumber : public FormItemText {
public:
	FormItemFAIndexUserAccessionNumber () :
		FormItemText ( "Accession Number", "", getName (), 20, 100, "1" ) {}
	static string getName () { return "accession_num"; }
};

class FormItemCreateUserDatabase : public FormItemCheckbox {
public:
	FormItemCreateUserDatabase () :
		FormItemCheckbox ( "Create User Database", "", getName (), true ) {}
	static string getName () { return "create_user_database"; }
};

class FormItemStartIndexNumber : public FormItemText {
public:
	FormItemStartIndexNumber () :
		FormItemText ( "Index Number Range", "", getName (), 8, 8, "1" ) {}
	static string getName () { return "start_index_number"; }
};

class FormItemEndIndexNumber : public FormItemText {
public:
	FormItemEndIndexNumber () :
		FormItemText ( "to", "", getName (), 8, 8, "1" ) {}
	static string getName () { return "end_index_number"; }
};

class FormItemAllIndicies : public FormItemCheckbox {
public:
	FormItemAllIndicies () :
		FormItemCheckbox ( "All", "", getName (), false ) {}
	static string getName () { return "all_indicies"; }
};

FAIndexForm::FAIndexForm ( const VectorConstParameterListPtr& params ) :
	proteinMWForm ( &fvj, params ),
	proteinPIForm ( &fvj, params ),
	searchResultsForm ( params, true )
{
	create ( params );
}
void FAIndexForm::createItems ()
{
	formItemMap [FormItemVersion::getName ()]	= new FormItemVersion ();
	formItemMap [FormItemInputFAIndexDatabase::getName()+"1"]= new FormItemInputFAIndexDatabase ();
	formItemMap [FormItemDNAToProtein::getName()]			= new FormItemDNAToProtein ();
	formItemMap [FormItemDeleteDNADatabase::getName()]		= new FormItemDeleteDNADatabase ();
	formItemMap [FormItemRandomDatabase::getName()]			= new FormItemRandomDatabase ();
	formItemMap [FormItemReverseDatabase::getName()]		= new FormItemReverseDatabase ();
	formItemMap [FormItemConcatDatabase::getName()]			= new FormItemConcatDatabase ();
	formItemMap [FormItemCreateDatabaseIndicies::getName()]	= new FormItemCreateDatabaseIndicies ( true );

	formItemMap [FormItemCreateSubDatabase::getName()]	= new FormItemCreateSubDatabase ();
	formItemMap [FormItemSubDatabaseID::getName()]		= new FormItemSubDatabaseID ();
	formItemMap [FormItemDatabase::getName ()+"2"]		= new FormItemDatabase ();
	formItemMap [FormItemTaxonomy::getName ()+"1"]		= new FormItemTaxonomy ();
	formItemMap [FormItemTaxonomyRemove::getName ()]	= new FormItemTaxonomyRemove ();
	formItemMap [FormItemTaxonomyNames::getName ()]		= new FormItemTaxonomyNames ();
	formItemMap [FormItemNames::getName ()]				= new FormItemNames ();
	formItemMap [FormItemAccessionNums::getName ()]		= new FormItemAccessionNums ();

	formItemMap [FormItemFAIndexUserDB::getName ()+"3"]			= new FormItemFAIndexUserDB ();
	formItemMap [FormItemFAIndexUserName::getName ()]			= new FormItemFAIndexUserName ();
	formItemMap [FormItemFAIndexUserSpecies::getName ()+"2"]	= new FormItemFAIndexUserSpecies ();
	formItemMap [FormItemFAIndexUserAccessionNumber::getName ()]= new FormItemFAIndexUserAccessionNumber ();

	formItemMap [FormItemCreateUserDatabase::getName ()]	= new FormItemCreateUserDatabase ();
	formItemMap [FormItemUserProteinSequence::getName ()]	= new FormItemUserProteinSequence ( true );

// Database summary report
	formItemMap [FormItemDNAReadingFrame::getName ()]	= new FormItemDNAReadingFrame ();
	formItemMap [FormItemStartIndexNumber::getName ()]	= new FormItemStartIndexNumber ();
	formItemMap [FormItemEndIndexNumber::getName ()]	= new FormItemEndIndexNumber ();
	formItemMap [FormItemAllIndicies::getName ()]		= new FormItemAllIndicies ();
	formItemMap [FormItemHideProteinSequence::getName ()]= new FormItemHideProteinSequence ( true );
}
void FAIndexForm::printHTML ( ostream& os )
{
	init_html ( os, "FA-Index", "faman.htm" );
	fvj.print ( os );
	tableStart ( os, false );
		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMStart ( os, "post", "faindex" );
				tableStart ( os, true );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", true );
							printHTMLFORMSubmit ( os, "Create Indicies For New Database" );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							os << "Newly Downloaded <a href=\"../html/instruct/allman.htm#database\">Database:</a>" << endl;
							os << "<br />" << endl;
							formItemMap [FormItemInputFAIndexDatabase::getName()+"1"]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemDNAToProtein::getName()]->printHTML ( os );
							formItemMap [FormItemDeleteDNADatabase::getName()]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemRandomDatabase::getName()]->printHTML ( os );
							formItemMap [FormItemReverseDatabase::getName()]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemConcatDatabase::getName()]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
				formItemMap [FormItemCreateDatabaseIndicies::getName()]->printHTMLHidden ( os );
				formItemMap [FormItemVersion::getName()]->printHTMLHidden ( os );
				formItemMap [FormItemTaxonomy::getName ()+"1"]->printHTMLHidden ( os );
				printHTMLFORMEnd ( os );
			tableHeaderEnd ( os );

			tableHeaderStart ( os, "", "center", true );
				printHTMLFORMStart ( os, "post", "faindex" );
				tableStart ( os, true );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", true );
							printHTMLFORMSubmit ( os, "Database Summary Report" );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							formItemMap [FormItemDatabase::getName()+"2"]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemDNAReadingFrame::getName()]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemStartIndexNumber::getName ()]->printHTML ( os );
							formItemMap [FormItemEndIndexNumber::getName ()]->printHTML ( os );
							formItemMap [FormItemAllIndicies::getName ()]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemHideProteinSequence::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
				formItemMap [FormItemVersion::getName()]->printHTMLHidden ( os );
				printHTMLFORMEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				printHTMLFORMStart ( os, "post", "faindex", false, true );
				tableStart ( os, true );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							printHTMLFORMSubmit ( os, "Create Subset Database" );
							formItemMap [FormItemSubDatabaseID::getName()]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemDatabase::getName()+"2"]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemTaxonomy::getName()+"1"]->printHTML ( os );
							formItemMap [FormItemTaxonomyNames::getName ()]->printHTML ( os );
							formItemMap [FormItemTaxonomyRemove::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
						tableHeaderStart ( os, "", "left", true );
							proteinMWForm.printHTML ( os );
							os << "<br />" << endl;
							proteinPIForm.printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemNames::getName ()]->printHTML ( os );
							formItemMap [FormItemAccessionNums::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
				formItemMap [FormItemCreateDatabaseIndicies::getName()]->printHTMLHidden ( os );
				formItemMap [FormItemCreateSubDatabase::getName()]->printHTMLHidden ( os );
				formItemMap [FormItemVersion::getName()]->printHTMLHidden ( os );
				printHTMLFORMEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				printHTMLFORMStart ( os, "post", "faindex" );
				tableStart ( os, true );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", true );
							printHTMLFORMSubmit ( os, "Create Subset Database with Indices from Saved Hits" );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							formItemMap [FormItemSubDatabaseID::getName()]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemDatabase::getName()+"2"]->printHTML ( os );
							os << "<br />" << endl;
							searchResultsForm.printHTMLFAIndex ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
				formItemMap [FormItemCreateDatabaseIndicies::getName()]->printHTMLHidden ( os );
				formItemMap [FormItemCreateSubDatabase::getName()]->printHTMLHidden ( os );
				formItemMap [FormItemVersion::getName()]->printHTMLHidden ( os );
				printHTMLFORMEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );

		tableRowStart ( os );
			tableHeaderStart ( os, "", "center", true, 2 );
				printHTMLFORMStart ( os, "post", "faindex" );
				tableStart ( os, true );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "center", true );
							printHTMLFORMSubmit ( os, "Create or Append to User Database" );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
					tableRowStart ( os );
						tableHeaderStart ( os, "", "left", true );
							formItemMap [FormItemFAIndexUserDB::getName()+"3"]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemFAIndexUserName::getName()]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemFAIndexUserSpecies::getName()+"2"]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemFAIndexUserAccessionNumber::getName()]->printHTML ( os );
							os << "<br />" << endl;
							formItemMap [FormItemUserProteinSequence::getName ()]->printHTML ( os );
						tableHeaderEnd ( os );
					tableRowEnd ( os );
				tableEnd ( os );
				formItemMap [FormItemCreateDatabaseIndicies::getName()]->printHTMLHidden ( os );
				formItemMap [FormItemCreateUserDatabase::getName()]->printHTMLHidden ( os );
				formItemMap [FormItemVersion::getName()]->printHTMLHidden ( os );
				printHTMLFORMEnd ( os );
			tableHeaderEnd ( os );
		tableRowEnd ( os );
	tableEnd ( os );

	printHTMLFooter ( os, "1995" );
}
