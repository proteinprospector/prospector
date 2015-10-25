/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_form_valid.cpp                                             *
*                                                                             *
*  Created    : March 12th 2008                                               *
*                                                                             *
*  Purpose    : Functions for validating HTML form input.                     *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2008-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <iostream>
#include <lu_html.h>
#include <lu_form_valid.h>

using std::ostream;
using std::string;
using std::endl;

// ? - zero or one occurrences
// * - zero or more occurrences
// + - one or more occurrences

string FormValidatingJavascript::positiveNonZeroIntegerRegEx	= "/^([1-9][0-9]*)$/";
string FormValidatingJavascript::positiveIntegerRegEx			= "/^(0|[1-9][0-9]*)$/";
//string FormValidatingJavascript::positiveDecimalRegEx			= "/(^(0?|[1-9][0-9]*)\\.(0*[1-9][0-9]*)$)|(^[1-9]+[0-9]*\\.0+$)|(^0\\.0+$)/";
//string FormValidatingJavascript::positiveIntegerOrDecimalRegEx	= "/(^(0|[1-9][0-9]*)$)|((^(0?|[1-9][0-9]*)\\.(0*[1-9][0-9]*)$)|(^[1-9]+[0-9]*\\.0+$)|(^0\\.0+$))/";
//string FormValidatingJavascript::signedIntegerOrDecimalRegEx	= "/(^[+]?0(\\.0+)?$)|(^([-+]?[1-9][0-9]*)$)|(^([-+]?((0?|[1-9][0-9]*)\\.(0*[1-9][0-9]*)))$)|(^[-+]?[1-9]+[0-9]*\\.0+$)/";
string FormValidatingJavascript::positiveFloatingPointRegEx		= "/^([0-9]*\\.?[0-9]+)$/";
string FormValidatingJavascript::signedFloatingPointRegEx		= "/^([-+]?[0-9]*\\.?[0-9]+)$/";

string FormValidatingJavascript::elementalModificationFormulaAllowBlankRegEx = "/^$|^((([1-9][0-9]*)?[A-Z][a-z]?([-]?[1-9][0-9]*)?)+( ([1-9][0-9]*)?[A-Z][a-z]?([-]?[1-9][0-9]*)?)*)$/";
string FormValidatingJavascript::elementalFormulaAllowBlankRegEx = "/^$|^((([1-9][0-9]*)?[A-Z][a-z]?([1-9][0-9]*)?)+( ([1-9][0-9]*)?[A-Z][a-z]?([1-9][0-9]*)?)*)$/";

string FormValidatingJavascript::positiveNonZeroIntegerAllowBlankRegEx	= "/^$|^([1-9][0-9]*)$/";
string FormValidatingJavascript::positiveFloatingPointAllowBlankRegEx	= "/^$|^([0-9]*\\.?[0-9]+)$/";
string FormValidatingJavascript::signedFloatingPointAllowBlankRegEx		= "/^$|^([-+]?[0-9]*\\.?[0-9]+)$/";

string FormValidatingJavascript::positiveFloatingPointOrExponentRegEx	= "/^([0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)$/";
string FormValidatingJavascript::signedFloatingPointOrExponentRegEx		= "/^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)$/";

string FormValidatingJavascript::sensibleFilenameRegEx			= "/^([a-zA-Z0-9][\\w #\\-\\+\\.\\(\\)\\[\\]\\{\\}]*)$/";
string FormValidatingJavascript::sensibleFilenameAllowBlankRegEx= "/^$|^([a-zA-Z0-9][\\w #\\-\\+\\.\\(\\)\\[\\]\\{\\}]*)$/";

string FormValidatingJavascript::msPatternRegEx = "/^([A-Z0-9\\[\\]\\{\\}\\.\\*\\\\\\^\\,\\$\\-]+)$/";

string FormValidatingJavascript::multipleSpaces = "/ {2,}/";
string FormValidatingJavascript::trailingSpace = "/ $/";
bool FormValidatingJavascript::validateItemPrinted = false;
bool FormValidatingJavascript::validateTextFunctionPrinted = false;

string FormValidatingJavascript::addNumberValidator ( const string& name, const string& alertName, const string& rxp )
{
	PairStringString pss;
	pss.first = alertName;
	pss.second = rxp;
	numberValidator [name] = pss;
	string jf = "validateItem( this.form.";
	jf += name;
	jf += ", ";
	jf += "'" + alertName + "'";
	jf += ", ";
	jf += rxp;
	jf += " )";
	return jf;
}
string FormValidatingJavascript::addTextValidator ( const string& name, const string& alertName, const string& rxp )
{
	PairStringString pss;
	pss.first = alertName;
	pss.second = rxp;
	textValidator [name] = pss;
	string jf = "validateTextItem( this.form.";
	jf += name;
	jf += ", ";
	jf += "'" + alertName + "'";
	jf += ", ";
	jf += rxp;
	jf += " )";
	return jf;
}
void FormValidatingJavascript::printUtilityFunctions ( ostream& os ) const
{
	if ( !validateItemPrinted ) {
		printValidateItemFunctions ( os );
		validateItemPrinted = true;
	}
	if ( !textValidator.empty () && !validateTextFunctionPrinted ) {
		printValidateTextFunction ( os );
		validateTextFunctionPrinted = true;
	}
}
void FormValidatingJavascript::printValidateItemFunctions ( ostream& os )
{
	os << "function validateItem( item, alertName, rxp ) {" << endl;	//Validates if regular expression true
	os << "\t" << "var val=item.value;" << endl;
	os << "\t" << "if(!rxp.test(val)){" << endl;
	os << "\t\t" << "alert( alertName + \" is invalid\");" << endl;
	os << "\t\t" << "item.value=\"\";" << endl;
	os << "\t\t" << "item.select();" << endl;
	os << "\t\t" << "return false;" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "return true;" << endl;
	os << "}" << endl;
	os << "function validateItem2( item, alertName, rxp ) {" << endl;	//Validates if regular expression false
	os << "\t" << "var val=item.value;" << endl;
	os << "\t" << "if(rxp.test(val)){" << endl;
	os << "\t\t" << "alert( alertName + \" is invalid\");" << endl;
	os << "\t\t" << "item.value=\"\";" << endl;
	os << "\t\t" << "item.select();" << endl;
	os << "\t\t" << "return false;" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "return true;" << endl;
	os << "}" << endl;
}
void FormValidatingJavascript::printValidateTextFunction ( ostream& os )
{
	os << "function validateTextItem( item, alertName, rxp ) {" << endl;
	os << "\t" << "if (validateItem2(item, alertName, " << multipleSpaces << ") && validateItem2(item, alertName, " << trailingSpace << ")) {" << endl;
	os << "\t\t" << "if (validateItem(item, alertName, rxp ) ) {" << endl;
	os << "\t\t\t" << "return true;" << endl;
	os << "\t\t" << "}" << endl;
	os << "\t\t" << "return false;" << endl;
	os << "\t" << "}" << endl;
	os << "\t" << "return false;" << endl;
	os << "}" << endl;
}
void FormValidatingJavascript::print ( ostream& os ) const
{
	startJavascript ( os );
	printUtilityFunctions ( os );
	printValidateFormFunction ( os );
	endJavascript ( os );
}
void FormValidatingJavascript::printValidateFormFunction ( ostream& os ) const
{
	os << "function validateForm( form ) {" << endl;
	for ( MapStringToPairStringStringConstIterator i = textValidator.begin () ; i != textValidator.end () ; i++ ) {
		os << "\t" << "if ( !validateTextItem(form." << (*i).first << ", \"" << (*i).second.first << "\", " << (*i).second.second << " ) ) {" << endl;
		os << "\t\t" << "return false;" << endl;
		os << "\t" << "}" << endl;
	}
	for ( MapStringToPairStringStringConstIterator j = numberValidator.begin () ; j != numberValidator.end () ; j++ ) {
		os << "\t" << "if ( !validateItem(form." << (*j).first << ", \"" << (*j).second.first << "\", " << (*j).second.second << " ) ) {" << endl;
		os << "\t\t" << "return false;" << endl;
		os << "\t" << "}" << endl;
	}
	os << "\t" << "return true;" << endl;
	os << "}" << endl;
}
