/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_form_valid.h                                               *
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

#ifndef __lu_form_valid_h
#define __lu_form_valid_h

#include <lgen_define.h>

class FormValidatingJavascript {
	MapStringToPairStringString numberValidator;
	MapStringToPairStringString textValidator;
	static bool validateItemPrinted;
	static bool validateTextFunctionPrinted;

	static std::string positiveNonZeroIntegerRegEx;
	static std::string positiveIntegerRegEx;
	//static std::string positiveDecimalRegEx;
	//static std::string positiveIntegerOrDecimalRegEx;
	//static std::string signedIntegerOrDecimalRegEx;
	static std::string positiveFloatingPointRegEx;
	static std::string signedFloatingPointRegEx;
	static std::string positiveFloatingPointOrExponentRegEx;
	static std::string signedFloatingPointOrExponentRegEx;

	static std::string elementalModificationFormulaAllowBlankRegEx;
	static std::string elementalFormulaAllowBlankRegEx;

	static std::string positiveNonZeroIntegerAllowBlankRegEx;
	static std::string positiveFloatingPointAllowBlankRegEx;
	static std::string signedFloatingPointAllowBlankRegEx;

	static std::string sensibleFilenameRegEx;
	static std::string sensibleFilenameAllowBlankRegEx;
	static std::string msPatternRegEx;
	static std::string multipleSpaces;
	static std::string trailingSpace;
	void printUtilityFunctions ( std::ostream& os ) const;
	static void printValidateItemFunctions ( std::ostream& os );
	static void printValidateTextFunction ( std::ostream& os );
	std::string addNumberValidator ( const std::string& name, const std::string& alertName, const std::string& rxp );
	std::string addTextValidator ( const std::string& name, const std::string& alertName, const std::string& rxp );
	void printValidateFormFunction ( std::ostream& os ) const;
public:
	FormValidatingJavascript () {}

	std::string addPositiveNonZeroIntegerValidator ( const std::string& name, const std::string& alertName )
		{ return addNumberValidator ( name, alertName, positiveNonZeroIntegerRegEx ); }
	std::string addPositiveIntegerValidator ( const std::string& name, const std::string& alertName )
		{ return addNumberValidator ( name, alertName, positiveIntegerRegEx ); }
	std::string addPositiveFloatingPointValidator ( const std::string& name, const std::string& alertName )
		{ return addNumberValidator ( name, alertName, positiveFloatingPointRegEx ); }
	std::string addSignedFloatingPointValidator ( const std::string& name, const std::string& alertName )
		{ return addNumberValidator ( name, alertName, signedFloatingPointRegEx ); }
	std::string addPositiveFloatingPointOrExponentValidator ( const std::string& name, const std::string& alertName )
		{ return addNumberValidator ( name, alertName, positiveFloatingPointOrExponentRegEx ); }
	std::string addSignedFloatingPointOrExponentValidator ( const std::string& name, const std::string& alertName )
		{ return addNumberValidator ( name, alertName, signedFloatingPointOrExponentRegEx ); }

	std::string addElementalModificationFormulaAllowBlankValidator ( const std::string& name, const std::string& alertName )
		{ return addTextValidator ( name, alertName, elementalModificationFormulaAllowBlankRegEx ); }
	std::string addElementalFormulaAllowBlankValidator ( const std::string& name, const std::string& alertName )
		{ return addTextValidator ( name, alertName, elementalFormulaAllowBlankRegEx ); }

	std::string addPositiveNonZeroIntegerAllowBlankValidator ( const std::string& name, const std::string& alertName )
		{ return addNumberValidator ( name, alertName, positiveNonZeroIntegerAllowBlankRegEx ); }
	std::string addPositiveFloatingPointAllowBlankValidator ( const std::string& name, const std::string& alertName )
		{ return addNumberValidator ( name, alertName, positiveFloatingPointAllowBlankRegEx ); }
	std::string addSignedFloatingPointAllowBlankValidator ( const std::string& name, const std::string& alertName )
		{ return addNumberValidator ( name, alertName, signedFloatingPointAllowBlankRegEx ); }

	std::string addFilenameValidator ( const std::string& name, const std::string& alertName )
		{ return addTextValidator ( name, alertName, sensibleFilenameRegEx ); }
	std::string addFilenameAllowBlankValidator ( const std::string& name, const std::string& alertName )
		{ return addTextValidator ( name, alertName, sensibleFilenameAllowBlankRegEx ); }
	std::string addMSPatternValidator ( const std::string& name, const std::string& alertName )
		{ return addTextValidator ( name, alertName, msPatternRegEx ); }

	void print ( std::ostream& os ) const;
};

#endif /* ! __lu_form_valid_h */
