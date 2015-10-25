/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_scomp_form.h                                               *
*                                                                             *
*  Created    : October 29th 2004                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_scomp_form_h
#define __lu_scomp_form_h

#include <ostream>
#include <lu_form_valid.h>
#include <lu_srch_form.h>

class FormItemSCFormat : public FormItemSelect {
	static const char* options1 [];
	static const char* options2 [];
public:
	FormItemSCFormat ( bool showPepXML = false ) :
		FormItemSelect ( "Format", "", getName (), showPepXML ? options2 : options1, "HTML" ) {}
	static std::string getName () { return "save_format"; }
};

class FormItemSCCheckboxes : public FormItemCheckbox {
public:
	FormItemSCCheckboxes () :
		FormItemCheckbox ( "Checkboxes", "", getName (), false ) {}
	static std::string getName () { return "report_checkboxes"; }
};

class FormItemSCReportType : public FormItemSelect {
	static const char* optionsCal [];
	static const char* optionsNormal [];
	static const char* optionsXLink [];
	static const char** getOptions ( bool cal, bool xlink );
	static std::string getSelect ( bool cal, bool xlink, const std::string& type );
public:
	FormItemSCReportType ( bool cal, bool xlink, const std::string& type = "" );
	static std::string getName () { return "report_type"; }
};


class SearchCompareForm : public ProspectorForm {
	PresentCompIonForm compIonForm;
	bool calForm;
	int numSearches;
	void printRawDataValidationJavascript ( std::ostream& os ) const;
protected:
	bool expectationFlag;
	bool xLinkFlag;
	StringVector instrumentList;
	FormValidatingJavascript fvj;
	CrosslinkingForm crosslinkingForm;
	void setOptions () {}
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
	void printHTMLDiscriminentItems ( std::ostream& os );
public:
	SearchCompareForm ( const VectorConstParameterListPtr& params, const std::string& reportDefault = "Protein" );
	virtual void printHTML ( std::ostream& os );
};

class SearchCompareFormCalibrationForm : public SearchCompareForm {
public:
	SearchCompareFormCalibrationForm ( const VectorConstParameterListPtr& params );
	void printHTML ( std::ostream& os );
};

#endif /* ! __lu_scomp_form_h */
