/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_data_form.h                                                *
*                                                                             *
*  Created    : December 9th 2004                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_data_form_h
#define __lu_data_form_h

#include <ostream>
#include <lg_string.h>
#include <lu_html_form.h>
#include <lu_pros_form.h>

class FormItemInstrumentName : public FormItemSelect {
public:
	FormItemInstrumentName ( int num = 1 );
	static std::string getName ( int num = 1 )
	{
		std::string n = ( num == 1 ) ? "" : gen_itoa ( num );
		return "instrument_name" + n;
	}
};

class DataFormFromViewer : public ProspectorForm {
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	DataFormFromViewer ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};
class DataForm : public ProspectorForm {
	std::string searchName;
	bool filterValue;

	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	DataForm ( const VectorConstParameterListPtr& params, const std::string& searchName );
	virtual void printHTML ( std::ostream& os );
};

class UploadDataForm : public ProspectorForm {
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	UploadDataForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSUploadDataForm : public ProspectorForm {
	std::string searchName;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSUploadDataForm ( const VectorConstParameterListPtr& params, const std::string& searchName );
	virtual void printHTML ( std::ostream& os );
};

class PasteAreaDataForm : public ProspectorForm {
	std::string searchNameValue;

	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	PasteAreaDataForm ( const VectorConstParameterListPtr& params, const std::string& searchNameValue );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_data_form_h */
