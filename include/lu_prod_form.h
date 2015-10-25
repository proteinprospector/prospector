/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prod_form.h                                                *
*                                                                             *
*  Created    : January 10th 2005                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_prod_form_h
#define __lu_prod_form_h

#include <ostream>
#include <lu_srch_form.h>
#include <lu_data_form.h>
#include <lu_form_valid.h>

class MSIsotopeForm : public ProspectorForm {
	FormValidatingJavascript fvj;
	UserAAForm userAAForm;
	SaveResultsForm saveResultsForm;
	int numSequences;
	static void printJavascriptFunctions ( std::ostream& os );
protected:
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSIsotopeForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class IonTypeForm : public ProspectorForm {
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	IonTypeForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
	void removeName ( ParameterList* p );
};

class MSProductForm : public ProspectorForm {
	FormValidatingJavascript fvj;
	IonTypeForm ionTypeForm;
	UserAAForm userAAForm;
	SaveResultsForm saveResultsForm;
	ProspectorForm* dataForm;
	int numSequences;
protected:
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSProductForm ( const VectorConstParameterListPtr& params, ProspectorForm* dataForm );
	virtual void printHTML ( std::ostream& os );
};

class MSProductPasteForm : public MSProductForm {
public:
	MSProductPasteForm ( const VectorConstParameterListPtr& params );
};

class MSProductUploadForm : public MSProductForm {
public:
	MSProductUploadForm ( const VectorConstParameterListPtr& params );
};

class FormItemMaxCharge : public FormItemSelect {
	static const char* options [];
public:
	FormItemMaxCharge () :
		FormItemSelect ( "Max Charge", "prodman.htm#Max_charge", getName (), options, "1" ) {}
	static std::string getName () { return "max_charge"; }
};

class FormItemCountPosZ : public FormItemSelect {
	static const char* options [];
public:
	FormItemCountPosZ () :
		FormItemSelect ( "", "", getName (), options, "Count Basic AA" ) {}
	static std::string getName () { return "count_pos_z"; }
};

class FormItemMaxInternalLen : public FormItemText {
public:
	FormItemMaxInternalLen ( FormValidatingJavascript* fvj );
	static std::string getName () { return "max_internal_len"; }
};

#endif /* ! __lu_prod_form_h */
