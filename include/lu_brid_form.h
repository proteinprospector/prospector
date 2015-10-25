/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_brid_form.h                                                *
*                                                                             *
*  Created    : January 7th 2005                                              *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_brid_form_h
#define __lu_brid_form_h

#include <ostream>
#include <lu_data_form.h>
#include <lu_form_valid.h>

class MSSingleForm : public ProspectorForm {
	FormValidatingJavascript* fvj;
	bool enzymeOption;
	SaveResultsForm saveResultsForm;
	PresentCompIonForm* compIonForm;
	static std::string defaultUserSequence;
protected:
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSSingleForm ( const VectorConstParameterListPtr& params, bool enzymeOption, bool tabDelimitedTextOption, FormValidatingJavascript* fvj );
	~MSSingleForm ();
	virtual void printHTML ( std::ostream& os );
	static const std::string getDefaultUserSequence ();
};

class MSBridgeForm : public ProspectorForm {
protected:
	FormValidatingJavascript fvj;
	MSSingleForm singleForm;
	UserAAForm userAAForm;
	VariableModsForm variableModsForm;
	MSToleranceForm msToleranceForm;
	CrosslinkingForm crosslinkingForm;
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSBridgeForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSBridgeFileForm : public MSBridgeForm {
protected:
	DataForm dataForm;
public:
	MSBridgeFileForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSBridgeStandardForm : public MSBridgeForm {
protected:
	PasteAreaDataForm dataForm;
public:
	MSBridgeStandardForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSDigestForm : public ProspectorForm {
	FormValidatingJavascript fvj;
	MSSingleForm singleForm;
	UserAAForm userAAForm;
	VariableModsForm variableModsForm;
protected:
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSDigestForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSNonSpecificForm : public ProspectorForm {
	FormValidatingJavascript fvj;
	MSSingleForm singleForm;
	UserAAForm userAAForm;
	PasteAreaDataForm dataForm;
	MSToleranceForm msToleranceForm;
protected:
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSNonSpecificForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};
#endif /* ! __lu_brid_form_h */
