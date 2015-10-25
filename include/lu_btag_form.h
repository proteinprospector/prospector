/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_btag_form.h                                                *
*                                                                             *
*  Created    : December 13th 2004                                            *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_btag_form_h
#define __lu_btag_form_h

#include <ostream>
#include <lu_srch_form.h>
#include <lu_data_form.h>
#include <lu_form_valid.h>

class BatchTagForm : public ProspectorForm {
protected:
	FormValidatingJavascript fvj;
	HeaderForm headerForm;
	SearchForm searchForm;
	VariableModsForm variableModsForm;
	ExtraUserModsForm extraUserModsForm;
	MassModificationForm massModificationForm;
	CrosslinkingForm crosslinkingForm;
	MatrixModeForm matrixModeForm;
	MSMSToleranceForm msmsToleranceForm;
	std::string formType;

	virtual void setValues ( const VectorConstParameterListPtr& params );
	virtual void createItems ();
	virtual void printHTMLFormInit ( std::ostream& os ) const = 0;
	virtual void printDataForm ( std::ostream& os ) = 0;
public:
	BatchTagForm ( const VectorConstParameterListPtr& params, const std::string& formType );
	virtual void printHTML ( std::ostream& os );
};

class BatchTagFullForm : public BatchTagForm {
	DataForm dataForm;
	virtual void printHTMLFormInit ( std::ostream& os ) const;
	virtual void printDataForm ( std::ostream& os );
public:
	BatchTagFullForm ( const VectorConstParameterListPtr& params, const std::string& formType );
};

class BatchTagWebForm : public BatchTagForm {
	UploadDataForm uploadDataForm;
	std::string user;
	bool setCookie;
	virtual void printHTMLFormInit ( std::ostream& os ) const;
	virtual void printDataForm ( std::ostream& os );
public:
	BatchTagWebForm ( const VectorConstParameterListPtr& params, bool setCookie );
};

#endif /* ! __lu_btag_form_h */
