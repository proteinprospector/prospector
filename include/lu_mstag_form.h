/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mstag_form.h                                               *
*                                                                             *
*  Created    : December 7th 2004                                             *
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

#ifndef __lu_mstag_form_h
#define __lu_mstag_form_h

#include <ostream>
#include <lu_srch_form.h>
#include <lu_data_form.h>
#include <lu_form_valid.h>

class MSTagIonTypeForm : public ProspectorForm {
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSTagIonTypeForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSSeqForm : public ProspectorForm {
protected:
	FormValidatingJavascript fvj;
	HeaderForm headerForm;
	SearchForm searchForm;
	MSMSToleranceForm msmsToleranceForm;
	ImmoniumIonForm immoniumIonForm;
	PasteAreaDataForm dataForm;

	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSSeqForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSTagForm : public ProspectorForm {
protected:
	FormValidatingJavascript fvj;
	HeaderForm headerForm;
	SearchForm searchForm;
	MSMSToleranceForm msmsToleranceForm;
	VariableModsForm variableModsForm;
	ExtraUserModsForm extraUserModsForm;
	MassModificationForm massModificationForm;
	CrosslinkingForm crosslinkingForm;
	MatrixModeForm matrixModeForm;
	bool filterPeaks;

	void printHTMLFirstTable ( std::ostream& os );
	void printHTMLSecondTable ( std::ostream& os );
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSTagForm ( const VectorConstParameterListPtr& params, bool filterPeaks );
	virtual void printHTML ( std::ostream& os ) = 0;
};

class MSTagFileForm : public MSTagForm {
protected:
	DataForm dataForm;
public:
	MSTagFileForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSTagFileFormFromViewer : public MSTagForm {
protected:
	DataFormFromViewer dataForm;
public:
	MSTagFileFormFromViewer ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSTagStandardForm : public MSTagForm {
protected:
	ImmoniumIonForm immoniumIonForm;
	MSTagIonTypeForm ionTypeForm;
	PasteAreaDataForm dataForm;
public:
	MSTagStandardForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_mstag_form_h */
