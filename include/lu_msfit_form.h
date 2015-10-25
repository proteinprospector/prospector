/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_msfit_form.h                                               *
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
*  Copyright (2004-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_msfit_form_h
#define __lu_msfit_form_h

#include <ostream>
#include <lu_srch_form.h>
#include <lu_data_form.h>
#include <lu_form_valid.h>

class MSFitForm : public ProspectorForm {
protected:
	FormValidatingJavascript fvj;
	HeaderForm headerForm;
	SearchForm searchForm;
	MatrixModeForm matrixModeForm;
	MSToleranceForm msToleranceForm;

	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	MSFitForm ( const VectorConstParameterListPtr& params );
	void printHTMLMSFitItems ( std::ostream& os );
	virtual void printHTML ( std::ostream& os ) = 0;
};

class MSFitStandardForm : public MSFitForm {
protected:
	PasteAreaDataForm dataForm;
public:
	MSFitStandardForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

class MSFitUploadForm : public MSFitForm {
protected:
	MSUploadDataForm dataForm;
public:
	MSFitUploadForm ( const VectorConstParameterListPtr& params );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_msfit_form_h */
