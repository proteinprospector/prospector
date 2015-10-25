/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pspot_form.h                                               *
*                                                                             *
*  Created    : December 15th 2004                                            *
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

#ifndef __lu_pspot_form_h
#define __lu_pspot_form_h

#include <ostream>
#include <lu_html_form.h>
#include <lu_pros_form.h>

class PeakSpotterForm : public ProspectorForm {
protected:
	StringVector spot_set_names_options;
	StringVector intensity_type_options;
	StringVector instrument_name_options;

	void setOptions ();
	void setValues ( const VectorConstParameterListPtr& params );
	void createItems ();
public:
	PeakSpotterForm ( const VectorConstParameterListPtr& params, const StringVector& spotSetNamesOptions );
	virtual void printHTML ( std::ostream& os );
};

#endif /* ! __lu_pspot_form_h */
