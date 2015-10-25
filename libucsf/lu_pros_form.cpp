/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pros_form.cpp                                              *
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
*  Copyright (2004-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_pros_form.h>
#include <lu_html_form.h>
#include <lu_param_list.h>
using std::ostream;
using std::vector;

void ProspectorForm::create ( const VectorConstParameterListPtr& params )
{
	if ( !params.empty () ) {
		const ParameterList* p = params [0];
		pFile = p->getStringValue ( "output_filename" );
	}
	setOptions ();
	createItems ();
	setValues ( params );
}
void ProspectorForm::showHiddenItems ( ostream& os ) const
{
	for ( MapStringFormItemPtrConstIterator i = formItemMap.begin () ; i != formItemMap.end () ; i++ ) {
		if ( i->second->getShown () == false ) i->second->printHTMLHidden ( os );
	}
}
void ProspectorForm::printCGI ( ostream& os ) const
{
	for ( MapStringFormItemPtrConstIterator i = formItemMap.begin () ; i != formItemMap.end () ; i++ ) {
		i->second->printCGI ( os );
	}
}
void ProspectorForm::printHTMLJavascriptHidden ( ostream& os ) const
{
	for ( MapStringFormItemPtrConstIterator i = formItemMap.begin () ; i != formItemMap.end () ; i++ ) {
		i->second->printHTMLJavascriptHidden ( os );
	}
}
