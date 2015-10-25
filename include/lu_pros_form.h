/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_pros_form.h                                                *
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

#ifndef __lu_pros_form_h
#define __lu_pros_form_h

#include <ostream>
#include <map>
#include <vector>

class ParameterList;
class FormItem;

typedef std::vector <const ParameterList*> VectorConstParameterListPtr;
typedef VectorConstParameterListPtr::size_type VectorConstParameterListPtrSizeType;

class ProspectorForm {
protected:
	std::string pFile;
	typedef std::map <std::string, FormItem*> MapStringFormItemPtr;
	typedef MapStringFormItemPtr::const_iterator MapStringFormItemPtrConstIterator;
	MapStringFormItemPtr formItemMap;
	virtual void setOptions () {};
	virtual void setValues ( const std::vector <const ParameterList*>& params ) {}
	virtual void createItems () = 0;
	void create ( const VectorConstParameterListPtr& params );
public:
	ProspectorForm () {}
	virtual void printHTML ( std::ostream& os ) = 0;
	void showHiddenItems ( std::ostream& os ) const;
	virtual void printCGI ( std::ostream& os ) const;
	virtual void printHTMLJavascriptHidden ( std::ostream& os ) const;
};

#endif /* ! __lu_pros_form_h */
