/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_html_form.h                                                *
*                                                                             *
*  Created    : September 16th 2002                                           *
*                                                                             *
*  Purpose    : Functions for printing out commonly used elements of HTML     *
*               forms.                                                        *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2015) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_html_form_h
#define __lu_html_form_h

#include <ostream>
#include <string>
#include <vector>
#include <lgen_define.h>
#include <lu_cgi_val.h>

class ParameterList;

void printHTMLFORMLabel ( std::ostream& os, const std::string& label, const std::string& manual );
void printHTMLFORMStart ( std::ostream& os, const std::string& method, const std::string& program, bool multipart = false, bool onSubmit = false, bool newWindow = false, const std::string& type = ".cgi" );
void printHTMLFORMEnd ( std::ostream& os );
void printHTMLFORMSubmit ( std::ostream& os, const std::string& label, const std::string& color = "009F00" );
void printHTMLFORMCheckbox ( std::ostream& os, const std::string& label, const std::string& manual, const std::string& name, bool checked = false, const std::string& value = "1", const std::string& onChangeFunction = "" );
void printHTMLFORMSelect ( std::ostream& os, const std::string& label, const std::string& manual, const std::string& name, const StringVector& value, const std::string& select, int size = 1, const std::string& onChangeFunction = "", const StringVector& theClass = StringVector () );
void printHTMLFORMSelectOptgroup ( std::ostream& os, const std::string& label, const std::string& manual, const std::string& name, const StringVector& value, const StringVector& optgroup, const std::string& select, int size );
void printHTMLFORMSelectMultiple ( std::ostream& os, const std::string& label, const std::string& manual, const std::string& name, const StringVector& value, const StringVector& select, int size = 1, const std::string& onChangeFunction = "", const StringVector& theClass = StringVector ());
void printHTMLFORMPassword ( std::ostream& os, const std::string& label, const std::string& manual, const std::string& name, int size, int maxlength, const std::string& value );
void printHTMLFooter ( std::ostream& os, const std::string& date );
void printHTMLCopyright ( std::ostream& os, const std::string& startYear );
void printHTMLFORMHidden ( std::ostream& os, const std::string& name, double value, int precision );
void printHTMLFORMHiddenSigFig ( std::ostream& os, const std::string& name, double value, int precision );
template <class T>
void printHTMLFORMHidden ( std::ostream& os, const std::string& name, const T& value )
{
	os << "<input";
	os << " ";
	os << "type=\"hidden\"";
	os << " ";
	os << "name=\"" << name << "\"";
	os << " ";
	os << "value=\"" << value << "\"";
	os << " />" << std::endl;
}
void printHTMLFORMJavascriptHidden2 ( std::ostream& os, const std::string& name, const std::string& value );
template <class T>
void printHTMLFORMJavascriptHidden ( std::ostream& os, const std::string& name, const T& value )
{
	os << "\tdocument.writeln ( \"";
	os << "<input";
	os << " ";
	os << "type=\\\"hidden\\\"";
	os << " ";
	os << "name=\\\"" << name << "\\\"";
	os << " ";
	os << "value=\\\"" << value << "\\\"";
	os << " />";
	os << "\" );" << std::endl;;
}
void printHTMLFORMHiddenContainer ( std::ostream& os, const std::string& name, const StringVector& value );
void printHTMLFORMJavascriptHiddenContainer ( std::ostream& os, const std::string& name, const StringVector& value );
void printHTMLFORMJavascriptHiddenContainer2 ( std::ostream& os, const std::string& name, const StringVector& value );
template <class T>
void printHTMLFORMHiddenVector ( std::ostream& os, const std::string& name, const T& value )
{
	os << "<input";
	os << " ";
	os << "type=\"hidden\"";
	os << " ";
	os << "name=\"" << name << "\"";
	os << " ";
	os << "value=\"";
	for ( typename T::size_type i = 0 ; i < value.size () ; i++ ) {
		os << value [i];
		os << "\n";
	}
	os << "\" />" << std::endl;
}
void printHTMLFORMCGI ( std::ostream& os, const std::string& name, const std::string& value );

class FormItem {
protected:
	std::string label;
	std::string manual;
	std::string _name;
	StringVector value;
	mutable bool shown;
public:
	FormItem ( const std::string& label, const std::string& manual, const std::string& name, const StringVector& value );
	virtual void printHTML ( std::ostream& os ) const = 0;
	bool getShown () const { return shown; }
	virtual void printHTMLHidden ( std::ostream& os ) const = 0;
	virtual void printCGI ( std::ostream& os ) const = 0;
	virtual void printHTMLJavascriptHidden ( std::ostream& os ) const = 0;
	virtual void setValue ( const ParameterList* p, const std::string& n = "" ) {}
	static StringVector getOptions ( const char** options )
	{
		StringVector sv;
		for ( int i = 0 ; options [i] != 0 ; i++ ) {
			sv.push_back ( options [i] );
		}
		return sv;
	}
	virtual void setValue ( const StringVector& sv ) { value = sv; }
	virtual void setValue ( bool c ) {}
	virtual StringVector getValue () { return value; }
	virtual std::string getValue ( int i ) { return value [i]; }
	virtual bool getChecked () const { return false; }
	virtual void setOptions ( const StringVector& sv ) {}
	virtual void setOptions ( const ParameterList* p, const std::string& n = "" ) {}
};

class FormItemCheckbox : public FormItem {
	bool checked;
	std::string changeFunction;
public:
	FormItemCheckbox ( const std::string& label, const std::string& manual, const std::string& name, bool checked, const std::string& val = "1", const std::string& changeFunction = "" );
	void printHTML ( std::ostream& os ) const;
	void printHTMLHidden ( std::ostream& os ) const;
	void printCGI ( std::ostream& os ) const;
	void printHTMLJavascriptHidden ( std::ostream& os ) const;
	void setValue ( const ParameterList* p, const std::string& n = "" );
	void setValue ( bool c );
	virtual bool getChecked () const { return checked; }
};

class FormItemText : public FormItem {
	int size;
	int maxlength;
	std::string validateFunction;
public:
	FormItemText ( const std::string& label, const std::string& manual, const std::string& name, int size, int maxlength, const std::string& val, const std::string& validateFunction = "" );
	void printHTML ( std::ostream& os ) const;
	void printHTMLHidden ( std::ostream& os ) const;
	void printCGI ( std::ostream& os ) const;
	void printHTMLJavascriptHidden ( std::ostream& os ) const;
	void setValue ( const ParameterList* p, const std::string& n = "" );
};

class FormItemPassword : public FormItem {
	int size;
	int maxlength;
public:
	FormItemPassword ( const std::string& label, const std::string& manual, const std::string& name, int size, int maxlength, const std::string& val );
	void printHTML ( std::ostream& os ) const;
	void printHTMLHidden ( std::ostream& os ) const;
	void printCGI ( std::ostream& os ) const;
	void printHTMLJavascriptHidden ( std::ostream& os ) const;
	void setValue ( const ParameterList* p, const std::string& n = "" );
};

class FormItemFile : public FormItem {
	int size;
	int maxlength;
public:
	FormItemFile ( const std::string& label, const std::string& manual, const std::string& name, int size, int maxlength, const std::string& val );
	void printHTML ( std::ostream& os ) const;
	void printHTMLHidden ( std::ostream& os ) const;
	void printCGI ( std::ostream& os ) const;
	void printHTMLJavascriptHidden ( std::ostream& os ) const;
	void setValue ( const ParameterList* p, const std::string& n = "" );
};

class FormItemMultiFile : public FormItem {
public:
	FormItemMultiFile ( const std::string& name );
	void printHTML ( std::ostream& os ) const;
	void printHTMLHidden ( std::ostream& os ) const;
	void printCGI ( std::ostream& os ) const;
	void printHTMLJavascriptHidden ( std::ostream& os ) const;
};

class FormItemSelect : public FormItem {
protected:
	StringVector optgroup;
	std::string select;
	StringVector theClass;
	int size;
	std::string changeFunction;
public:
	FormItemSelect ( const std::string& label, const std::string& manual, const std::string& name, const StringVector& value, const StringVector& optgroup, const std::string& select, int size = 1, const std::string& changeFunction = "" );
	FormItemSelect ( const std::string& label, const std::string& manual, const std::string& name, const StringVector& value, const std::string& select, int size = 1, const std::string& changeFunction = "", const StringVector& theClass = StringVector () );
	FormItemSelect ( const std::string& label, const std::string& manual, const std::string& name, const char** options, const std::string& select, int size = 1, const std::string& changeFunction = "" );
	void printHTML ( std::ostream& os ) const;
	void printHTMLHidden ( std::ostream& os ) const;
	void printCGI ( std::ostream& os ) const;
	void printHTMLJavascriptHidden ( std::ostream& os ) const;
	void setValue ( const ParameterList* p, const std::string& n = "" );
	void setValue ( const StringVector& sv );
	virtual std::string getValue ( int i ) { return select; }
	void setOptions ( const StringVector& sv );
};

class FormItemSelectMultiple : public FormItem {
protected:
	StringVector select;
	StringVector theClass;
	int size;
	std::string changeFunction;
public:
	FormItemSelectMultiple ( const std::string& label, const std::string& manual, const std::string& name, const StringVector& value, const StringVector& select, int size = 1, const std::string& changeFunction = "" );
	FormItemSelectMultiple ( const std::string& label, const std::string& manual, const std::string& name, const StringVector& value, const StringVector& select, int size, const StringVector& theClass );
	FormItemSelectMultiple ( const std::string& label, const std::string& manual, const std::string& name, const char** options, const StringVector& select, int size = 1 );
	void printHTML ( std::ostream& os ) const;
	void printHTMLHidden ( std::ostream& os ) const;
	void printCGI ( std::ostream& os ) const;
	void printHTMLJavascriptHidden ( std::ostream& os ) const;
	void setValue ( const ParameterList* p, const std::string& n = "" );
	void setValue ( const StringVector& sv );
	void setOptions ( const ParameterList* p, const std::string& n = "" );
};

class FormItemTextArea : public FormItem {
protected:
	int size;
	int maxlength;
public:
	FormItemTextArea ( const std::string& label, const std::string& manual, const std::string& name, int size, int maxlength, const StringVector& value );
	FormItemTextArea ( const std::string& label, const std::string& manual, const std::string& name, int size, int maxlength, const char** v );
	virtual void printHTML ( std::ostream& os ) const;
	void printHTMLHidden ( std::ostream& os ) const;
	void printCGI ( std::ostream& os ) const;
	void printHTMLJavascriptHidden ( std::ostream& os ) const;
	void setValue ( const ParameterList* p, const std::string& n = "" );
};

class FormItemTextAreaSimple : public FormItemTextArea {
public:
	FormItemTextAreaSimple ( const std::string& label, const std::string& manual, const std::string& name, int size, int maxlength, const StringVector& value );
	FormItemTextAreaSimple ( const std::string& label, const std::string& manual, const std::string& name, int size, int maxlength, const char** v );
	void printHTML ( std::ostream& os ) const;
};

#endif /* ! __lu_html_form_h */
