/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_xml.h                                                      *
*                                                                             *
*  Created    : February 3rd 2003                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_xml_h
#define __lu_xml_h

#include <iostream>
#include <string>
#include <lgen_define.h>

void printXMLHeader ( std::ostream& os );
void printXMLStylesheet ( std::ostream& os, const std::string& scriptType, const std::string& scriptName );
void printXMLVersion ( std::ostream& os );

std::string getXMLFileDate ( const std::string& fname );

void openXMLTag ( std::ostream& os, int ntab, const std::string& name, const VectorPairStringString& attr = VectorPairStringString (), bool close = false );
void openXMLValueTag ( std::ostream& os, int ntab, const std::string& name, const VectorPairStringString& attr = VectorPairStringString (), bool close = false );
void closeXMLTag ( std::ostream& os, int ntab, const std::string& name );
void closeXMLTag ( std::ostream& os, const std::string& name );

namespace XMLParser {
std::string getStringValue ( const std::string& str, const std::string& tag );
std::string getStringValue ( const std::string& str, const std::string& tag, const std::string& defaultValue );
bool getBoolValue ( const std::string& str, const std::string& tag );
char getCharValue ( const std::string& str, const std::string& tag );
int getIntValue ( const std::string& str, const std::string& tag );
double getDoubleValue ( const std::string& str, const std::string& tag );
double getDoubleValue ( const std::string& str, const std::string& tag, const std::string& defaultValue );
double getDoubleValue ( const std::string& str, const std::string& tag, double defaultValue );
void getMZI ( const std::string& str, double& mOverZ, int& charge, double& intensity );
StringVector getStringVectorValue ( const std::string& str, const std::string& tag );
}

class XMLStringList {
	const std::string& str;
	std::string nameStart;
	int len;
	std::string nameEnd;
	int start;
public:
	XMLStringList ( const std::string& str, const std::string& tag );
	bool getNext ( std::string& s );
};

class XMLIStreamList {
	std::istream& ist;
	std::string nameStart;
	std::string nameEnd;
public:
	XMLIStreamList ( std::istream& ist, const std::string& tag );
	bool getNext ( std::string& s );
	bool getNextBlock ( std::string& s );
};

class XMLOutputItem {
protected:
	std::string tagName;
	VectorPairStringString attr;
public:
	XMLOutputItem ( const std::string& tagName, const std::string& id ) :
		tagName ( tagName )
	{
		if ( !id.empty () ) attr.push_back ( makePairStringString ( "id", id ) );
	}
	virtual void print ( std::ostream& os, int ntab ) const = 0;
	void printOpenTag ( std::ostream& os, int ntab ) const;
	void printCloseTag ( std::ostream& os, int ntab ) const;
};
typedef std::vector <XMLOutputItem*> VectorXMLOutputItemPtr;
typedef VectorXMLOutputItemPtr::size_type VectorXMLOutputItemPtrSizeType;

/*

Example XMLOutputAttrItem

<enzymatic_search_constraint enzyme="trypsin_k" max_num_internal_cleavages="2" min_number_termini="1"/>

*/

class XMLOutputAttrItem : public XMLOutputItem {
public:
	XMLOutputAttrItem ( const std::string& tagName, const std::string& id = "" ) :
		XMLOutputItem ( tagName, id ) {}
	void print ( std::ostream& os, int ntab ) const;
};

/* 
Example XMLOutputContainerItem

<Provider id="PROVIDER">

......

</Provider>

Attributes can be added by subclass

*/ 
class XMLOutputContainerItem : public XMLOutputItem {
protected:
	VectorXMLOutputItemPtr subItems;
public:
	XMLOutputContainerItem ( const std::string& tagName, const VectorXMLOutputItemPtr& subItems = VectorXMLOutputItemPtr (), const std::string& id = "" ) :
		XMLOutputItem ( tagName, id ),
		subItems ( subItems ) {}
	~XMLOutputContainerItem ();
	void print ( std::ostream& os, int ntab ) const;
	void addSubItem ( XMLOutputItem* si )
	{
		subItems.push_back ( si );
	}
	void updateSubItems ( const VectorXMLOutputItemPtr& sit )
	{
		subItems = sit;
	}
};

/*

Example XMLOutputValueItem

<seq>MGLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHGTVV</seq>

*/

class XMLOutputValueItem : public XMLOutputItem {
protected:
	std::string value;
public:
	XMLOutputValueItem ( const std::string& tagName, const std::string& id, const std::string& value ) :
		XMLOutputItem ( tagName, id ),
		value ( value ) {}
	XMLOutputValueItem ( const std::string& tagName, const std::string& value ) :
		XMLOutputItem ( tagName, "" ),
		value ( value ) {}
	void print ( std::ostream& os, int ntab ) const;
};

#endif /* ! __lu_xml_h */
