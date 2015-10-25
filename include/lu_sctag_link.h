/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_sctag_link.h                                               *
*                                                                             *
*  Created    : September 9th 2005                                            *
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

#ifndef __lu_sctag_link_h
#define __lu_sctag_link_h

class SpecID;
class SCMSTagLink {
	void putCGI ( std::ostream& os ) const;
public:
	SCMSTagLink () {}
	void write ( std::ostream& os, const std::string& searchKey, const SpecID& specID, const std::string& str ) const;
	void printHTML ( std::ostream& os ) const;
};

#endif /* ! __lu_sctag_link_h */
