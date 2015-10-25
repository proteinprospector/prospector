/******************************************************************************
*                                                                             *
*  Library    : liboracle                                                     *
*                                                                             *
*  Filename   : lo_ss_list.h                                                  *
*                                                                             *
*  Created    :                                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef __lo_ss_list_h
#define __lo_ss_list_h

#include <map>
#include <string>
#include <ostream>
#include <lgen_define.h>
#include <lo_init.h>

class SpotSetList : public OracleStatement {
	MapStringToUInt spotSets;
public:
	SpotSetList ( OracleConnection& oc  );
	virtual ~SpotSetList ();
	void print ( std::ostream& os ) const;
	unsigned int getID ( const std::string& spotSet ) const;
	StringVector getList () const;
};


#endif /* ! __lo_ss_list_h */
