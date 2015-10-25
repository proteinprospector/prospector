/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_acc_num.h                                                  *
*                                                                             *
*  Created    : April 4th 2001                                                *
*                                                                             *
*  Purpose    : Functions for doing accession number searches.                *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_acc_num_h
#define __lu_acc_num_h

#include <string>

class AccessionNumberMap {
public:
	virtual ~AccessionNumberMap () = 0;
	virtual int getIndexNumber ( const std::string& accessionNumber ) const = 0;
	virtual bool accessionNumberUnique ( const std::string& accessionNumber ) const = 0;
	virtual void reset () const = 0;
};

AccessionNumberMap* getAccessionNumberMap ( const std::string& database );

#endif /* ! __lu_acc_num_h */
