/******************************************************************************
*                                                                             *
*  Library    : liboracle                                                     *
*                                                                             *
*  Filename   : lo_jr_view.h                                                  *
*                                                                             *
*  Created    : March 30th 2006                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2006-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef __lo_jr_view_h
#define __lo_jr_view_h

#include <lgen_define.h>
#include <lo_init.h>

class SpotSetJobRuns : public OracleStatement {
	UIntVector jobRunIDList;
public:
	SpotSetJobRuns ( OracleConnection& oc, unsigned int spotSetID );
	virtual ~SpotSetJobRuns ();
	UIntVector getJobRunIDList () const { return jobRunIDList; }
};


#endif /* ! __lo_jr_view_h */
