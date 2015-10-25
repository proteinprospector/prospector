/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_count_scan.h                                               *
*                                                                             *
*  Created    : August 3rd 2005                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2005-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_count_scan_h
#define __lu_count_scan_h

#include <string>

class CountScans {
protected:
	int count;
public:
	virtual ~CountScans ();
	virtual int getCount () const { return count; }
};

class PPExpatCountScans;

class XMLCountScans : public CountScans {
	PPExpatCountScans* ppesc;
public:
	XMLCountScans ( PPExpatCountScans* ppesc, const std::string& filename );
	~XMLCountScans ();
	int getCount () const;
};

class MascotCountScans : public CountScans {
public:
	MascotCountScans ( const std::string& filename );
};

class MS2CountScans : public CountScans {
public:
	MS2CountScans ( const std::string& filename );
};

class APLCountScans : public CountScans {
public:
	APLCountScans ( const std::string& filename );
};

class PPCountMSScans : public CountScans {
public:
	PPCountMSScans ( const std::string& filename );
};

class PPCountMSMSScans : public CountScans {
public:
	PPCountMSMSScans ( const std::string& filename );
};

class SpaceSeparatedSpectraCountScans : public CountScans {
public:
	SpaceSeparatedSpectraCountScans ( const std::string& filename );
};

#endif /* ! __lu_count_scan_h */
