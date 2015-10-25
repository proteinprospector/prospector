/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mgf.h                                                      *
*                                                                             *
*  Created    : April 20th 2007                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2007-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_mgf_h
#define __lu_mgf_h

#include <lgen_define.h>

class MGFInstance;

typedef std::vector <std::pair <std::string, MGFInstance*> > VectorPairStringMGFInstancePtr;
typedef VectorPairStringMGFInstancePtr::size_type VectorPairStringMGFInstancePtrSizeType;

class MGFInfo {
	std::vector <std::pair <std::string, MGFInstance*> > mgfi;
	MGFInstance* currentInstancePtr;
	MGFInfo ();
	MGFInstance* getMGFInstancePtr ( const std::string& title );
public:
	~MGFInfo ();
	static MGFInfo& instance ();
	bool getTitleParams ( const std::string& line, std::string& spot, int& run, std::string& msmsInfo );
	std::string getMSMSInfoFromScans ( const IntVector& scans ) const;
	bool isSpottingPlate ( const std::string& line, bool& spottingPlate );
	void reset ();
};

#endif /* ! __lu_mgf_h */
