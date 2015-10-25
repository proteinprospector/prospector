/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_r_plot.h                                                   *
*                                                                             *
*  Created    : March 8th 2006                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2006-2010) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_r_plot_h
#define __lu_r_plot_h

#include <ostream>
#include <string>

class RPlot {
	static std::string rCommand;
	static bool rFlag;
	static bool keepRDataFile;
	std::string rScriptName;
	std::string dataFileFullPath;
	static std::string getRCommand ();
public:
	RPlot ( const std::string& rScriptName );
	virtual ~RPlot ();
	std::string getDataFileFullPath () const { return dataFileFullPath; }
	void printImage ( std::ostream& os, const std::string& imageType ) const;
	void printImageAndLink ( std::ostream& os ) const;
	static bool getRFlag () { return rFlag; }
};

#endif /* ! __lu_r_plot_h */
