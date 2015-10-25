/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prog_par.h                                                 *
*                                                                             *
*  Created    : September 10th 2001                                           *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_prog_par_h
#define __lu_prog_par_h

#include <string>
#include <lu_tofil_par.h>

class ParameterList;

class MSProgramParameters {
private:
	std::string searchName;
	std::string reportTitle;
	std::string outputType;
	ToFileParameters toFileParameters;
	bool displayGraph;
public:
	MSProgramParameters ( const ParameterList* params );
	virtual ~MSProgramParameters ();

	std::string getSearchName () const { return searchName; }
	std::string getReportTitle () const { return reportTitle; }
	std::string getOutputType () const { return outputType; }
	ToFileParameters getToFileParameters () const { return toFileParameters; }
	bool getToFile () const { return toFileParameters.getToFile (); }
	std::string getOutputFilename () const { return toFileParameters.getOutputFilename (); }
	bool isRandomSearch () const { return toFileParameters.isRandomSearch (); }
	std::string getOutputPath ( const std::string& suffix ) const { return toFileParameters.getOutputPath ( suffix ); }
	bool getDisplayGraph () const { return displayGraph; }
	void printOutputFileLink ( std::ostream& os, const std::string& linkText, const std::string& fileSuffix ) const
		{ toFileParameters.printOutputFileLink ( os, linkText, fileSuffix ); }
};
void initialiseProspector ();

#endif /* ! __lu_prog_par_h */
