/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_program.h                                                  *
*                                                                             *
*  Created    : October 22nd 2001                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2012) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_program_h
#define __lu_program_h

#include <ostream>
#include <lu_prog_par.h>

class MSProgram {
	const MSProgramParameters& params;
	virtual void printBodyHTML ( std::ostream& os ) {};
	virtual void saveHits () const {};
	void printParamsHTML ( std::ostream& os ) const;
	virtual void printParamsBodyHTML ( std::ostream& os ) const {};
	static const ParameterList* paramList;
public:
	MSProgram ( const MSProgramParameters& params );
	virtual ~MSProgram ();
	void outputResults ( bool showParams = true );
	void printHTML ( std::ostream& os, bool showParams );
	void printXMLTop ( std::ostream& os );
	virtual void printBodyTabDelimitedText ( std::ostream& os ) {};
	virtual void printBodyMGFSpectrum ( std::ostream& os ) {};
	virtual void printBodyXML ( std::ostream& os ) {};
	void printXMLBottom ( std::ostream& os );
	void printXML ( std::ostream& os );
	void printTabDelimitedText ( std::ostream& os );
	void printMGFSpectrum ( std::ostream& os );
	static void setParams ( const ParameterList* p ) { paramList = p; };
};

#endif /* ! __lu_program_h */
