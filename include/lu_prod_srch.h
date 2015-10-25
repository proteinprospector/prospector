/******************************************************************************
*                                                                             *
*  Program    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_prod_srch.h                                                *
*                                                                             *
*  Created    : June 18th 2001                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_prod_srch_h
#define __lu_prod_srch_h

#include <lu_prod_par.h>
#include <lu_program.h>

class BiemannFragments;
class PeakContainer;
class PepInfo;

class MSProductSearch : public MSProgram {
	MSProductParameters& prodParams;
	std::vector <PepInfo*> pInfo;
	BiemannFragments* biemannFragments;
	std::vector <PeakContainer*> peaks;
	double precursorMOverZ;
	bool rawDataOption;
	void printRawDataValidationJavascript ( std::ostream& os ) const;
#ifdef RAW_DATA
	void writeSpectrumReport ( std::ostream& os, double mOverZ, const SpecID& specID, const ParameterList* params );
#endif
	void printIsotopeLink ( std::ostream& os, ElementalFormula& elemComp, int maxCharge ) const;
	void printPrecursor ( std::ostream& os, int maxCharge, double peptideMassMi, double peptideMassAv, int precision ) const;
	void printJavascriptFunctions ( std::ostream& os );
	void setLinkSearchVisualizationFlags ( const std::string& val, bool& div_ls );
public:
	MSProductSearch ( MSProductParameters& params );
	~MSProductSearch ();
	void printBodyHTML ( std::ostream& os );
	void printBodyTabDelimitedText ( std::ostream& os );
	void printBodyMGFSpectrum ( std::ostream& os );
	void printForm ( std::ostream& os );
	void printBodyXML ( std::ostream& os );
	void printParamsBodyHTML ( std::ostream& os ) const { prodParams.printHTML ( os ); };
};

#endif /* ! __lu_prod_srch_h */
