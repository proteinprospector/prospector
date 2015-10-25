/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_fas_enz.h                                                  *
*                                                                             *
*  Created    : April 9th 2001                                                *
*                                                                             *
*  Purpose    : Enzymatic and chemical digest cleavage functions.             *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2001-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_fas_enz_h
#define __lu_fas_enz_h

#include <string>
#include <vector>
#include <lgen_define.h>
#include <lgen_sort.h>

class ParameterList;

class DigestTable {
	StringVector names;
	MapStringToString breakMaskTable;
	MapStringToString excludeMaskTable;
	MapStringToChar digestSpecificityTable;
	static void errorHandler ( const std::string& enzymeName );
	DigestTable ();
public:
	static DigestTable& instance ();
	std::string getBreakMask ( const std::string& enzymeName ) const;
	std::string getExcludeMask ( const std::string& enzymeName ) const;
	char getSpecificity ( const std::string& enzymeName ) const;
	StringVector getNames () const { return names; }
};

class EnzymeParameters {
	std::string enzyme;
	std::string allowNonSpecific;
	int missedCleavages;
	bool endTerminus;
	char strippingTerminal;
	int startStrip;
	int endStrip;
	StringVector compIons;
	std::string compMaskType;
	unsigned int compMask;
	static std::string initEnzyme ( const ParameterList* params, bool noNoEnzyme );

	static std::string ENZYME;
	static std::string ALLOW_NON_SPECIFIC;
	static std::string MISSED_CLEAVAGES;
	static std::string END_TERMINUS;
	static std::string STRIPPING_TERMINAL;
	static std::string START_STRIP;
	static std::string END_STRIP;
	static std::string COMP_ION;
	static std::string COMP_MASK_TYPE;

	static const std::string NO_ENZYME;
	static std::string ENZYME_DEFAULT;
	static std::string ALLOW_NON_SPECIFIC_DEFAULT;
	static int MISSED_CLEAVAGES_DEFAULT;
	static bool END_TERMINUS_DEFAULT;
	static char STRIPPING_TERMINAL_DEFAULT;
	static int START_STRIP_DEFAULT;
	static int END_STRIP_DEFAULT;
	static std::string COMP_MASK_TYPE_DEFAULT;
public:
	EnzymeParameters ( const ParameterList* params, bool noNoEnzyme = false );
	std::string getEnzyme () const { return enzyme; }
	bool getNoEnzyme () const { return enzyme == NO_ENZYME; }
	std::string getAllowNonSpecific () const { return allowNonSpecific; }
	int getMissedCleavages () const { return missedCleavages; }
	bool getEndTerminus () const { return endTerminus; }
	char getStrippingTerminal () const { return strippingTerminal; }
	int getStartStrip () const { return startStrip; }
	int getEndStrip () const { return endStrip; }
	unsigned int getCompMask () const { return compMask; }
	std::string getCompMaskType () const { return compMaskType; }
	void printHTML ( std::ostream& os ) const;
	void putCGI ( std::ostream& os ) const;
	void putNoEnzymeCGI ( std::ostream& os ) const;
	void putHiddenFormEntryJavascript ( std::ostream& os ) const;
	void putNoEnzymeHiddenFormEntryJavascript ( std::ostream& os ) const;

	static void copyToCGI ( std::ostream& os, const ParameterList* params, bool noNoEnzyme = false );
	static void copyToHiddenFormEntry ( std::ostream& os, const ParameterList* params, bool noNoEnzyme = false );
	static void copyToHiddenFormJavascriptEntry ( std::ostream& os, const ParameterList* params, bool noNoEnzyme = false );
};

void init_fasta_enzyme_function ( const std::string& enzyme );
char get_enzyme_terminal_specificity ();
DoubleVector& get_cleaved_masses ( const std::string& protein, const IntVector& cleavage_index );
DoubleVector& get_cleaved_masses_to_limit ( const std::string& protein, const IntVector& cleavage_index, double limit );

#ifdef LUCSF_FAS_ENZ_MAIN
#define LUCSF_FAS_ENZ_EXTERN
#else
#define LUCSF_FAS_ENZ_EXTERN extern
#endif
LUCSF_FAS_ENZ_EXTERN IntVector& ( *enzyme_fragmenter ) ( const std::string& peptideFormula );
LUCSF_FAS_ENZ_EXTERN bool cnbr_digest;

#endif /* ! __lu_fas_enz_h */
