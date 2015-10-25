/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_mass.h                                                     *
*                                                                             *
*  Created    : July 20th 1996                                                *
*                                                                             *
*  Purpose    : Include file for:                                             *
*                                                                             *
*               lu_mass_pep.cpp                                               *
*                                                                             *
*               Initialises some global arrays concerned with mass            *
*               calculations.                                                 *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2011) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_mass_h
#define __lu_mass_h

#include <map>
#include <string>
class ConstMod;
typedef std::map <std::string, ConstMod*> MapStringConstModPtr;				// Avoids declaration of const_mod.h
typedef std::map <char, ConstMod*> MapCharConstModPtr;						// Avoids declaration of const_mod.h
typedef MapCharConstModPtr::const_iterator MapCharConstModPtrConstIterator;	// Avoids declaration of const_mod.h
#include <lu_formula.h>

typedef int MassType;
const double MASS_TYPE_MULTIPLIER = 10000.0;	//Should be 1 for floating point, 10000 for int
//typedef double MassType;
//const double MASS_TYPE_MULTIPLIER = 1.0;	//Should be 1 for floating point, 10000 for int
typedef std::vector <MassType> MassTypeVector;
typedef std::vector <MassTypeVector> MassTypeVectorVector;

const int AA_ARRAY_SIZE = 128;
const double MONOISOTOPIC_NOMINAL_MASS_OFFSET = 1.00050792;
const double AVERAGE_NOMINAL_MASS_OFFSET = 1.00114130;
const double AVERAGE_TO_MONOISOTOPIC_FACTOR = 0.999367342;
const double MONOISOTOPIC_TO_AVERAGE_FACTOR = 1.000633058;
const double ELECTRON_REST_MASS = 0.00054858;

/******************************************************************************
*               lu_mass_pep.cpp                                               *
******************************************************************************/

class ModCodeMap {
	static std::string oneLetterCodes;
	static std::string lowerCaseCodes;
	static const char* threeLetterCodes [];
	static const char* lowerCasePSIMods [];
	MapCharToString codeMap;
	MapCharToString codeMap2;
	MapStringToChar codeMap3;
	ModCodeMap ();
public:
	static ModCodeMap& instance ();
	std::string getModString ( char db, char mut );
	std::string getLowerMod ( char aa ) { return codeMap2 [aa]; }
	bool getMutationFromPSIMod ( const std::string& psiMod, char& mutation );
};

void checkConstantAndVariableModCompatibility ( const ParameterList* params );

class AAInitInfo {
	StringVector constModNames;
	std::string nTermName;
	std::string cTermName;
	std::string nLossName;
	MapStringConstModPtr constMods;
	MapCharConstModPtr cm;
	std::string cation;
	std::string instrumentName;
	std::vector <ElementalFormulaInfo> userAAElemForms;
public:
	AAInitInfo ( const ParameterList* params );
	void printHTML ( std::ostream& os ) const;
	void putCGI ( std::ostream& os ) const;
	void putHiddenFormEntryJavascript ( std::ostream& os ) const;

	ElementalFormula getNTermModFormula () const;
	ElementalFormula getCTermModFormula () const;
	ElementalFormulaVector getUserAAElemForms () const
	{
		ElementalFormulaVector elemForm;
		for ( ElementalFormulaVectorSizeType i = 0 ; i < userAAElemForms.size () ; i++ ) {
			elemForm.push_back ( userAAElemForms [i].getElementalFormula () );
		}
		return elemForm;
	}
	MapCharToString getConstModMap () const;
	MapStringConstModPtr getConstMods () const { return constMods; }
	ConstMod* getConstMod ( char c ) const
	{
		MapCharConstModPtrConstIterator cur = cm.find ( c );
		if ( cur != cm.end () )
			return (*cur).second;
		else
			return 0;
	}
	std::string getCation () const { return cation; }
	std::string getInstrumentName () const { return instrumentName; }

	static void copyToCGI ( std::ostream& os, const ParameterList* params, bool cysUnmodified = false );
	static void copyToHiddenFormEntry ( std::ostream& os, const ParameterList* params, bool cysUnmodified = false );
	static void copyToHiddenFormJavascriptEntry ( std::ostream& os, const ParameterList* params, bool cysUnmodified = false );
	std::string getNTermName () const { return nTermName; }
	std::string getCTermName () const { return cTermName; }
	std::string getNLossName () const { return nLossName; }
};

class MassInfo {
	std::string massType;
	bool monoisotopicFlag;

	static std::string MASS_TYPE;

	static std::string MASS_TYPE_DEFAULT;
	static bool MONOISOTOPIC_FLAG_DEFAULT;
public:
	MassInfo ( const ParameterList* params );
	void printHTML ( std::ostream& os ) const;
	void putCGI ( std::ostream& os ) const;
	void putCGIFragment ( std::ostream& os ) const;
	void putHiddenFormEntry ( std::ostream& os ) const;
	static void copyToCGI ( std::ostream& os, const ParameterList* params );
	static void copyToHiddenFormEntry ( std::ostream& os, const ParameterList* params );
	static void copyToHiddenFormJavascriptEntry ( std::ostream& os, const ParameterList* params );

	std::string getMassType () const { return massType; }
	bool getMonoisotopicFlag () const { return monoisotopicFlag; }
	bool getMonoParentAverageFragments () const { return massType == "Par(mi)Frag(av)"; }
	bool getAverageParentMonoFragments () const { return massType == "Par(av)Frag(mi)"; }
};

void initialise_amino_acid_weights ( const MapStringConstModPtr& constMods, const std::vector <ElementalFormula>& userAAElemFormulae, bool monoisotopicFlag );
void modified_mass_convert ( bool monoisotopicFlag );
unsigned int string_to_mask ( const std::string& str );
unsigned int string_to_mask ( const char* str, int len );
double getNTermMassShift ( const std::string& modAA );
double getCTermMassShift ( const std::string& modAA );
double getNeutralLossMassShift ( const std::string& modAA );
void modifyAAInfo ( char aaCode, const std::string& modAA );
double getAminoAcidWt ( const std::string& modAA );
int sumBasicResidues ( const std::string& fragment );
int sumBasicResidues ( const StringVector& fragment );
int sumResidues ( const std::string& fragment, unsigned int mask );
double getMaxDSubstituent ();
double getMaxWSubstituent ();

/* Exported data from the lu_mass_pep module */

#ifdef LUCSF_MASS_PEP_MAIN
#define LUCSF_MASS_PEP_EXTERN
#else
#define LUCSF_MASS_PEP_EXTERN extern
#endif

LUCSF_MASS_PEP_EXTERN double *amino_acid_wt;
LUCSF_MASS_PEP_EXTERN double *substituent_da;
LUCSF_MASS_PEP_EXTERN double *substituent_wa;
LUCSF_MASS_PEP_EXTERN double *substituent_db;
LUCSF_MASS_PEP_EXTERN double *substituent_wb;
LUCSF_MASS_PEP_EXTERN MassTypeVector aaArrayMT;
LUCSF_MASS_PEP_EXTERN MassTypeVector substituentDaMT;
LUCSF_MASS_PEP_EXTERN MassTypeVector substituentDbMT;
LUCSF_MASS_PEP_EXTERN MassTypeVector substituentWaMT;
LUCSF_MASS_PEP_EXTERN MassTypeVector substituentWbMT;
LUCSF_MASS_PEP_EXTERN double pk_c_term [AA_ARRAY_SIZE];
LUCSF_MASS_PEP_EXTERN double pk_n_term [AA_ARRAY_SIZE];
LUCSF_MASS_PEP_EXTERN double pk_acidic_sc [AA_ARRAY_SIZE];
LUCSF_MASS_PEP_EXTERN double pk_basic_sc [AA_ARRAY_SIZE];
LUCSF_MASS_PEP_EXTERN double ( *massConvert ) (const char*);
LUCSF_MASS_PEP_EXTERN double n_terminus_wt;
LUCSF_MASS_PEP_EXTERN double c_terminus_wt;
LUCSF_MASS_PEP_EXTERN double cation_wt;
LUCSF_MASS_PEP_EXTERN double terminal_wt;
LUCSF_MASS_PEP_EXTERN double internal_n_terminus_wt;
LUCSF_MASS_PEP_EXTERN char* terminal_formula;

LUCSF_MASS_PEP_EXTERN unsigned int aa_composition_mask [AA_ARRAY_SIZE];

LUCSF_MASS_PEP_EXTERN double cnbr_homoserine_lactone_mod;

LUCSF_MASS_PEP_EXTERN double c2_h2_n_o;
LUCSF_MASS_PEP_EXTERN double s_o_c_h4;
LUCSF_MASS_PEP_EXTERN double h3_p_o4;
LUCSF_MASS_PEP_EXTERN double o_h;
LUCSF_MASS_PEP_EXTERN double h2_o;
LUCSF_MASS_PEP_EXTERN double n_h3;
LUCSF_MASS_PEP_EXTERN double h1;
LUCSF_MASS_PEP_EXTERN double h2;
LUCSF_MASS_PEP_EXTERN double c_o;

LUCSF_MASS_PEP_EXTERN unsigned int phosphorylation_mask;
LUCSF_MASS_PEP_EXTERN unsigned int oxidized_m_mask;

LUCSF_MASS_PEP_EXTERN double a_tag_offset;
LUCSF_MASS_PEP_EXTERN double a_h2o_tag_offset;
LUCSF_MASS_PEP_EXTERN double a_nh3_tag_offset;
LUCSF_MASS_PEP_EXTERN double a_h3po4_tag_offset;
LUCSF_MASS_PEP_EXTERN double b_tag_offset;
LUCSF_MASS_PEP_EXTERN double b_h2o_tag_offset;
LUCSF_MASS_PEP_EXTERN double b_plus_h2o_tag_offset;
LUCSF_MASS_PEP_EXTERN double b_nh3_tag_offset;
LUCSF_MASS_PEP_EXTERN double b_soch4_tag_offset;
LUCSF_MASS_PEP_EXTERN double b_h3po4_tag_offset;
LUCSF_MASS_PEP_EXTERN double cPlus2DaTagOffset;
LUCSF_MASS_PEP_EXTERN double cPlus1DaTagOffset;
LUCSF_MASS_PEP_EXTERN double c_tag_offset;
LUCSF_MASS_PEP_EXTERN double cMinus1DaTagOffset;
LUCSF_MASS_PEP_EXTERN double d_tag_offset;

LUCSF_MASS_PEP_EXTERN double v_tag_offset;
LUCSF_MASS_PEP_EXTERN double w_tag_offset;
LUCSF_MASS_PEP_EXTERN double x_tag_offset;
LUCSF_MASS_PEP_EXTERN double y_tag_offset;
LUCSF_MASS_PEP_EXTERN double y_nh3_tag_offset;
LUCSF_MASS_PEP_EXTERN double y_h2o_tag_offset;
LUCSF_MASS_PEP_EXTERN double y_soch4_tag_offset;
LUCSF_MASS_PEP_EXTERN double y_h3po4_tag_offset;
LUCSF_MASS_PEP_EXTERN double Y_tag_offset;
LUCSF_MASS_PEP_EXTERN double z_tag_offset;
LUCSF_MASS_PEP_EXTERN double zPlus1DaTagOffset;
LUCSF_MASS_PEP_EXTERN double zPlus2DaTagOffset;
LUCSF_MASS_PEP_EXTERN double zPlus3DaTagOffset;

LUCSF_MASS_PEP_EXTERN double nominal_mass_offset;

int* get_protein_int_mass_array ( const char* protein );
double* get_protein_double_mass_array ( const char* protein );

#endif /* ! __lu_mass_h */
