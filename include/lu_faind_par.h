/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_faind_par.h                                                *
*                                                                             *
*  Created    : June 13th 2001                                                *
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

#ifndef __lu_faind_par_h
#define __lu_faind_par_h

#include <lu_mass.h>
#include <lu_prog_par.h>
#include <lu_pre_srch.h>

class FAIndexParameters : public MSProgramParameters {
	AAInitInfo aaInitInfo;
	std::string database;
	bool createSubDatabase;
	bool createUserDatabase;
	bool createDatabaseIndicies;
public:
	FAIndexParameters ( const ParameterList* params );
	std::string getDatabase () const { return database; }
	bool getCreateSubDatabase () const { return createSubDatabase; }
	bool getCreateUserDatabase () const { return createUserDatabase; }
	bool getCreateDatabaseIndicies () const { return createDatabaseIndicies; }
};

class FAIndexNormalParameters : public FAIndexParameters {
	bool dnaToProtein;
	bool randomDatabase;
	bool reverseDatabase;
	bool concatDatabase;
	bool deleteDNADatabase;
public:
	FAIndexNormalParameters ( const ParameterList* params );
	bool getDNAToProtein () const { return dnaToProtein; }
	bool getRandomDatabase () const { return randomDatabase; }
	bool getReverseDatabase () const { return reverseDatabase; }
	bool getConcatDatabase () const { return concatDatabase; }
	bool getDeleteDNADatabase () const { return deleteDNADatabase; }
};

class FAIndexSubsetDatabaseParameters : public FAIndexParameters {
	std::string subDatabaseID;
	PreSearchInfo preSearchInfo;
public:
	FAIndexSubsetDatabaseParameters ( const ParameterList* params );
	std::string getSubDatabaseID () const { return subDatabaseID; }
	IntVector getIndicies ( FastaServer* fsPtr );
};

class FAIndexDatabaseSummaryParameters : public FAIndexParameters {
	int dnaReadingFrame;
	int startIndexNumber;
	int endIndexNumber;
	bool allIndicies;
	bool hideProteinSequence;
public:
	FAIndexDatabaseSummaryParameters ( const ParameterList* params );
	int getDNAReadingFrame () const { return dnaReadingFrame; }
	int getStartIndexNumber () const { return startIndexNumber; }
	int getEndIndexNumber () const { return endIndexNumber; }
	bool getAllIndicies () const { return allIndicies; }
	bool getHideProteinSequence () const { return hideProteinSequence; }
};

class FAIndexCreateOrAppendParameters : public FAIndexParameters {
	std::string userProteinSequence;
	std::string nameField;
	std::string accessionNum;
	std::string species;
public:
	FAIndexCreateOrAppendParameters ( const ParameterList* params );
	std::string getUserProteinSequence () const { return userProteinSequence; }
	std::string getNameField () const { return nameField; }
	std::string getAccessionNum () const { return accessionNum; }
	std::string getSpecies () const { return species; }
};

#endif /* ! __lu_faind_par_h */
