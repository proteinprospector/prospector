/******************************************************************************
*                                                                             *
*  Library    : libgen                                                        *
*                                                                             *
*  Filename   : lgen_reg_exp.h                                                *
*                                                                             *
*  Created    : June 20th 1996                                                *
*                                                                             *
*  Purpose    : Regular expression search functions.                          *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1996-2013) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lgen_reg_exp_h
#define __lgen_reg_exp_h

#include <string>

class RegularExpression {
protected:
	static const int MAX_BUF;
	static const int CBRA;
	static const int CCHR;
	static const int CDOT;
	static const int CCL;
	static const int CXCL;
	static const int CDOL;
	static const int CCEOF;
	static const int CKET;
	static const int CBRC;
	static const int CLET;
	static const int CBACK;
	static const int NCCL;

	static const int STAR;
	static const int RNGE;

	static const int NBRA;

	static const char bittab [];

	char* expbuf;
	int circf;
	const char* loc1;
	const char* loc2;
	const char* locs;
	const char* braslist[9/*NBRA*/];
	const char* braelist[9/*NBRA*/];
	int	size;
	int low;
	bool sameString;
	const char* multiStr;

	int nerr;

	char* compile ( const char* instr, char* ep, char* endbuf, int seof );
	virtual int step ( const char* p1, const char* p2 );
	virtual int advance ( const char* lp, const char* ep );
	void getrnge ( const char* str );
	void reg_err ( int error );
public:
	RegularExpression ( const std::string& regexp, bool caseInsensitive = false );
	~RegularExpression ();
	bool isPresent ( const char* str );
	bool isPresentMulti ( const char* str );
	bool isPresentMultiOverlap ( const char* str );
	std::string getMatch () const;
	const char* getMatchEndPointer () const;
	char getPrevious () const;
	char getNext () const;
	short getStartAA () const;
	const char* getLoc1 () const { return loc1; }
	int getNErr () const { return nerr; }
	void resetMulti () { sameString = false; }
};


class RegularExpressionWithErrors : public RegularExpression {
	int maxerr;
	int step ( const char* p1, const char* p2 );
	int advance ( const char* lp, const char* ep );
public:
	RegularExpressionWithErrors ( const std::string& regexp, int maxErrors );
	const char* errorLocation ( const char* string, bool reset );
};

#endif /* ! __lgen_reg_exp_h */
