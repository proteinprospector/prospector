/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_delim.h                                                    *
*                                                                             *
*  Created    : April 1st 2004                                                *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2004-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_delim_h
#define __lu_delim_h

#include <ostream>
#include <string>

void delimitedRowStart ( std::ostream& os );
void delimitedRowEnd ( std::ostream& os );
void delimitedHeader ( std::ostream& os, const std::string& str );
void delimitedCellSigFig ( std::ostream& os, double number, int sigFig );
void delimitedCell ( std::ostream& os, double number, int precision );
void delimitedCellRange ( std::ostream& os, double number1, double number2, int precision );
void delimitedCell ( std::ostream& os, const std::string& str );
void delimitedCell ( std::ostream& os, char number );
void delimitedCell ( std::ostream& os, int number );
void delimitedEmptyCell ( std::ostream& os );
void delimitedEmptyNCells ( std::ostream& os, int num );

#endif /* ! __lu_delim_h */
