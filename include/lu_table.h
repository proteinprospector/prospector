/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_table.h                                                    *
*                                                                             *
*  Created    : March 25th 2003                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2003-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_table_h
#define __lu_table_h

#include <ostream>
#include <string>

void tableStart ( std::ostream& os, bool border = false, const std::string& align = "", const std::string& cellspacing = "", const std::string& id = "" );
void tableEnd ( std::ostream& os );
void tableRowStart ( std::ostream& os );
void tableRowStartDiv ( std::ostream& os, const std::string& id, bool open );
void tableRowEnd ( std::ostream& os );
void tableEmptyRow ( std::ostream& os );
void tableEmptyNRows ( std::ostream& os, int n );
void tableHeader ( std::ostream& os, const std::string& str, const std::string& styleID = "", const std::string& align = "", bool nowrap = false, int colspan = 0, int rowspan = 0 );
void tableHeaderStart ( std::ostream& os, const std::string& styleID = "", const std::string& align = "", bool nowrap = false, int colspan = 0, int rowspan = 0 );
void tableHeaderEnd ( std::ostream& os );
void tableData ( std::ostream& os, const std::string& str, const std::string& styleID = "" );
void tableDataStart ( std::ostream& os, const std::string& styleID = "", const std::string& align = "", bool nowrap = false, int colspan = 0, int rowspan = 0 );
void tableDataEnd ( std::ostream& os );
void tableCellStart ( std::ostream& os, const std::string& styleID = "", const std::string& align = "", bool nowrap = false, int colspan = 0, int rowspan = 0  );
void tableCellEnd ( std::ostream& os );
void tableCellSigFig ( std::ostream& os, double number, int sigFig, bool bold = false, const std::string& styleID = "", int colspan = 0, int rowspan = 0 );
void tableCell ( std::ostream& os, double number, int precision, bool bold = false, const std::string& styleID = "", int colspan = 0, int rowspan = 0 );
void tableCellRange ( std::ostream& os, double number1, double number2, int precision, bool bold = false, const std::string& styleID = "" );
void tableCell ( std::ostream& os, const std::string& str, bool nobr = false, bool bold = false, const std::string& styleID = "", int colspan = 0, int rowspan = 0 );
void tableCell ( std::ostream& os, char number, bool bold = false, const std::string& styleID = "", int colspan = 0, int rowspan = 0 );
void tableCell ( std::ostream& os, int number, bool bold = false, const std::string& styleID = "", int colspan = 0, int rowspan = 0 );
void tableCell ( std::ostream& os, unsigned int number, bool bold = false, const std::string& styleID = "", int colspan = 0, int rowspan = 0 );
void tableEmptyCell ( std::ostream& os, const std::string& styleID = "", int colspan = 0, int rowspan = 0 );
void tableEmptyNCells ( std::ostream& os, int n, const std::string& styleID = "" );
void divStart ( std::ostream& os, const std::string& name, bool open );
void divEnd ( std::ostream& os );
void spanStart ( std::ostream& os, const std::string& name, bool open );
void spanEnd ( std::ostream& os );

#endif /* ! __lu_table_h */
