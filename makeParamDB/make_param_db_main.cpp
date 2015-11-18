/******************************************************************************
*                                                                             *
*  Program    : makeParamDB                                                   *
*                                                                             *
*  Filename   : make_param_db_main.cpp                                        *
*                                                                             *
*  Created    : November 5th 2004                                             *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2014-2014) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#include <lu_param_db.h>
int main ( int argc, char** argv )
{
	ParamDBWrite pdbw ( "test.sqlite" );
	pdbw.create ();
	return 1;
}
