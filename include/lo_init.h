/******************************************************************************
*                                                                             *
*  Library    : liboracle                                                     *
*                                                                             *
*  Filename   : lo_init.h                                                     *
*                                                                             *
*  Created    :                                                               *
*                                                                             *
*  Purpose    :                                                               *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (2002-2007) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/
#ifndef __lo_init_h
#define __lo_init_h

#include <string>
#include <vector>
#include <oci.h>

class OracleConnection {
	OCIEnv* pEnv;		// The OCI environment pointer
    OCISvcCtx* pSvc;	// The OCI service pointer
	OCIError* pErr;		// The OCI error handle
public:
	OracleConnection ( const std::string& server, const std::string& username, const std::string& password );
	~OracleConnection ();
	OCIEnv* get_pEnv () const { return pEnv; }
    OCISvcCtx* get_pSvc () const { return pSvc; }
	OCIError* get_pErr () const { return pErr; }
};

class OracleStatement {
  	OCIEnv* pEnv;	// The OCI environment pointer
	OCISvcCtx* pSvc;
	OCIError* pErr;

	OCIStmt* pSelect;
	OCIStmt* cursor;
	int bindIndex;
	OCIDefine** pDefn;
	OCIBind** pBind;
	int i;
	ub4 nPrefetch;
	std::vector <sb2> pInd;
public:
	OracleStatement ( OracleConnection& oc, const char* pszSQL );
	virtual ~OracleStatement ();
	void BindByName ();
	void BindByName ( const char* placeholder, dvoid* valuep, sb4 value_sz, ub2 dty );
	void BindByPos ();
	void AttrSet ( const int prefetchRows );
	void StmtExecute ();
	sword StmtFetch ( const unsigned int blockSize = 1 );
	void printAttributeList ();
	void DescriptorAlloc ( OCILobLocator** blob );
	void DescriptorFree ( OCILobLocator* blob );
	template <class T>
		void DefineByPosition ( T* value, ub2 dty )
	{
		i++;
		checkerr ( OCIDefineByPos ( cursor, &pDefn[i], pErr, i, value, sizeof(T), dty, (dvoid *) 0, (ub2 *) 0, (ub2 *) 0, OCI_DEFAULT ) );
	}
	void DefineArrayOfStruct ();
	sword GetRowCount ( ub4* nRows );
	void reset () { i = 0; }
	int checkerr ( sword nStatus );
	void writeBlob ( OracleConnection& oc, const std::string& filename, OCILobLocator* blob );
	void close ();
};

#endif /* ! __lo_init_h */
