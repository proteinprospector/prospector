/******************************************************************************
*                                                                             *
*  Library    : libucsf                                                       *
*                                                                             *
*  Filename   : lu_faindex.h                                                  *
*                                                                             *
*  Created    : October 27th 1997                                             *
*                                                                             *
*  Purpose    : Functions for creating the database index files.              *
*                                                                             *
*  Author(s)  : Peter Baker                                                   *
*                                                                             *
*  This file is the confidential and proprietary product of The Regents of    *
*  the University of California.  Any unauthorized use, reproduction or       *
*  transfer of this file is strictly prohibited.                              *
*                                                                             *
*  Copyright (1997-2009) The Regents of the University of California.         *
*                                                                             *
*  All rights reserved.                                                       *
*                                                                             *
******************************************************************************/

#ifndef __lu_faindex_h
#define __lu_faindex_h

class FastaServer;

void create_new_database_from_indicies_list ( FastaServer* fs, const std::string& database, const std::string& subDatabaseID, const IntVector& indicies );
void create_database_files ( const std::string& file_name, bool parallel );
void add_single_database_entry ( const std::string& filename, const std::string& protein, const std::string& name, const std::string& species, const std::string& accession_number );
std::string dna_database_to_protein_database ( const std::string& database, bool delete_dna_database );
std::string convertDatabaseToRandomDatabase ( const std::string& database, bool reverseFlag, bool concatFlag );

#endif /* ! __lu_faindex_h */
