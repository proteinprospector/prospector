##################################################################################
#                                                                                #
#  Library    : libucsf                                                          #
#                                                                                #
#  Filename   : libucsf.linux.basic.make                                         #
#                                                                                #
#  Created    : September 11th 2001                                              #
#                                                                                #
#  Purpose    : LINUX makefile for libucsf.                                      #
#                                                                                #
#  Author(s)  : Peter Baker                                                      #
#                                                                                #
#  This file is the confidential and proprietary product of The Regents of       #
#  the University of California.  Any unauthorized use, reproduction or          #
#  transfer of this file is strictly prohibited.                                 #
#                                                                                #
#  Copyright (2001-2015) The Regents of the University of California.            #
#                                                                                #
#  All rights reserved.                                                          #
#                                                                                #
##################################################################################

LIBRARY_NAME = ../lib/libucsf.a
COMPILER=g++
ARCHIVER=ar
RANLIB=ranlib
OPTIONS= -O2 -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
ADD_OPTIONS=
INCLUDEDIRS=-I../include

ARCHIVE_COMMAND = $(ARCHIVER) r $(LIBRARY_NAME)

OBJ=lu_aa_calc.o \
	lu_acc_link.o \
	lu_acc_num.o \
	lu_add_user.o \
	lu_ambiguity.o \
	lu_app_gr.o \
	lu_bdg_par.o \
	lu_bdg_srch.o \
	lu_blib.o \
	lu_brid_form.o \
	lu_btag_form.o \
	lu_btag_run.o \
	lu_btag_submit.o \
	lu_btsel_form.o \
	lu_cgi_val.o \
	lu_charge.o \
	lu_check_db.o \
	lu_chem_sc.o \
	lu_comb_perm.o \
	lu_comp_form.o \
	lu_comp_par.o \
	lu_comp_rgex.o \
 	lu_comp_srch.o \
	lu_composit.o \
	lu_const_mod.o \
	lu_cookie.o \
	lu_count_scan.o \
	lu_coverage.o \
	lu_data.o \
	lu_data_form.o \
	lu_db_entry.o \
	lu_db_srch.o \
	lu_dbst_form.o \
	lu_dbst_srch.o \
	lu_del_proj.o \
	lu_delim.o \
	lu_df_info.o \
	lu_dig_par.o \
	lu_dig_srch.o \
	lu_disc_sc.o \
	lu_distribution.o \
	lu_dlsel_form.o \
	lu_elem_comp.o \
	lu_expec_par.o \
	lu_export_proj.o \
	lu_faind_form.o \
	lu_faind_par.o \
	lu_faindex.o \
	lu_fas_db.o \
	lu_fas_enz.o \
	lu_fas_ind.o \
	lu_fas_sp.o \
	lu_file_split.o \
	lu_file_type.o \
	lu_filter_form.o \
	lu_filter_par.o \
	lu_filter_srch.o \
	lu_fit_par.o \
	lu_fit_srch.o \
	lu_form_valid.o \
	lu_formula.o \
	lu_frag_info.o \
	lu_frag_mtch.o \
	lu_fragmentation.o \
	lu_get_file.o \
	lu_get_imm.o \
	lu_get_link.o \
	lu_histogram.o \
	lu_hit.o \
	lu_hom_form.o \
	lu_hom_par.o \
	lu_hom_srch.o \
	lu_html_form.o \
	lu_html_out.o \
	lu_import_proj.o \
	lu_indicies.o \
	lu_inst.o \
	lu_iso_dist.o \
	lu_iso_par.o \
	lu_iso_srch.o \
	lu_links.o \
	lu_login_form.o \
	lu_mask.o \
	lu_mass_aa.o \
	lu_mass_conv.o \
	lu_mass_elem.o \
	lu_mass_frag.o \
	lu_mass_pep.o \
	lu_mass_seq.o \
	lu_mat_score.o \
	lu_mgf.o \
	lu_msp.o \
	lu_mod_frag.o \
	lu_msf.o \
	lu_msfit_form.o \
	lu_mstag_form.o \
	lu_mut_mtrx.o \
	lu_mzidentml.o \
	lu_nspec_par.o \
	lu_nspec_srch.o \
	lu_param_db.o \
	lu_param_list.o \
	lu_parent.o \
	lu_patt_form.o \
	lu_patt_par.o \
	lu_patt_srch.o \
	lu_pep_xml.o \
	lu_pi.o \
	lu_pk_filter.o \
	lu_pp_param.o \
	lu_pre_srch.o \
	lu_prod_form.o \
	lu_prod_par.o \
	lu_prod_srch.o \
	lu_prog.o \
	lu_prog_par.o \
	lu_program.o \
	lu_proj_file.o \
	lu_proj_form.o \
	lu_pros_form.o \
	lu_prot_srch.o \
	lu_pspot_form.o \
	lu_quan_multi.o \
	lu_quan_ratio.o \
	lu_r_plot.o \
	lu_rep_links.o \
	lu_repos_info.o \
	lu_repository.o \
	lu_scomp_form.o \
	lu_script.o \
	lu_scsel_form.o \
	lu_sctag_link.o \
	lu_seq_exp.o \
	lu_sim_ent.o \
	lu_sing_fit.o \
	lu_spec_id.o \
	lu_species.o \
	lu_spep_srch.o \
	lu_sqlite.o \
	lu_srch_form.o \
	lu_srch_par.o \
	lu_t2d.o \
	lu_table.o \
	lu_tag_frag.o \
	lu_tag_par.o \
	lu_tag_seq.o \
	lu_tag_srch.o \
	lu_tofil_par.o \
	lu_tol.o \
	lu_unk_rexp.o \
	lu_usermod.o \
	lu_version.o \
	lu_viewer_form.o \
	lu_viewer_par.o \
	lu_viewer_srch.o \
	lu_xml.o \
	lu_xml_data.o

.SUFFIXES:	.cpp

.cpp.o:	
	$(COMPILER) $(OPTIONS) $(ADD_OPTIONS) $(INCLUDEDIRS) -c $<

libucsf.a:$(OBJ)
	$(ARCHIVE_COMMAND) $(OBJ)
	$(RANLIB) $(LIBRARY_NAME)
clean:	
	/bin/rm -f *.o 	
