#
# linux (for mpich2 set MPI_COMPILER=g++, MPILIB=mpich2)
#
STATIC=
ADD_OPTIONS=
C_COMPILER=gcc
COMPILER=g++
MPI_COMPILER=mpic++
ARCHIVER=ar
RANLIB=ranlib
PSILIB=
LLIB=-ldl
MPILIB=
LIBDIRS="-L../lib -L/usr/lib64/mysql"
#

#
# mingw on Windows
#
# STATIC=-static-libgcc
# ADD_OPTIONS=-DVIS_C
# C_COMPILER=x86_64-w64-mingw32-gcc
# COMPILER=x86_64-w64-mingw32-g++
# MPI_COMPILER=x86_64-w64-mingw32-g++
# ARCHIVER=x86_64-w64-mingw32-ar
# RANLIB=x86_64-w64-mingw32-ranlib
# PSILIB=-lpsapi
# LLIB=
# MPILIB=-lmpi
# LIBDIRS="-L../lib"
#

L_COMP_OPTIONS=ARCHIVER=$(ARCHIVER) RANLIB=$(RANLIB) ADD_OPTIONS=$(ADD_OPTIONS)
C_LIB_COMP_OPTIONS=COMPILER=$(C_COMPILER) $(L_COMP_OPTIONS)
LIB_COMP_OPTIONS=COMPILER=$(COMPILER) $(L_COMP_OPTIONS)
MPI_LIB_COMP_OPTIONS=COMPILER=$(MPI_COMPILER) $(L_COMP_OPTIONS)

A_COMP_OPTIONS=LIBDIRS=$(LIBDIRS) STATIC=$(STATIC) ADD_OPTIONS=$(ADD_OPTIONS)
APP_COMP_OPTIONS=COMPILER=$(COMPILER) $(A_COMP_OPTIONS)
MPI_APP_COMP_OPTIONS=COMPILER=$(MPI_COMPILER) $(A_COMP_OPTIONS)

all:
	make -C libzip -f libzip.linux.make $(C_LIB_COMP_OPTIONS)
	make -C libsqlite -f libsqlite.linux.make $(C_LIB_COMP_OPTIONS)
	make -C libexpat -f libexpat.linux.make $(LIB_COMP_OPTIONS)
	make -C libdbase -f libdbase.linux.make $(LIB_COMP_OPTIONS)
	make -C libgen -f libgen.linux.make $(LIB_COMP_OPTIONS)
	make -C libnrec -f libnrec.linux.make $(LIB_COMP_OPTIONS)
	make -C libprocess -f libsingle.linux.full.make $(LIB_COMP_OPTIONS)
	rm -rf libprocess/*.o
	make -C libprocess -f libmulti.linux.full.make $(MPI_LIB_COMP_OPTIONS)
	make -C libraw -f libraw.linux.make $(LIB_COMP_OPTIONS)
	make -C libucsf -f libucsf.linux.full.make $(LIB_COMP_OPTIONS)

	make -C btag-daemon -f btag_daemon.linux.make $(APP_COMP_OPTIONS)
	make -C faindex -f faindex.linux.full.make $(APP_COMP_OPTIONS)
	make -C jobStatus -f jobStatus.linux.make $(APP_COMP_OPTIONS)
	make -C kill -f kill.linux.full.make $(APP_COMP_OPTIONS) PSILIB=$(PSILIB)
	make -C makeParamDB -f make_param_db.linux.full.make $(APP_COMP_OPTIONS)
	make -C msbatch -f msbatch.linux.make $(APP_COMP_OPTIONS)
	make -C msdisplay -f msdisplay.linux.full.make $(APP_COMP_OPTIONS)
	make -C msform -f msform.linux.full.make $(APP_COMP_OPTIONS)
	make -C mssearch -f mssearch.linux.full.make $(APP_COMP_OPTIONS) LLIB=$(LLIB)
	rm -rf mssearch/*.o
	make -C mssearch -f mssearchmpi.linux.full.make $(MPI_APP_COMP_OPTIONS) MPILIB=$(MPILIB) LLIB=$(LLIB)
	make -C readDB -f readDB.linux.make $(APP_COMP_OPTIONS)
	make -C searchCompare -f search_compare.linux.full.make $(APP_COMP_OPTIONS) LLIB=$(LLIB)

allbasic:
	make -C libzip -f libzip.linux.make $(C_LIB_COMP_OPTIONS)
	make -C libsqlite -f libsqlite.linux.make $(C_LIB_COMP_OPTIONS)
	make -C libexpat -f libexpat.linux.make $(LIB_COMP_OPTIONS)
	make -C libgen -f libgen.linux.make $(LIB_COMP_OPTIONS)
	make -C libnrec -f libnrec.linux.make $(LIB_COMP_OPTIONS)
	make -C libprocess -f libsingle.linux.make $(LIB_COMP_OPTIONS)
	make -C libucsf -f libucsf.linux.basic.make $(LIB_COMP_OPTIONS)

	make -C faindex -f faindex.linux.basic.make $(APP_COMP_OPTIONS)
	make -C kill -f kill.linux.basic.make $(APP_COMP_OPTIONS) PSILIB=$(PSILIB)
	make -C makeParamDB -f make_param_db.linux.basic.make $(APP_COMP_OPTIONS)
	make -C msform -f msform.linux.basic.make $(APP_COMP_OPTIONS)
	make -C mssearch -f mssearch.linux.basic.make $(APP_COMP_OPTIONS) LLIB=$(LLIB)

allcl:
	make -C libzip -f libzip.linux.make $(C_LIB_COMP_OPTIONS)
	make -C libsqlite -f libsqlite.linux.make $(C_LIB_COMP_OPTIONS)
	make -C libexpat -f libexpat.linux.make $(LIB_COMP_OPTIONS)
	make -C libgen -f libgen.linux.make $(LIB_COMP_OPTIONS)
	make -C libnrec -f libnrec.linux.make $(LIB_COMP_OPTIONS)
	make -C libprocess -f libsingle.linux.make $(LIB_COMP_OPTIONS)
	make -C libraw -f libraw.linux.make $(LIB_COMP_OPTIONS)
	make -C libucsf -f libucsf.linux.cl.make $(LIB_COMP_OPTIONS)

	make -C faindex -f faindex.linux.cl.make $(APP_COMP_OPTIONS)
	make -C kill -f kill.linux.cl.make $(APP_COMP_OPTIONS) PSILIB=$(PSILIB)
	make -C makeParamDB -f make_param_db.linux.cl.make $(APP_COMP_OPTIONS)
	make -C msdisplay -f msdisplay.linux.cl.make $(APP_COMP_OPTIONS)
	make -C msform -f msform.linux.cl.make $(APP_COMP_OPTIONS)
	make -C mssearch -f mssearch.linux.cl.make $(APP_COMP_OPTIONS) LLIB=$(LLIB)
	make -C searchCompare -f search_compare.linux.cl.make $(APP_COMP_OPTIONS) LLIB=$(LLIB)

allclmpi:
	make -C libzip -f libzip.linux.make $(C_LIB_COMP_OPTIONS)
	make -C libsqlite -f libsqlite.linux.make $(C_LIB_COMP_OPTIONS)
	make -C libexpat -f libexpat.linux.make $(LIB_COMP_OPTIONS)
	make -C libgen -f libgen.linux.make $(LIB_COMP_OPTIONS)
	make -C libnrec -f libnrec.linux.make $(LIB_COMP_OPTIONS)
	make -C libprocess -f libsingle.linux.make $(LIB_COMP_OPTIONS)
	rm -rf libprocess/*.o
	make -C libprocess -f libmulti.linux.make $(MPI_LIB_COMP_OPTIONS)
	make -C libraw -f libraw.linux.make $(LIB_COMP_OPTIONS)
	make -C libucsf -f libucsf.linux.cl.make $(LIB_COMP_OPTIONS)

	make -C faindex -f faindex.linux.cl.make $(APP_COMP_OPTIONS)
	make -C kill -f kill.linux.cl.make $(APP_COMP_OPTIONS) PSILIB=$(PSILIB)
	make -C makeParamDB -f make_param_db.linux.cl.make $(APP_COMP_OPTIONS)
	make -C msdisplay -f msdisplay.linux.cl.make $(APP_COMP_OPTIONS)
	make -C msform -f msform.linux.cl.make $(APP_COMP_OPTIONS)
	make -C mssearch -f mssearch.linux.cl.make $(APP_COMP_OPTIONS) LLIB=$(LLIB)
	rm -rf mssearch/*.o
	make -C mssearch -f mssearchmpi.linux.cl.make $(MPI_APP_COMP_OPTIONS) MPILIB=$(MPILIB) PSILIB=$(PSILIB) LLIB=$(LLIB)
	make -C searchCompare -f search_compare.linux.cl.make $(APP_COMP_OPTIONS) LLIB=$(LLIB)

clean:
	rm -rf libzip/*.o
	rm -rf libsqlite/*.o
	rm -rf libexpat/*.o
	rm -rf libdbase/*.o
	rm -rf libgen/*.o
	rm -rf libnrec/*.o
	rm -rf libprocess/*.o
	rm -rf libraw/*.o
	rm -rf libucsf/*.o

	rm -rf btag-daemon/*.o
	rm -rf faindex/*.o
	rm -rf jobStatus/*.o
	rm -rf kill/*.o
	rm -rf makeParamDB/*.o
	rm -rf msbatch/*.o
	rm -rf msdisplay/*.o
	rm -rf msform/*.o
	rm -rf mssearch/*.o
	rm -rf readDB/*.o
	rm -rf searchCompare/*.o
	rm -rf bin/*
	rm -f lib/libzip.a
	rm -f lib/libsqlite.a
	rm -f lib/libexpat.a
	rm -f lib/libdbase.a
	rm -f lib/libgen.a
	rm -f lib/libnrec.a
	rm -f lib/libsingle.a
	rm -f lib/libmulti.a
	rm -f lib/libraw.a
	rm -f lib/libucsf.a
