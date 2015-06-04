#=======================================================================
# User definitions
#

MAIN=timcom
export NETCDF=/home/share/cesm/software/esm-soft-150306

#======================================================================
# Libraries
#
ifndef NETCDF
nonetcdf:
	@echo ""
	@echo "!!!NETCDF variable does not set!!!"
	@echo ""
	@echo "  You need to specify the NETCDF root directory as below:"
	@echo "    for bash: export NETCDF=<netcdf root directory>"
	@echo "    for  csh: setenv NETCDF <netcdf root directory>"
	@echo ""
endif
INCS = -I$(NETCDF)/include
LIBS = -L$(NETCDF)/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -ldl -lm -lz -lcurl -lpnetcdf

include make.defs
F90FILES = $(AFILES) $(OFILES) $(RFILES) $(SFILES) $(TFILES) $(MFILES) $(GFILES) $(PFILES)



#=======================================================================
# Standard definitions
#

OFILESF90 = $(F90FILES:.f90=.o)
OFILESC = $(CFILES:.c=.o)


VPATH = obj
.SUFFIXES:
.SUFFIXES: .out .o .F90 .f90 .f .c .inc .h


#=======================================================================
# Flags 
#
FLAGS = $(INCS) 
CFLAGS = -O2

#=======================================================================
# Compiler Option
#
include ../../Compiler/make.compiler


#=======================================================================
# Targets and dependencies
#
default: opt
all: $(MAIN) 
debug:
	@make all "FLAGS = -C -g $(GENERAL_FLAGS) -D_DEBUG_"
#prof:
#	@make all "LIBS = -p $(LIBS)"
opt:
	@make all "FLAGS = $(GENERAL_FLAGS)"

$(MAIN): mod obj $(OFILESF90) $(OFILESC)
	@echo "$(LD) $(OFILESF90) $(OFILESC) $(LIBS)"
	@cd obj; $(LD) $(LDOPTION) $(OFILESF90) $(OFILESC) $(LIBS) -o $(@); mv $(@) ../;cd ../../;ln -fs src/$(@) ./

mod:
	mkdir mod
obj:
	mkdir obj
clean:
	rm -rf obj mod $(MAIN)
distclean:
	@make clean;rm -rf make.defs *.f90 *.c *.h

#=======================================================================
# Compilation rules
#
.F90.o:
	$(F90) $(FLAGS) $(FPPPFLAGS) $(LIBS)  -c $*.F90 -o obj/$*.o
.f90.o:
	$(F90) $(FLAGS) $(LIBS) -c $*.f90 -o obj/$*.o
.f.o:
	$(F77) $(FLAGS) $(LIBS) -c $*.f -o obj/$*.o
.c.o:
	$(CC) $(CFLAGS) $(LIBS) -c $*.c -o obj/$*.o

