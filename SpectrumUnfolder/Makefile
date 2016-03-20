#
# General Makefile for the OfflineUser package
#
#

# Replace the wildcard expression with .cc file list if you do
# not want to compile all .cc files in this directory
#
TOPDIR= $(shell pwd)

BINDIR  = $(TOPDIR)/bin
LIBDIR  = $(TOPDIR)/lib
SRCDIR  = $(TOPDIR)/src
INCDIR  = $(TOPDIR)/include
OBJDIR  = $(TOPDIR)/obj

USER_SRCS = $(wildcard $(SRCDIR)/*.cc)
files := $(foreach USER_SRCS,$(SRCDIR),$(USER_SRCS))

HEADERS_DICT = $(INCDIR)/MathUtils.h $(INCDIR)/SpectrumUtils.h $(INCDIR)/ForwardFolder.h

OBJS = $(USER_SRCS:.cc=.o)

MAIN_UNFOLD = UnfoldSpectrum.cc
MAIN_UNFOLD_OBJ = $(MAIN_UNFOLD:.cc=.o)

## Get platform 32/64 bit
LBITS   = $(shell getconf LONG_BIT)

## Get SVN revision number
SVN_REVISION_VERSION = -D'SVN_REVISION_VERSION="$(shell svnversion -n .)"'

# Set executable a name
EXE = libSpectrumUnfolder.so
UNFOLD_EXE = UnfoldSpectrum
#
#############################################################

## You should not need to change anything below this line ###

.PHONY: all depend clean


######################################
###  CPPFLAGS & CXXFLAGS  ############
######################################

CPPFLAGS = -I$(INCDIR)

ifeq ($(LBITS),64)
  # do 64 bit stuff here
	CPPFLAGS += -I/usr/include -pthread -m64
	CXXFLAGS = -std=c++11 -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m64 
	##CXXFLAGS = -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m64 
	SYSLIBDIR = /usr/lib/x86_64-linux-gnu
else
  # do 32 bit stuff here
	CPPFLAGS += -I/usr/include -pthread -m32
  CXXFLAGS = -std=c++11 -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m32 
	##CXXFLAGS = -O2 -Wall -Wunused -Wuninitialized -fPIC -pthread -m32 
	SYSLIBDIR = /usr/lib
endif

SOFLAGS = -fPIC -ggdb3 -Wall -shared

ROOT_CPPFLAGS = -I$(ROOTSYS)/include
CPPFLAGS += $(ROOT_CPPFLAGS) 

###########################
###  LDFLAGS   ############
###########################
ROOT_LDFLAGS = $(shell root-config --libs) -L$(ROOTSYS)/lib -lSpectrum -lMathMore -lMinuit
SYSTEM_LDFLAGS = -L$(SYSLIBDIR) -lrt 

LDFLAGS = $(ROOT_LDFLAGS) 

################################################################


all: GETOBJS PRINTINFO $(EXE) PUTOBJS $(UNFOLD_EXE)

PRINTINFO: 
	@echo 'Compiling $(SVN_REVISION_VERSION) trunk version on a $(LBITS) bit machine' \

REMOVEOBJS:
	- rm -f *.o $(SRCDIR)/*.o 

GETOBJS:
	@echo "Put object and lib files again to $(SRCDIR) dir"
	- mv -f $(OBJDIR)/*.o $(SRCDIR) 2>/dev/null
	- mv -f $(LIBDIR)/*.so $(TOPDIR) 2>/dev/null
	- mv -f $(LIBDIR)/ClassDictionary_rdict.pcm $(TOPDIR) 2>/dev/null


PUTOBJS:
	@echo "Moving object files to $(OBJDIR) dir and libraries to $(LIB) dir"
	- mv -f $(SRCDIR)/*.o $(OBJDIR) 2>/dev/null
	- mv -f *.o $(OBJDIR) 2>/dev/null
	- mv -f $(EXE) $(LIBDIR) 2>/dev/null
	- mv -f ClassDictionary_rdict.pcm $(LIBDIR) 2>/dev/null

ClassDictionary.cc: $(HEADERS_DICT) LinkDef.h
	rootcint -f $@ -c -p $(CPPFLAGS) $^
	
$(EXE): $(OBJS) ClassDictionary.o
	@echo "Compiler is $(CXX) or $(CC), options are $(CXXFLAGS), generating $@ shared library..."
	@$(CXX) $(CXXFLAGS) $(SOFLAGS) $(CPPFLAGS) $^ -o $@ $(LDFLAGS) 
	

$(UNFOLD_EXE): $(MAIN_UNFOLD_OBJ)
	@echo "Building unfold spectrum ..."
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(MAIN) $(LDFLAGS) -L$(LIBDIR) -lSpectrumUnfolder

#############################################################
# gcc can generate the dependency list

depend: Make-depend

Make-depend: $(USER_SRCS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $^ > $@

clean:
	- rm -f $(EXE) $(UNFOLD_EXE) *.o $(OBJDIR)/*.o $(SRCDIR)/*.o ClassDictionary.cc ClassDictionary.h $(LIBDIR)/ClassDictionary_rdict.pcm $(BINDIR)/*.o *.so $(LIBDIR)/*.so *.ps core Make-depend


#############################################################
# 'make run' will run the thing



#############################################################
# the lines below are for running with debugger 'make run_gdb'

.INTERMEDIATE: gdb.cmdl

# batch mode gdb needs a file with commands
gdb.cmdl:
	echo "r -b bootstrap.xml" > $@

run_gdb: gdb.cmdl $(EXE)
	gdb -batch -x $< ./$(EXE) && touch $@

#include Make-depend
