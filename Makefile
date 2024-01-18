#Some Makefile for CLASS.
#Julien Lesgourgues, 28.11.2011
#Nils Schöneberg, Matteo Lucca, 27.02.2019
#Boris Bolliet since 2015 - for the class_sz parts.


MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
OUTDIR = $(MDIR)/output

.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base
	if ! [ -e $(OUTDIR) ]; then mkdir $(OUTDIR) ; fi;

vpath %.c source:tools:main:test
vpath %.o build
vpath .base build

########################################################
###### LINES TO ADAPT TO YOUR PLATEFORM ################
########################################################

# your C compiler:
# CC       = gcc-12
# on Mac M1:
CC       = /usr/bin/clang

#CC       = gcc  -Wunused-variable
#CC       = icc
#CC       = pgcc

# your tool for creating static libraries:
AR        = ar rv

# Your python interpreter.
# In order to use Python 3, you can manually
# substitute python3 to python in the line below, or you can simply
# add a compilation option on the terminal command line:
# "PYTHON=python3 make all" (Thanks to Marius Millea for python3 compatibility)
PYTHON ?= python
#PYTHON ?= /Users/boris/opt/anaconda2/bin/python

# your optimization flag
# on Mac M1
OPTFLAG = -O4 -ffast-math #-arch x86_64
#OPTFLAG = -O3 # on v2.10.3
#OPTFLAG = -Ofast -ffast-math #-march=native
#OPTFLAG = -fast

# your openmp flag (comment for compiling without openmp)
#OMPFLAG   = -fopenmp
# on Mac M1
OMPFLAG   = -Xclang -fopenmp
#OMPFLAG   = -mp -mp=nonuma -mp=allcores -g
#OMPFLAG   = -openmp

# all other compilation flags
CCFLAG = -g -fPIC
LDFLAG = -g -fPIC

#on Mac M1
LDFLAG += -lomp

# leave blank to compile without HyRec, or put path to HyRec directory
# (with no slash at the end: e.g. "external/RecfastCLASS")
HYREC = external/HyRec2020
RECFAST = external/RecfastCLASS
HEATING = external/heating

#
# GSL          = gsl
# GSL_INC_PATH = /usr/local/include/
# GSL_LIB_PATH = /usr/local/lib/
# LIBSGSL = -L$(GSL_LIB_PATH) -l$(GSL)

########################################################
###### IN PRINCIPLE THE REST SHOULD BE LEFT UNCHANGED ##
########################################################

# pass current working directory to the code
CCFLAG += -D__CLASSDIR__='"$(MDIR)"'

# where to find include files *.h
INCLUDES = -I../include -I/usr/local/include/ -I/Users/boris/opt/miniconda3/include/
HEADERFILES = $(wildcard ./include/*.h)

# automatically add external programs if needed. First, initialize to blank.
EXTERNAL =

vpath %.c $(RECFAST)
#CCFLAG += -DRECFAST
INCLUDES += -I../$(RECFAST)
EXTERNAL += wrap_recfast.o
HEADERFILES += $(wildcard ./$(RECFAST)/*.h)

vpath %.c $(HEATING)
#CCFLAG += -DHEATING
INCLUDES += -I../$(HEATING)
EXTERNAL += injection.o noninjection.o
HEADERFILES += $(wildcard ./$(HEATING)/*.h)

# update flags for including HyRec
ifneq ($(HYREC),)
vpath %.c $(HYREC)
CCFLAG += -DHYREC
#LDFLAGS += -DHYREC
INCLUDES += -I../$(HYREC)
EXTERNAL += hyrectools.o helium.o hydrogen.o history.o wrap_hyrec.o energy_injection.o
HEADERFILES += $(wildcard ./$(HYREC)/*.h)
endif

%.o:  %.c .base $(HEADERFILES)
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o  evolver_ndf15.o arrays.o parser.o quadrature.o hyperspherical.o common.o trigonometric_integrals.o r8lib.o class_sz_tools.o class_sz_custom_profiles.o Patterson.o fft.o

SOURCE = input.o background.o thermodynamics.o perturbations.o primordial.o nonlinear.o transfer.o spectra.o lensing.o distortions.o class_sz.o class_sz_clustercounts.o

INPUT = input.o

PRECISION = precision.o

BACKGROUND = background.o

THERMO = thermodynamics.o

PERTURBATIONS = perturbations.o

TRANSFER = transfer.o

PRIMORDIAL = primordial.o

SPECTRA = spectra.o

NONLINEAR = nonlinear.o

LENSING = lensing.o

DISTORTIONS = distortions.o

OUTPUT = output.o

CLASS_SZ = class_sz_driver.o

TEST_LOOPS = test_loops.o

TEST_LOOPS_OMP = test_loops_omp.o

TEST_SPECTRA = test_spectra.o

TEST_TRANSFER = test_transfer.o

TEST_NONLINEAR = test_nonlinear.o

TEST_PERTURBATIONS = test_perturbations.o

TEST_THERMODYNAMICS = test_thermodynamics.o

TEST_BACKGROUND = test_background.o

TEST_HYPERSPHERICAL = test_hyperspherical.o

C_TOOLS =  $(addprefix tools/, $(addsuffix .c,$(basename $(TOOLS))))
C_SOURCE = $(addprefix source/, $(addsuffix .c,$(basename $(SOURCE) $(OUTPUT))))
C_TEST = $(addprefix test/, $(addsuffix .c,$(basename $(TEST_DEGENERACY) $(TEST_LOOPS) $(TEST_TRANSFER) $(TEST_NONLINEAR) $(TEST_PERTURBATIONS) $(TEST_THERMODYNAMICS))))
C_MAIN = $(addprefix main/, $(addsuffix .c,$(basename $(CLASS_SZ))))
C_ALL = $(C_MAIN) $(C_TOOLS) $(C_SOURCE)
H_ALL = $(addprefix include/, common.h svnversion.h $(addsuffix .h, $(basename $(notdir $(C_ALL)))))
PRE_ALL = cl_ref.pre clt_permille.pre
INI_ALL = explanatory.ini lcdm.ini
MISC_FILES = Makefile CPU psd_FD_single.dat myselection.dat myevolution.dat README bbn/sBBN.dat external_Pk/* cpp
PYTHON_FILES = python/classy.pyx python/setup.py python/cclassy.pxd python/test_class.py

all: class_sz libclass.a classy_sz

libclass.a: $(TOOLS) $(SOURCE) $(EXTERNAL)
	$(AR)  $@ $(addprefix build/, $(TOOLS) $(SOURCE) $(EXTERNAL))

class_sz: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(CLASS_SZ)
	#$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -g -o class $(addprefix build/,$(notdir $^)) -lm -L/home/runner/work/SOLikeT/SOLikeT/gsl-2.6/lib -lgsl -lgslcblas
	 # $(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -g -o class $(addprefix build/,$(notdir $^)) -L/usr/local/lib -L/Users/boris/miniconda/lib -lgsl -lgslcblas -lfftw3 -lm
	 $(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -g -o class_sz $(addprefix build/,$(notdir $^)) -lgsl -lgslcblas -lfftw3 -lm -L/Users/boris/opt/miniconda3/lib

	 # $(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -g -o class $(addprefix build/,$(notdir $^)) -L/usr/local/lib -L/Users/boris/opt/anaconda3/lib -L/opt/homebrew/lib -lgsl -lgslcblas -lfftw3 -lm
	 # $(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -g -o class $(addprefix build/,$(notdir $^)) -L/usr/local/lib -lgsl -lgslcblas -lm

test_loops: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_LOOPS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_loops_omp: $(TOOLS) $(SOURCE) $(EXTERNAL) $(OUTPUT) $(TEST_LOOPS_OMP)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix build/,$(notdir $^)) -lm

test_spectra: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_SPECTRA)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_transfer: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_TRANSFER)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_nonlinear: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_NONLINEAR)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_perturbations: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_PERTURBATIONS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_thermodynamics: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_THERMODYNAMICS)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_background: $(TOOLS) $(SOURCE) $(EXTERNAL) $(TEST_BACKGROUND)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o  $@ $(addprefix build/,$(notdir $^)) -lm

test_hyperspherical: $(TOOLS) $(TEST_HYPERSPHERICAL)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o test_hyperspherical $(addprefix build/,$(notdir $^)) -lm


tar: $(C_ALL) $(C_TEST) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)
	tar czvf class_sz.tar.gz $(C_ALL) $(H_ALL) $(PRE_ALL) $(INI_ALL) $(MISC_FILES) $(HYREC) $(PYTHON_FILES)

classy_sz: libclass.a python/classy.pyx python/cclassy.pxd
ifdef OMPFLAG
	cp python/setup.py python/autosetup.py
else
	grep -v "lgomp" python/setup.py > python/autosetup.py
endif
	cd python; export CC=$(CC); $(PYTHON) autosetup.py install || $(PYTHON) autosetup.py install --user
	rm python/autosetup.py

clean: .base
	rm -rf $(WRKDIR);
	rm -f libclass.a
	rm -f $(MDIR)/python/classy.c
	rm -rf $(MDIR)/python/build
	rm -f python/autosetup.py
