# ~~~~ Neha
# modified on 15 July 2014, Neha Khetan
# system specifications:  Linux cathale 3.9.10-100.fc17.x86_64 #1 SMP Sun Jul 14 01:31:27 UTC 2013 x86_64 x86_64 x86_64 GNU/Linux

# ............... Processor	: IntelXenon(R): i(ntel® Xeon(R) CPU E5645 @ 2.40GHz × 12 )
# ............... System	: Linux
# ............... machine	: 64 bit: x86_64
# ............... platform	: Linux-3.6.11-1.fc17.x86_64-x86_64-with-fedora-17-Beefy_Miracle
# ............... compiler	: g++ (GCC) 4.7.2 20120921 (Red Hat 4.7.2-2)

### defining the MACHINE as linXenon: for linux and Xenon


### --------   Modified lines: Line # from (197-204) has been commented out from that of the original
###	       and replaced with line (# 208 and 210)
# ~~~~~~~~~~~~~~~~~~~~~~  End of changes made by neha

### CA
#RCS: $Id: makefile.mk,v 2.38 2005/04/25 08:05:09 nedelec Exp $
#include file for all makefiles of cytosim,
#Francois Nedelec 2000, Dietrich Foethke 2004.Chaitanya Athale 2008

#-----------------Define the machine and the compiler here:--------------------
#note: not all combinations of machines and compilers do exist in the makefile!
#available machines are:  linopt linppc linath linp3 linp4 linpm bibo macG4 macG5 win macintel4
#available compilers are: g++ icpc xlc++




# machine type
MACHINE=linXenon

# compiler
CPP_COMPILER=g++

# fortran compiler
FORTRAN_COMPILER=g77


#-----------------define the compile mode: G)eneric F)ast, P)rofiling, D)ebug

MODE=F

#-----------------common compile options for icc, gcc, fast and profiling

ICCOPT=-wd1572
GCCOPT=
FASTOPT=-fno-rtti -ffast-math -fno-trapping-math

#-----------------paths to the math-kernel libraries---------------------------

MKL=/opt/intel/mkl/lib/32

#-----------------flags for the fortran-compiler (valid for all machines):

FORTRAN_FLAGS=-O3
G2C=/usr/lib/gcc-lib/i386-redhat-linux/3.2.2

#-----------------Static linking for sim on Linux:

LL_static= $(MKL)/libmkl_lapack.a $(MKL)/libmkl_ia32.a $(MKL)/libguide.a /usr/lib/libpthread.a /usr/lib/libpthread_nonshared.a $(G2C)/libg2c.a -static

LL_static= lapack.a /usr/lib/libpthread.a /usr/lib/libpthread_nonshared.a $(G2C)/libg2c.a -static


#------------------------------------------------------------------------------
#--------------------------MACHINE-SPECIFIC OPTIONS----------------------------


#-----------------extension for copying the binaries to the bin directory
#                 (usually empty, only needed for windows so far)

BINEXT=

#-----------------Linux Athlon, g++ compiler:

ifeq ($(MACHINE),linath)

  CDg++FLAGS=$(GCCOPT) -O0 -g3 -ggdb -Wall
  CPg++FLAGS=$(GCCOPT) -O0 -fno-inline -pg -march=athlon -malign-double $(FASTOPT) 
  CCg++FLAGS=$(GCCOPT) -O0 -fprofile-arcs -ftest-coverage -march=athlon $(FASTOPT)
  CFg++FLAGS=$(GCCOPT) -O3 -finline-functions -march=athlon -malign-double $(FASTOPT)

  #---MKL lapack & blas
  #LL= -L$(MKL) -lmkl_lapack -lmkl_p4 -lguide -lg2c
  #---for standard blas and lapack:
  LL= -llapack -lblas -lg2c
  LGL= -L/usr/lib -L/usr/X11R6/lib -lglut -lGL -lGLU -lXt -lX11 -lXi -lXmu

endif

#-----------------Linux Pentium3 with intel/g++ compiler:

ifeq ($(MACHINE),linp3)

  CDicpcFLAGS=-g -O0 -inline_debug_info
  CPicpcFLAGS=-O2 -unroll -pg
  CFicpcFLAGS=-O3 -tpp6 -unroll

  CDg++FLAGS=-O0 -g3 -ggdb -Wall
  CPg++FLAGS=-O0 -fno-inline -pg -march=i686 -malign-double $(FASTOPT)
  CCg++FLAGS=-O0 -fprofile-arcs -ftest-coverage -march=i686  $(FASTOPT)
  CFg++FLAGS=-O3 -finline-functions -march=i686 -malign-double $(FASTOPT)

  #---MKL lapack & blas
  #LL= -L$(MKL) -lmkl_lapack -lmkl_p3 -lguide -lg2c
  #---for standard blas and lapack:
  LL= -llapack -lblas -lg2c
  LGL=-L/usr/lib -L/usr/X11R6/lib -lglut -lGL -lGLU -lXt -lX11 -lXi -lXmu

endif

#-----------------Linux Pentium4, intel/g++ compiler:

ifeq ($(MACHINE),linp4)

  CDicpcFLAGS=$(ICCOPT) -g
  CPicpcFLAGS=$(ICCOPT) -O2 -unroll -pg
  CFicpcFLAGS=$(ICCOPT) -O3 -tpp7 -unroll -fno-rtti

  CDg++FLAGS=$(GCCOPT) -O0 -g3 -ggdb -Wall
  CPg++FLAGS=$(GCCOPT) -O0 -fno-inline -pg -march=pentium4 -malign-double $(FASTOPT) 
  CCg++FLAGS=$(GCCOPT) -O0 -fprofile-arcs -ftest-coverage -march=pentium4 $(FASTOPT)
  CFg++FLAGS=$(GCCOPT) -O3 -finline-functions -march=pentium4 -malign-double $(FASTOPT)

  #---MKL version 6 and 7 lapack & blas
  #LL= -L$(MKL) -lmkl_lapack -lmkl_p4 -lguide
  #---MKL version 5 lapack & blas
  #LL= -L$(MKL) -lmkl_lapack -lmkl_p4 -lguide -lg2c
  #---for standard blas and lapack:
  LL= -llapack -lblas -lg2c
  LGL= -L/usr/lib -L/usr/X11R6/lib -lglut -lGL -lGLU -lXt -lX11 -lXi -lXmu

endif

#-----------------Linux Pentium M, intel/g++ compiler:

ifeq ($(MACHINE),linpm)

  CDicpcFLAGS=$(ICCOPT) -g
  CPicpcFLAGS=$(ICCOPT) -O2 -unroll -pg
  CFicpcFLAGS=$(ICCOPT) -O3 -tpp7 -unroll -fno-rtti

  CDg++FLAGS=$(GCCOPT) -O0 -g3 -ggdb -Wall
#  CPg++FLAGS=$(GCCOPT) -O0 -fno-inline -pg -march=pentium-m -malign-double $(FASTOPT) 
#  CCg++FLAGS=$(GCCOPT) -O0 -fprofile-arcs -ftest-coverage -march=pentium-m $(FASTOPT)
#  CFg++FLAGS=$(GCCOPT) -O3 -finline-functions -march=pentium-m -malign-double $(FASTOPT)
  CPg++FLAGS=$(GCCOPT) -O0 -fno-inline -pg -march=pentium3 -msse2 -malign-double $(FASTOPT) 
  CCg++FLAGS=$(GCCOPT) -O0 -fprofile-arcs -ftest-coverage -march=pentium3 -msse2 $(FASTOPT)
  CFg++FLAGS=$(GCCOPT) -O3 -finline-functions -march=pentium3 -msse2 -malign-double $(FASTOPT)

  #---MKL version 6 and 7 lapack & blas
  #LL= -L$(MKL) -lmkl_lapack -lmkl_p3 -lguide
  #---MKL version 5 lapack & blas
  #LL= -L$(MKL) -lmkl_lapack -lmkl_p3 -lguide -lg2c
  #---for standard blas and lapack:
  LL= -llapack -lblas -lg2c
  LGL= -L/usr/lib -L/usr/X11R6/lib -lglut -lGL -lGLU -lXt -lX11 -lXi -lXmu

endif#-----------------Linux PPC, g++/IBM compiler:

ifeq ($(MACHINE),linppc)

  CDg++FLAGS=-O0 -g3 -ggdb -Wall
  CPg++FLAGS=-O0 -fno-inline -pg -mcpu=970 $(FASTOPT) 
  CCg++FLAGS=-O0 -fprofile-arcs -ftest-coverage $(FASTOPT)
  CFg++FLAGS=-O3 -finline-functions -mcpu=970 $(FASTOPT)

  CDxlc++FLAGS=-O0 -g -qflag=w
  CFxlc++FLAGS=-O3 -qtune=ppc970

  #BLAS_PPC=/home/foethke/cb1/src/LAPACK
  ATLAS_PPC=/usr/local/ATLAS/lib/Linux_PPCG4AltiVec_2


  #---for standard blas and lapack: this lib are only on ibm-node57
  LL= -llapack -lblas -lg2c
  #---for ATLAS blas
  #LL= -llapack $(ATLAS_PPC)/libf77blas.a $(ATLAS_PPC)/libatlas.a -lg2c

  LGL= -L/usr/lib -L/usr/X11R6/lib -lglut -lGL -lGLU -lXt -lX11 -lXi -lXmu

endif

#-----------------Linux AMD Opteron, g++ compiler:

ifeq ($(MACHINE),linXenon)

  CDg++FLAGS=-O0 -g3 -ggdb -Wall
  CPg++FLAGS=-O0 -fno-inline -pg -march=opteron -m64 $(FASTOPT)
  CCg++FLAGS=-O0 -fprofile-arcs -ftest-coverage $(FASTOPT)
  CFg++FLAGS=-O3 -finline-functions -march=opteron -m64 $(FASTOPT)

 
  #### changes made my nk on 15 July 2014: 
  ## the following lines have been commented out (197-204) and replaced with line (# 208 and 210).
 
   #--- nk: BLAS_OPT=/home/foethke/cb1/src/LAPACK
  #ATLAS_OPT=/usr/local/ATLAS/lib/Linux_PPCG4AltiVec_2
  #---for standard blas and lapack:
  # LL= $(BLAS_OPT)/lapack_OPT.a $(BLAS_OPT)/blas_OPT.a -lg2c
  #---for ATLAS blas
  #LL= -llapack $(ATLAS_PPC)/libf77blas.a $(ATLAS_PPC)/libatlas.a -lg2c
  #LINK= $(STATIC) -L/g/nedelec/opt/netlib/linux_em64t -llapack -lblas -lgfortran 
  #LGL= -L/usr/lib -L/usr/X11R6/lib -lglut -lGL -lGLU -lXt -lX11 -lXi -lXmu
  
  
  
  # for LAPACK
  LL= -llapack -lblas -lgfortran
  LGL= -L/usr/lib -L/usr/lib64/X11 -lglut -lGL -lGLU -lXt -lX11 -lXi -lXmu
  
  

endif

#-----------------Linux Pentium4, intel compiler:

ifeq ($(MACHINE),bibo)

  CDicpcFLAGS=-
  CPicpcFLAGS=-O2 -unroll -p
  CFicpcFLAGS=-O3 -tpp7 -unroll $(FASTOPT)

  LL = -L$(MKL) -lmkl_lapack -lmkl_p4 -lguide /usr/local/gcc-3.2.3/lib/libg2c.so
  LGL= -L/usr/lib -L/usr/X11R6/lib -lglut -lGL -lGLU -lXt -lX11 -lXi -lXmu

endif

#-----------------MacG4 ( -mcpu=7450 )

ifeq ($(MACHINE),macG4)

  CDg++FLAGS=-g3 -Wall -mdynamic-no-pic
  CPg++FLAGS=-g3 -Wall -arch ppc -mcpu=7450 -O3 -falign-loops=16 $(FASTOPT) -mdynamic-no-pic
  CFg++FLAGS=-arch ppc -mcpu=7450 -falign-loops=16 -fast $(FASTOPT)
  LL  = -framework vecLib
  LGL = -framework GLUT -framework openGL -framework Foundation -framework AGL

endif

#-----------------MacG5 ( -mcpu=970 )

ifeq ($(MACHINE),macG5)

  CDg++FLAGS=-g3 -Wall -mdynamic-no-pic $(FASTOPT)
  CPg++FLAGS=-Wall -arch ppc -O3 -fast $(FASTOPT)
  #CPg++FLAGS=-g3 -Wall -arch ppc -fast
  #CFg++FLAGS=-g3 -ggdb3 -O1 -arch ppc -fast $(FASTOPT)
  CFg++FLAGS=-O3 -arch ppc -fast $(FASTOPT)
  LL  = -framework vecLib
  LGL = -framework GLUT -framework openGL -framework Foundation -framework AGL

endif

#-----------------windows PC ( -mcpu=evil_bill )

ifeq ($(MACHINE),win)

  #--- use a home-made blas/lapack for distribution:
  # put the libraries libblas_win32.a and liblapack_win32 in /cytosim/lib


  CDg++FLAGS=-g 
  CPg++FLAGS=-pg -O3 -march=i686 -ffast-math -malign-double
  CFg++FLAGS=-O3 -finline-functions -march=i686 -ffast-math -malign-double
  LL= -llapack_win32 -lblas_win32 -lgen -lg2c

  LGL=-lglut32 -lglu32 -lopengl32
  BINEXT=.exe

endif



#-----------------MAC Intel ( -mcpu=macintel )

ifeq ($(MACHINE),macintel)

  CDicpcFLAGS=$(ICCOPT) -g
  CPicpcFLAGS=$(ICCOPT) -O2 -unroll -pg
  CFicpcFLAGS=$(ICCOPT) -O3 -tpp7 -unroll -fnolib-rtti

  CDg++FLAGS=$(GCCOPT) -O0 -g3 -ggdb -Wall
  CPg++FLAGS=$(GCCOPT) -O0 -fno-inline -pg -march=pentium4 -malign-double $(FASTOPT) 
  CCg++FLAGS=$(GCCOPT) -O0 -fprofile-arcs -ftest-coverage -march=pentium4 $(FASTOPT)
  CFg++FLAGS=$(GCCOPT) -O3 -finline-functions -march=pentium4 -malign-double $(FASTOPT)

  #---for standard blas and lapack:
  LL= -llapack -lblas -lgcc  
  
  #--- openGL
  LGL= -L/usr/lib -L/usr/X11R6/lib -lglut -lGL -lGLU -lXt -lX11 -lXi -lXmu -framework GLUT -framework openGL -framework AGL


endif

#-----------------MAC Intel quad core( -mcpu=macintel4 )

ifeq ($(MACHINE),macintel4)

  CDicpcFLAGS=$(ICCOPT) -g
  CPicpcFLAGS=$(ICCOPT) -O2 -unroll -pg
  CFicpcFLAGS=$(ICCOPT) -O3 -tpp7 -unroll -fnolib-rtti

  CDg++FLAGS=$(GCCOPT) -O0 -g3 -ggdb -Wall
  CPg++FLAGS=$(GCCOPT) -O0 -fno-inline -pg -march=core2 -malign-double $(FASTOPT) 
  CCg++FLAGS=$(GCCOPT) -O0 -fprofile-arcs -ftest-coverage -march=pentium4 $(FASTOPT)
  CFg++FLAGS=$(GCCOPT) -O3 -finline-functions -march=core2 $(FASTOPT)
  #-malign-double $(FASTOPT)

  #---for standard blas and lapack:
  LL= -llapack -lblas -lgcc  
  
  #--- openGL
  LGL= -L/usr/lib -L/usr/X11R6/lib -lglut -lGL -lGLU -lXt -lX11 -lXi -lXmu -framework GLUT -framework openGL -framework AGL


endif


#-------------------------------------------
#from 2009 version of cytosim

ifeq ($(MACHINE),mac)

    FlagsD = $(DEBGOPT) $(WARNOPT)
    FlagsP = $(DEBGOPT) -g -O3 $(FASTOPT) $(ALIGNOPT)
    FlagsF = -O3 -mdynamic-no-pic $(FASTOPT) $(ALIGNOPT)
  
    LINK   = -framework vecLib  -framework Foundation
    LINKGL = -framework GLUT -framework openGL -framework AGL
    
endif


#-----------------end of machine-specific options------------------------------
#------------------------------------------------------------------------------




#-----------------lists of objects to be made----------------------------------
SPACEOBJS=space.o space_square.o space_sphere.o space_oval.o \
          space_roundsquare.o space_periodic.o space_strip.o \
          space_ellipse.o space_cylinder.o space_cylinderZ.o space_tee.o \
          space_banana.o space_boomerang.o \
          space_inflate.o space_deflate.o space_combine.o \

BASEOBJS=iomessages.o iowrapper.o exceptions.o random.o smath.o \
	 vecteur1.o vecteur2.o vecteur3.o matrix1.o matrix2.o matrix3.o \
	 node.o node_list.o name_list.o rasterizer.o map.o \
	 matrix.o matsparse.o matsparsesym.o  matsparsebandsym.o \
	 matsym.o matblock.o matsparsesym2.o \
	 quaternion.o gmres.o conjgradient.o pointsonsphere.o \
	 parameter.o parameter_typed.o parameter_list.o poly_roots.o\
	 $(SPACEOBJS)


SIMOBJS=sim_param.o object.o object_list.o pointset.o simobject.o\
	point_exact.o point_interpolated.o \
        solid.o fiber.o microtub.o point_microtub.o\
	hand.o grafted.o complex.o microtub_organizer.o aster.o nucleus.o\
        aster_list.o complex_list.o grafted_list.o microtub_list.o \
        nucleus_list.o solid_list.o sim.o lattice.o

PLAYOBJS=reader.o play_param.o glextensions.o globjects3d.o offscreen.o play.o

TESTOBJS=testglut.o testglut2.o testglut2d.o testmap.o testparam.o testspace.o\
         testquaternion.o testrandom.o testrasterizer.o testreader.o  \
         testsolve.o testsolved.o testnucleus.o testsphere.o testpolyroots.o

TOOLOBJS=analyse0.o analyse1.o analyse2.o analyse3.o analyse4.o \
         datagen.o datagen1.o datagen2.o datagen3.o datagen4.o \
         datagen5.o datagen6.o datagen7.o datagen8.o datagen9.o \
         frametool.o meanlength.o readwrite.o readwritedata.o \
         reducedata.o createstart.o analyse_nucpos.o

#-----------------lists of all the source-files for base, sim and play---------

BASEHEADERS=$(wildcard $(BASEDIR)/*.h)
BASESOURCES=$(wildcard $(BASEDIR)/*.cc)

SIMHEADERS=$(wildcard $(SIMDIR)/*.h)
SIMSOURCES=$(wildcard $(SIMDIR)/*.cc)

PLAYHEADERS=$(wildcard $(PLAYDIR)/*.h)
PLAYSOURCES=$(wildcard $(PLAYDIR)/*.cc)


ALLSIMHEADERS=$(SIMHEADERS) $(BASEHEADERS)
ALLSIMSOURCES=$(SIMSOURCES) $(BASESOURCES)

ALLPLAYHEADERS=$(ALLSIMHEADERS) $(PLAYHEADERS)
ALLPLAYSOURCES=$(ALLSIMSOURCES) $(PLAYSOURCES)

#-----------------lists of all the programs to be made-------------------------

#the targets in sim
ALLSIMPROGS=sim sim1 sim2

#the targets in play
ALLPLAYPROGS=play

#the test-programs
#SIMPLE_TESTS=
#SIMPLE_GLTESTS=testglut2 testglut2d
#ALLTESTS=$(SIMPLE_TESTS) $(SIMPLE_GLTESTS)
ALLTESTS=\
testglut testglut2 testglut2d testreader testsolve testsolved testnucleus \
testparam testquaternion testrasterizer testspace testmap testrandom \
testsphere testpolyroots

#the tools
ALLTOOLS=frametool readwrite readwritedata meanlength reducedata \
datagen datagen1 datagen2 datagen3 datagen4 datagen5 datagen6 datagen7 \
datagen8 datagen9 analyse0 analyse1 analyse2 analyse3 analyse4 createstart \
analyse_nucpos

#everything
EVERYTHING=$(ALLSIMPROGS) $(ALLPLAYPROGS) $(ALLTESTS) $(ALLTOOLS)


#==============================================================================
#              no user-serviceable parts beyond this point!
#==============================================================================


#-----------------prepare the compiler-flags and the linkage options $(L)------

ifeq ($(LL),)
  $(error No linkage-options defined for $$(MACHINE)=$(MACHINE)!)
endif

L=$(LL)


#-----------------set the compilers and the compiler-flags---------------------

ifeq ($(C$(MODE)$(CPP_COMPILER)FLAGS),)
$(error No compiler-options defined for the combination $$(MACHINE)=$(MACHINE)\
 and $$(CPP_COMPILER)=$(CPP_COMPILER)!)
endif

CXX=$(CPP_COMPILER)
CXXFLAGS=$(C$(MODE)$(CPP_COMPILER)FLAGS) $(INCLUDE)

FC=$(FORTRAN_COMPILER)
FFLAGS=$(FORTRAN_FLAGS)


#-----------------Rules--------------------------------------------------------

%.o: %.cc
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $(BUILDDIR)/$@ $<

%.o: %.f
	$(FC) -c $(FFLAGS) -o $(BUILDDIR)/$@ $<

