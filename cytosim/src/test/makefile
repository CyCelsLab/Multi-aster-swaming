#RCS: $Id: makefile,v 2.5 2005/04/07 16:17:02 foethke Exp $
#makefile for the test-programs, Francois Nedelec 2000, Dietrich Foethke 2004.

ROOT=../..
SRCROOT=..
BASEDIR=$(SRCROOT)/base
SIMDIR=$(SRCROOT)/sim
PLAYDIR=$(SRCROOT)/play
TESTDIR=.
BUILDDIR=../../build
BINDIR=../../bin
LIBDIR=../../lib


INCLUDE=-I$(BASEDIR) -I$(SIMDIR) -I$(PLAYDIR)

include $(ROOT)/makefile.mk



#-----------------set the path for the prerequisites for .o files--------------

vpath %.o $(BUILDDIR)
vpath %.a $(LIBDIR)

#-----------------lists of all the files for the test-programs-----------------

TESTHEADERS=$(wildcard $(TESTDIR)/*.h)
TESTSOURCES=$(wildcard $(TESTDIR)/*.cc)

ALLTESTHEADERS=$(ALLPLAYHEADERS) $(TESTHEADERS)
ALLTESTSOURCES=$(ALLPLAYSOURCES) $(TESTSOURCES)

#----------patterns to grep the object and lib files from the prerequisites----

OBJPATTERN=$(patsubst %.o,$(BUILDDIR)/%.o, $(filter %.o, $(notdir $^)))
LIBPATTERN=$(patsubst lib%.a,-l%, $(filter %.a, $(notdir $^)))

#-----------------rules--------------------------------------------------------


#The variables $(SIMPLE_TESTS), $(SIMPLE_GLTESTS) and $(ALLTESTS)
#are set in $(ROOT)/makefile.mk.
#If you add a target in this makefile, please add it to the corresponding 
#variables, too!


.PHONY: all
all: $(ALLTESTS)



$(SIMPLE_TESTS): %: %.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
	@cp -f $@$(BINEXT) $(BINDIR)
$(SIMPLE_GLTESTS): %: %.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) $(LGL) -o $@
	@cp -f $@$(BINEXT) $(BINDIR)



testglut::
	@cd $(BASEDIR) && $(MAKE) random.o vecteur3.o quaternion.o pointsonsphere.o
testglut:: dep random.o vecteur3.o quaternion.o testglut.o pointsonsphere.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) $(LGL) -o $@
testglut::
	@cp -f $@$(BINEXT) $(BINDIR)

testglut2: dep testglut2.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) $(LGL) -o $@
	@cp -f $@$(BINEXT) $(BINDIR)

testglut2d: dep testglut2d.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) $(LGL) -o $@
	@cp -f $@$(BINEXT) $(BINDIR)

testreader::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
	@cd $(PLAYDIR) && $(MAKE) reader.o
testreader:: dep libcytobase.a libcytosim.a reader.o testreader.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -L$(LIBDIR) -lcytosim -lcytobase $(L) -o $@
testreader::
	@cp -f $@$(BINEXT) $(BINDIR)

testsolve::
	@cd $(BASEDIR) && $(MAKE) random.o vecteur2.o matrix2.o
testsolve:: dep random.o vecteur2.o testsolve.o matrix2.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) $(LGL) -o $@
testsolve::
	@cp -f $@$(BINEXT) $(BINDIR)

testsolved::
	@cd $(BASEDIR) && $(MAKE) random.o smath.o vecteur2.o
testsolved:: dep random.o smath.o vecteur2.o testsolved.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) $(L) $(LGL) -o $@
testsolved::
	@cp -f $@$(BINEXT) $(BINDIR)

testnucleus::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
testnucleus:: dep libcytobase.a libcytosim.a testnucleus.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -L$(LIBDIR) -lcytosim -lcytobase $(L) $(LGL) -o $@
testnucleus::
	@cp -f $@$(BINEXT) $(BINDIR)

testparam::
	@cd $(BASEDIR) && $(MAKE) parameter.o iomessages.o parameter_typed.o parameter_list.o
testparam:: dep testparam.o parameter.o iomessages.o parameter_typed.o parameter_list.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
testparam::
	@cp -f $@$(BINEXT) $(BINDIR)

testquaternion::
	@cd $(BASEDIR) && $(MAKE) quaternion.o matrix3.o random.o
testquaternion:: dep quaternion.o testquaternion.o matrix3.o random.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) $(L) -o $@
testquaternion::
	@cp -f $@$(BINEXT) $(BINDIR)

testrasterizer::
	@cd $(BASEDIR) && $(MAKE) random.o vecteur3.o quaternion.o rasterizer_visual.o
testrasterizer:: dep random.o rasterizer_visual.o vecteur3.o quaternion.o testrasterizer.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) $(LGL) -o $@
testrasterizer::
	@cp -f $@$(BINEXT) $(BINDIR)

testspace::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(PLAYDIR) && $(MAKE) glextensions.o
testspace:: dep libcytobase.a testspace.o glextensions.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -L$(LIBDIR) -lcytobase $(L) $(LGL) -o $@
testspace::
	@cp -f $@$(BINEXT) $(BINDIR)

testmap::
	@cd $(BASEDIR) && $(MAKE) random.o iomessages.o
testmap:: dep testmap.o random.o iomessages.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
testmap::
	@cp -f $@$(BINEXT) $(BINDIR)

testrandom::
	@cd $(BASEDIR) && $(MAKE) random.o smath.o
testrandom:: dep testrandom.o random.o smath.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
testrandom::
	@cp -f $@$(BINEXT) $(BINDIR)

testsphere::
	@cd $(BASEDIR) && $(MAKE) random.o pointsonsphere.o
testsphere:: dep testsphere.o random.o pointsonsphere.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
testsphere::
	@cp -f $@$(BINEXT) $(BINDIR)

testpolyroots::
	@cd $(BASEDIR) && $(MAKE) poly_roots.o
testpolyroots:: dep testpolyroots.o poly_roots.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
testpolyroots::
	@cp -f $@$(BINEXT) $(BINDIR)

#NEWS: 08/17/04: We now create libraries from all the objects in base and sim.
#Now some of the test-programs depend on libcytobase.a and libcytosim.a and we
#therefore need empty rules for these two libs in addition to the rules for
#the single object files.

libcytobase.a:

libcytosim.a:

#This rule tells make to do nothing if a base-, sim- or play-object is not
#there. But since some of the test-programs depend on $(BASEOBJS), $(SIMOBJS)
#and/or $(PLAYOBJS), they are remade if one of these objects changed for
#whatever reason.

$(BASEOBJS): 

rasterizer_visual.o:

$(SIMOBJS):

$(PLAYOBJS):


.PHONY: clean
clean:
	rm -f $(addprefix $(BUILDDIR)/, $(TESTOBJS)) core dep $(ALLTESTS)
#	rm -f $(BUILDDIR)/*.o core dep $(ALLTESTS)

.PHONY: sterile
sterile: clean
	@cd $(BASEDIR) && $(MAKE) clean
	@cd $(SIMDIR) && $(MAKE) clean
	@cd $(PLAYDIR) && $(MAKE) clean

dep: $(ALLTESTSOURCES) $(ALLTESTHEADERS)
	gcc -MM $(INCLUDE)            *.cc > dep


-include dep
