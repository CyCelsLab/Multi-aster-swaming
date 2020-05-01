#RCS: $Id: makefile,v 2.8 2005/02/17 16:24:30 foethke Exp $
#makefile for the tools, Francois Nedelec 2000, Dietrich Foethke 2004.

ROOT=../..
SRCROOT=..
BASEDIR=$(SRCROOT)/base
SIMDIR=$(SRCROOT)/sim
PLAYDIR=$(SRCROOT)/play
TOOLDIR=.
BUILDDIR=../../build
BINDIR=../../bin
LIBDIR=../../lib


INCLUDE=-I$(BASEDIR) -I$(SIMDIR) -I$(PLAYDIR)

include $(ROOT)/makefile.mk



#-----------------set the path for the prerequisites for .o files--------------

vpath %.o $(BUILDDIR)
vpath %.a $(LIBDIR)

#-----------------lists of all the files for the test-programs-----------------

TOOLHEADERS=$(wildcard $(TOOLDIR)/*.h)
TOOLSOURCES=$(wildcard $(TOOLDIR)/*.cc)

ALLTOOLHEADERS=$(ALLPLAYHEADERS) $(TOOLHEADERS)
ALLTOOLSOURCES=$(ALLPLAYSOURCES) $(TOOLSOURCES)

#----------patterns to grep the object and lib files from the prerequisites----

OBJPATTERN=$(patsubst %.o,$(BUILDDIR)/%.o, $(filter %.o, $(notdir $^)))
LIBPATTERN=$(patsubst lib%.a,-l%, $(filter %.a, $(notdir $^)))
COMPILER_COMMAND=$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
COMPILER_COMMAND_L=$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -L$(LIBDIR) $(LIBPATTERN) $(L) -o $@

#-----------------rules--------------------------------------------------------


#The variable $(ALLTOOLS) is set in $(ROOT)/makefile.mk.
#If you add a target in this makefile, please add it to $(ALLTOOLS), too!


.PHONY: all
all: $(ALLTOOLS)



frametool: dep frametool.o
	$(COMPILER_COMMAND)
	@cp -f $@ $(BINDIR)

readwrite::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
readwrite:: dep libcytosim.a libcytobase.a readwrite.o
	$(COMPILER_COMMAND_L)
readwrite::
	@cp -f $@$(BINEXT) $(BINDIR)


createstart::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
createstart:: dep libcytosim.a libcytobase.a createstart.o
	$(COMPILER_COMMAND_L)
createstart::
	@cp -f $@$(BINEXT) $(BINDIR)


meanlength::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
meanlength:: dep libcytosim.a libcytobase.a meanlength.o
	$(COMPILER_COMMAND_L)
meanlength::
	@cp -f $@$(BINEXT) $(BINDIR)

readwritedata::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o
readwritedata:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o readwritedata.o
	$(COMPILER_COMMAND)
readwritedata::
	@cp -f $@$(BINEXT) $(BINDIR)

reducedata::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o
reducedata:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o reducedata.o random.o
	$(COMPILER_COMMAND)
reducedata::
	@cp -f $@$(BINEXT) $(BINDIR)

analyse0::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
	@cd $(PLAYDIR) && $(MAKE) reader.o
analyse0:: dep libcytosim.a libcytobase.a reader.o analyse0.o
	$(COMPILER_COMMAND_L)
analyse0::
	@cp -f $@$(BINEXT) $(BINDIR)


analyse1::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
	@cd $(PLAYDIR) && $(MAKE) reader.o
analyse1:: dep libcytosim.a libcytobase.a reader.o analyse1.o
	$(COMPILER_COMMAND_L)
analyse1::
	@cp -f $@$(BINEXT) $(BINDIR)

analyse2::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
	@cd $(PLAYDIR) && $(MAKE) reader.o
analyse2:: dep libcytosim.a libcytobase.a reader.o analyse2.o
	$(COMPILER_COMMAND_L)
analyse2::
	@cp -f $@$(BINEXT) $(BINDIR)

analyse3::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
	@cd $(PLAYDIR) && $(MAKE) reader.o
analyse3:: dep libcytosim.a libcytobase.a reader.o analyse3.o
	$(COMPILER_COMMAND_L)
analyse3::
	@cp -f $@$(BINEXT) $(BINDIR)

analyse4::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
	@cd $(PLAYDIR) && $(MAKE) reader.o
analyse4:: dep libcytosim.a libcytobase.a reader.o analyse4.o
	$(COMPILER_COMMAND_L)
analyse4::
	@cp -f $@$(BINEXT) $(BINDIR)

analyse_nucpos::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
	@cd $(PLAYDIR) && $(MAKE) reader.o
analyse_nucpos:: dep libcytosim.a libcytobase.a reader.o analyse_nucpos.o
	$(COMPILER_COMMAND_L)
analyse_nucpos::
	@cp -f $@$(BINEXT) $(BINDIR)

datagen::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o iomessages.o random.o
	@cd $(SIMDIR) && $(MAKE) sim_param.o
datagen:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o datagen.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
datagen::
	@cp -f $@$(BINEXT) $(BINDIR)

datagen1::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o iomessages.o random.o
	@cd $(SIMDIR) && $(MAKE) sim_param.o
datagen1:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o datagen1.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
datagen1::
	@cp -f $@$(BINEXT) $(BINDIR)

datagen2::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o iomessages.o random.o
	@cd $(SIMDIR) && $(MAKE) sim_param.o
datagen2:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o datagen2.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
datagen2::
	@cp -f $@$(BINEXT) $(BINDIR)

datagen3::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o iomessages.o random.o
	@cd $(SIMDIR) && $(MAKE) sim_param.o
datagen3:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o datagen3.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
datagen3::
	@cp -f $@$(BINEXT) $(BINDIR)

datagen4::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o iomessages.o random.o
	@cd $(SIMDIR) && $(MAKE) sim_param.o
datagen4:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o datagen4.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
datagen4::
	@cp -f $@$(BINEXT) $(BINDIR)

datagen5::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o iomessages.o random.o
	@cd $(SIMDIR) && $(MAKE) sim_param.o
datagen5:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o datagen5.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
datagen5::
	@cp -f $@$(BINEXT) $(BINDIR)

datagen6::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o iomessages.o random.o
	@cd $(SIMDIR) && $(MAKE) sim_param.o
datagen6:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o datagen6.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
datagen6::
	@cp -f $@$(BINEXT) $(BINDIR)

datagen7::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o iomessages.o random.o
	@cd $(SIMDIR) && $(MAKE) sim_param.o
datagen7:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o datagen7.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
datagen7::
	@cp -f $@$(BINEXT) $(BINDIR)

datagen8::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o iomessages.o random.o
	@cd $(SIMDIR) && $(MAKE) sim_param.o
datagen8:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o datagen8.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
datagen8::
	@cp -f $@ $(BINDIR)

datagen9::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o iomessages.o random.o
	@cd $(SIMDIR) && $(MAKE) sim_param.o
datagen9:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o datagen9.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
datagen9::
	@cp -f $@ $(BINDIR)

#this is for future use
#please add datagen10 to $(ALLTOOLS) and datagen10.o to $(TOOLOBJS)
#in $(SRCROOT)/makefile.mk if you want to use it
datagen10::
	@cd $(BASEDIR) && $(MAKE) parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o
datagen10:: dep parameter.o parameter_typed.o parameter_list.o sim_param.o iomessages.o random.o datagen10.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(OBJPATTERN) -o $@
datagen10::
	@cp -f $@ $(BINDIR)


#NEWS: 08/17/04: We now create libraries from all the objects in base and sim.
#Now some of the tools depend on libcytobase.a and libcytosim.a and we
#therefore need empty rules for these two libs in addition to the rules for
#the single object files.

libcytobase.a:

libcytosim.a:

#These rules tell make to do nothing if a base-, sim- or play-object is not
#there. But since some of the tools depend on $(BASEOBJS), $(SIMOBJS) and\or
#$(PLAYOBJS), they are remade if one of these objects changed for whatever
#reason.

$(BASEOBJS):

$(SIMOBJS):

$(PLAYOBJS):


.PHONY: clean
clean:
	rm -f $(addprefix $(BUILDDIR)/, $(TOOLOBJS)) core dep $(ALLTOOLS)
#	rm -f $(BUILDDIR)/*.o core dep $(ALLTOOLS)

.PHONY: sterile
sterile: clean
	@cd $(BASEDIR) && $(MAKE) clean
	@cd $(SIMDIR) && $(MAKE) clean
	@cd $(PLAYDIR) && $(MAKE) clean

dep: $(ALLTOOLSOURCES) $(ALLTOOLHEADERS)
	gcc -MM $(INCLUDE)            *.cc > dep


-include dep
