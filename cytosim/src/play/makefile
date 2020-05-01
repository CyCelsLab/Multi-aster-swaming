#RCS: $Id: makefile,v 2.2 2005/02/17 16:24:30 foethke Exp $
#makefile for play, Francois Nedelec 2000, Dietrich Foethke 2004.

ROOT=../..
SRCROOT=..
BASEDIR=$(SRCROOT)/base
SIMDIR=$(SRCROOT)/sim
PLAYDIR=.
BUILDDIR=../../build
BINDIR=../../bin
LIBDIR=../../lib


INCLUDE=-I$(BASEDIR) -I$(SIMDIR) -I$(PLAYDIR)

include $(ROOT)/makefile.mk



PLAYOBJSPATH=$(PLAYOBJS:%.o=$(BUILDDIR)/%.o)


#-----------------set the path for the prerequisites for .o files--------------

vpath %.o $(BUILDDIR)
vpath %.a $(LIBDIR)

#-----------------rules--------------------------------------------------------


#If you add a target in this makefile, please add it to $(ALLPLAYPROGS) in
#$(ROOT)/makefile.mk, too!

play::
	@cd $(BASEDIR) && $(MAKE)
	@cd $(SIMDIR) && $(MAKE) libcytosim.a
play:: dep libcytobase.a libcytosim.a $(PLAYOBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PLAYOBJSPATH) -L$(LIBDIR) -lcytosim -lcytobase $(L) $(LGL) -o play
play::
	@cp -f $@$(BINEXT) $(BINDIR)


#NEWS: 08/17/04: We now create libraries from all the objects in base and sim.
#Now play only depends on libcytobase.a and libcytosim.a and we therefore
#need empty rules for these two libs rather then all the single object files.

libcytobase.a:

libcytosim.a:

#This rule tells make to do nothing if a base- or sim-object is not there.
#But since play depends on $(BASEOBJS) and $(SIMOBJS), play is remade if one of
#these objects changed for whatever reason.

#$(BASEOBJS):

#$(SIMOBJS):


.PHONY: clean
clean:
	rm -f $(addprefix $(BUILDDIR)/, $(PLAYOBJS)) core dep play
#	rm -f $(BUILDDIR)/*.o $(LIBDIR)/*.a core dep play

.PHONY: sterile
sterile: clean
	@cd $(BASEDIR) && $(MAKE) clean
	@cd $(SIMDIR) && $(MAKE) clean

dep: $(ALLPLAYSOURCES) $(ALLPLAYHEADERS)
	gcc -MM $(INCLUDE)            *.cc > dep

-include dep
