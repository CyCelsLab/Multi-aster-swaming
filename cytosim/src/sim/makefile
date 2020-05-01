#RCS: $Id: makefile,v 2.2 2005/02/17 16:24:30 foethke Exp $
#makefile for sim, Francois Nedelec 2000, Dietrich Foethke 2004.

ROOT=../..
SRCROOT=..
BASEDIR=$(SRCROOT)/base
SIMDIR=.
PLAYDIR=$(SRCROOT)/play
BUILDDIR=../../build
BINDIR=../../bin
LIBDIR=../../lib


INCLUDE=-I$(BASEDIR) -I$(SIMDIR)

include $(ROOT)/makefile.mk



#-----------------set the path for the prerequisites for .o files--------------

vpath %.o $(BUILDDIR)
vpath %.a $(LIBDIR)

#-----------------rules--------------------------------------------------------


#If you add a target in this makefile, please add it to $(ALLSIMPROGS) in
#$(ROOT)/makefile.mk, too!

sim::
	@cd $(BASEDIR) && $(MAKE)
sim:: dep libcytosim.a libcytobase.a main.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BUILDDIR)/main.o -L$(LIBDIR) -lcytosim -lcytobase $(L) -o sim
sim::
	@cp -f $@$(BINEXT) $(BINDIR)



sim1::
	@cd $(BASEDIR) && $(MAKE)
sim1:: dep libcytosim.a libcytobase.a main.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BUILDDIR)/main.o -L$(LIBDIR) -lcytosim -lcytobase $(L) -o sim1
sim1:: 
	@cp -f $@$(BINEXT) $(BINDIR)



sim2::
	@cd $(BASEDIR) && $(MAKE)
sim2:: dep libcytosim.a libcytobase.a main.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(BUILDDIR)/main.o -L$(LIBDIR) -lcytosim -lcytobase $(L) -o sim2
sim2::
	@cp -f $@$(BINEXT) $(BINDIR)


#NEWS: 08/17/04: We now create libraries from all the objects in base and sim.
#Now sim only depends on libcytobase.a and libcytosim.a and we therefore
#need an empty rule for libcytobase.a. rather then all the single $(BASEOBJS).

libcytobase.a:

#This rule tells make to do nothing if a base-object is not there.
#But since sim depends on $(BASEOBJS), sim is remade if one of
#these objects changed for whatever reason.

#$(BASEOBJS):



#make a library out of all sim object files
libcytosim.a: $(SIMOBJS)
	ar rcs $(LIBDIR)/libcytosim.a $(addprefix $(BUILDDIR)/, $(notdir $^))



.PHONY: clean
clean:
	rm -f $(addprefix $(BUILDDIR)/, $(SIMOBJS) main.o) $(LIBDIR)/libcytosim.a core dep sim sim1 sim2
#	rm -f $(BUILDDIR)/*.o $(LIBDIR)/*.a core dep sim sim1 sim2

.PHONY: sterile
sterile: clean
	@cd $(BASEDIR) && $(MAKE) clean


dep: $(ALLSIMSOURCES) $(ALLSIMHEADERS)
	gcc -MM $(INCLUDE)            *.cc > dep


-include dep
