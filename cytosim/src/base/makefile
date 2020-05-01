#RCS: $Id: makefile,v 2.2 2005/02/17 16:24:30 foethke Exp $
#makefile for the base-objects, Francois Nedelec 2000, Dietrich Foethke 2004.

ROOT=../..
SRCROOT=..
BASEDIR=.
SIMDIR=$(SRCROOT)/sim
PLAYDIR=$(SRCROOT)/play
BUILDDIR=../../build
LIBDIR=../../lib


INCLUDE=-I$(BASEDIR)

include $(ROOT)/makefile.mk


#-----------------set the path for the prerequisites for .o and .a files--------------

vpath %.o $(BUILDDIR)
vpath %.a $(LIBDIR)

#-----------------rules--------------------------------------------------------


.PHONY: all
all: baseobjs libcytobase.a

.PHONY: baseobjs
baseobjs: dep $(BASEOBJS)


rasterizer_visual.o: rasterizer.cc rasterizer.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DVISUAL_DEBUG -c rasterizer.cc -o $(BUILDDIR)/rasterizer_visual.o


dr_burgersm_rock.o: dr_burgersm_rock.f 

rock2.o: rock2.f

#make a library out of all base object files
libcytobase.a: $(BASEOBJS)
	ar rcs $(LIBDIR)/libcytobase.a $(addprefix $(BUILDDIR)/, $(notdir $^))


.PHONY: clean
clean:
	rm -f $(addprefix $(BUILDDIR)/, $(BASEOBJS) rasterizer_visual.o) $(LIBDIR)/libcytobase.a core dep
#	rm -f $(BUILDDIR)/*.o $(LIBDIR)/*.a core dep


dep: $(BASESOURCES) $(BASEHEADERS)
	gcc -MM $(INCLUDE) *.cc > dep
#	gcc -MM            *.f  >> dep



-include dep
