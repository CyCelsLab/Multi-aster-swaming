#RCS: $Id: makefile,v 2.7 2005/02/17 16:24:30 foethke Exp $
#top-level makefile for cytosim, Francois Nedelec 2000, Dietrich Foethke 2003.

ROOT=.
SRCROOT=src
BASEDIR=$(SRCROOT)/base
SIMDIR=$(SRCROOT)/sim
PLAYDIR=$(SRCROOT)/play
TESTDIR=$(SRCROOT)/test
TOOLDIR=$(SRCROOT)/tools
BUILDDIR=build
LIBDIR=lib
BINDIR=bin
DOCDIR=doc

FULLSRCROOT=$(CURDIR)/src
FULLBASEDIR=$(FULLSRCROOT)/base
FULLSIMDIR=$(FULLSRCROOT)/sim
FULLPLAYDIR=$(FULLSRCROOT)/play
FULLTOOLDIR=$(FULLSRCROOT)/tools
FULLBUILDDIR=$(CURDIR)/build
FULLBINDIR=$(CURDIR)/bin


include $(ROOT)/makefile.mk


#-----------------rules--------------------------------------------------------


#-----------------binary programs----------------------------------------------

.PHONY: sim
sim:
	@cd $(SIMDIR) && $(MAKE)

.PHONY: sim1
sim1:
	@cd $(SIMDIR) && $(MAKE) sim1

.PHONY: sim2
sim2:
	@cd $(SIMDIR) && $(MAKE) sim2

.PHONY: play
play:
	@cd $(PLAYDIR) && $(MAKE)

.PHONY: test
test:
	@cd $(TESTDIR) && $(MAKE)

.PHONY: tools
tools:
	@cd $(TOOLDIR) && $(MAKE)

.PHONY: base
base:
	@cd $(BASEDIR) && $(MAKE)

.PHONY: all
all: sim play test tools


.PHONY: $(ALLTESTS)
$(ALLTESTS):
	@cd $(TESTDIR) && $(MAKE) $@

.PHONY: $(ALLTOOLS)
$(ALLTOOLS):
	@cd $(TOOLDIR) && $(MAKE) $@


#-----------------single object files------------------------------------------


.PHONY: $(BASEOBJS)
$(BASEOBJS):
	@cd $(BASEDIR) && $(MAKE) $@

.PHONY: $(SIMOBJS)
$(SIMOBJS):
	@cd $(SIMDIR) && $(MAKE) $@

.PHONY: $(TESTOBJS)
$(TESTOBJS):
	@cd $(TESTDIR) && $(MAKE) $@

.PHONY: $(TOOLOBJS)
$(TOOLOBJS):
	@cd $(TOOLDIR) && $(MAKE) $@


#-----------------maintenance--------------------------------------------------  


.PHONY: tar
tar: clean binclean
	tar cf cytosim.tar --exclude \*.out --exclude \*tar\* \
        --exclude $(DOCDIR)/\* --exclude cytosim \
        --exclude $(BUILDDIR)/\* --exclude CVS *

.PHONY: tarmac
tarmac: clean binclean
	tar cf cytosim_mac.tar --exclude \*.out --exclude \*tar\* \
        --exclude $(DOCDIR)/\* --exclude cytosim/build *

.PHONY: tarsrc
tarsrc: clean
	tar cf cytosim_src.tar src makefile makefile.mk --exclude CVS

.PHONY: tarbin
tarbin: all
	tar cf cytosim_bin.tar --exclude \*.out --exclude \*tar\* \
        --exclude $(DOCDIR)/\* --exclude cytosim --exclude CVS *

.PHONY: tarbz2
tarbz2: clean binclean
	tar cjf cytosim.tar.bz2 --exclude \*.out --exclude \*tar\* \
        --exclude $(DOCDIR)/\* --exclude cytosim --exclude CVS *

.PHONY: clean
clean:
	@cd $(BASEDIR) && $(MAKE) clean
	@cd $(SIMDIR)  && $(MAKE) clean
	@cd $(PLAYDIR) && $(MAKE) clean
	@cd $(TESTDIR) && $(MAKE) clean
	@cd $(TOOLDIR) && $(MAKE) clean

.PHONY: binclean
binclean:
	cd $(BINDIR); \rm -f $(EVERYTHING)

.PHONY: sterile
sterile: clean binclean
	rm -f cytosim.tar cytosim_mac.tar cytosim_bin.tar cytosim.tar.bz2
	rm -f $(BUILDDIR)/*.o $(LIBDIR)/*.a
