CC= cc
CFLAGS= -O -g -Wall
# CFLAGS= -g -Wall
ICFITSIO=-I/usr/include -I/usr/include/cfitsio
# ICFITSIO=$(ATLAS_HOME)/include
# LCFITSIO=-L/usr/lib -L/usr/lib64
# LCFITSIO=$(ATLAS_HOME)/lib

PDIR= ../tphot
IPSF= -I$(PDIR)
POBJ= $(PDIR)/nitfread.o $(PDIR)/rwfits.o $(PDIR)/psf2d.o $(PDIR)/fitlm.o $(PDIR)/linsolve.o $(PDIR)/sortc.o


all: strker

.c.o: 
	$(CC) $(CFLAGS) $(ICFITSIO) -I$(PDIR) -c $<

strker: strker.o tpfn.o 
	$(CC) $(CFLAGS) $(LCFITSIO) strker.o tpfn.o /atlas/src/trunk/util/sortlib/tsort.o $(POBJ) -lcfitsio -lm -o strker

strker06: strker06.o tpfn.o 
	$(CC) $(CFLAGS) $(LCFITSIO) strker06.o tpfn.o /atlas/src/trunk/util/sortlib/tsort.o $(POBJ) -lcfitsio -lm -o strker06

endfindtest: endfindtest.o tpfn.o 
	$(CC) $(CFLAGS) $(LCFITSIO) endfindtest.o tpfn.o /atlas/src/trunk/util/sortlib/tsort.o $(POBJ) -lcfitsio -lm -o endfindtest

tpfn: tpfn.o
	$(CC) $(CFLAGS) $(LCFITSIO) tpfn.o $(POBJ) -lcfitsio -lm -o tpfn

# Vestigial
strkonv: strkonv.o
	$(CC) $(CFLAGS) $(LCFITSIO) strkonv.o -lcfitsio -lm -o strkonv

install:
	install -m 0775 strker $(ATLAS_HOME)/bin/
#	install -m 0775 strker06 $(ATLAS_HOME)/bin/
#	install -m 0775 endfindtest $(ATLAS_HOME)/bin/
	install -m 0664 strker.man $(ATLAS_HOME)/man/man1/strker.1

clean:
	rm -f *.o strker


