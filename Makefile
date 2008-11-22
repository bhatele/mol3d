CHARMC		= ../../../bin/charmc 
OPTS		= -O3  -DUSE_SECTION_MULTICAST -DRUN_LIVEVIZ

all: mol3d

projections: mol3d.prj

mol3d: Patch.o Compute.o
	$(CHARMC) $(OPTS) -module CkMulticast -language charm++ -o mol3d Patch.o \
	Compute.o

mol3d.prj: Patch.o Compute.o
	$(CHARMC) $(OPTS) -language charm++ -tracemode projections -o mol3d.prj Patch.o Compute.o

Patch.o: Patch.C Patch.h Patch.decl.h common.h
	$(CHARMC) $(OPTS) -o Patch.o Patch.C

Patch.decl.h: Patch.ci
	$(CHARMC) Patch.ci

Compute.o: Compute.C Compute.h Patch.decl.h common.h
	$(CHARMC) $(OPTS) -o Compute.o Compute.C

clean:
	rm -f *.decl.h *.def.h *.o mol3d mol3d.prj charmrun
