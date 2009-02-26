CHARMC		= $(HOME)/charm/bin/charmc 
OPTS		= -O3  -DUSE_SECTION_MULTICAST 

all: mol3d

projections: mol3d.prj

mol3d: Main.o Patch.o Compute.o ConfigList.o 
	$(CHARMC) $(OPTS) -module CkMulticast -language charm++ -o mol3d Main.o Patch.o \
	Compute.o dcdplugin.o PluginIOMgr.o jsplugin.o ConfigList.o 

mol3d.prj: Main.o Patch.o Compute.o ConfigList.o
	$(CHARMC) $(OPTS) -module CkMulticast -tracemode projections -language charm++ -o mol3d.prj Main.o Patch.o \
	Compute.o dcdplugin.o PluginIOMgr.o jsplugin.o ConfigList.o

Main.o: Main.C Main.h Main.decl.h common.h
	$(CHARMC) $(OPTS) -o Main.o Main.C

Patch.o: Patch.C Patch.h Main.decl.h common.h
	$(CHARMC) $(OPTS) -o Patch.o Patch.C

Main.decl.h: Main.ci
	$(CHARMC) Main.ci

Compute.o: Compute.C Compute.h Main.decl.h common.h
	$(CHARMC) $(OPTS) -o Compute.o Compute.C

ConfigList.o: ConfigList.C ConfigList.h
	$(CHARMC) $(OPTS) -o ConfigList.o ConfigList.C

clean:
	rm -f *.decl.h *.def.h Main.o Patch.o Compute.o ConfigList.o mol3d mol3d.prj charmrun
