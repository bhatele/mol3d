include ../compile/Makefile.config

TARGET_LIB	= $(MOL3D)/compile/libCharm.a
OBJS		= Main.o Patch.o Compute.o
OBJS_DECL	= mol3d.decl.h
INC		= -I$(MAKE_MOL3D_PHYSICS) -I$(MAKE_MOL3D_NAMDSRC)

mol3d_charm: $(OBJS) $(OBJS_DECL)
	$(CHARMC) -o $(TARGET_LIB) $(OBJS) -language charm++ \
	-module CkMulticast $(OPTS)

Main.o: Main.C Main.h mol3d.decl.h $(MAKE_MOL3D_PHYSICS)/defs.h $(MAKE_MOL3D_PHYSICS)/sqrtTable.h $(MAKE_MOL3D_NAMDSRC)/PluginIOMgr.h $(MAKE_MOL3D_NAMDSRC)/ConfigList.h \
	$(MAKE_MOL3D_NAMDSRC)/Parameters.h $(MAKE_MOL3D_NAMDSRC)/SimParameters.h
	$(CHARMC) $(OPTS) $(INC) -o Main.o Main.C

Patch.o: Patch.C Patch.h mol3d.decl.h $(MAKE_MOL3D_PHYSICS)/defs.h 
	$(CHARMC) $(OPTS) $(INC) -o Patch.o Patch.C

mol3d.decl.h: mol3d.ci
	$(CHARMC) mol3d.ci

Compute.o: Compute.C Compute.h mol3d.decl.h $(MAKE_MOL3D_PHYSICS)/defs.h $(MAKE_MOL3D_PHYSICS)/nonbonded.h
	$(CHARMC) $(OPTS) $(INC) -o Compute.o Compute.C

clean:
	rm -f *.decl.h *.def.h *.o
