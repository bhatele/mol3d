include ../compile/Makefile.config

TARGET_LIB = $(MOL3D)/compile/libNamdsrc.a
OBJS	   = jsplugin.o PluginIOMgr.o ConfigList.o strlib.o common.o \
	     GromacsTopFile.o MStream.o Parameters.o Communicate.o \
	     ProcessorPrivate.o InfoStream.o parm.o
OBJS_DECL  =

mol3d_namdsrc: $(OBJS) $(OBJS_DECL)
	$(CHARMC) -o $(TARGET_LIB) $(OBJS) -language charm++ $(OPTS)

jsplugin.o: jsplugin.c largefiles.h hash.h hash.c fastio.h endianswap.h molfile_plugin.h vmdplugin.h
	$(CHARMC) $(OPTS) -DVMDPLUGIN=molfile_jsplugin -o jsplugin.o jsplugin.c

PluginIOMgr.o: PluginIOMgr.C PluginIOMgr.h molfile_plugin.h vmdplugin.h libmolfile_plugin.h
	$(CHARMC) $(OPTS) -o PluginIOMgr.o PluginIOMgr.C

ConfigList.o: ConfigList.C InfoStream.h ConfigList.h common.h strlib.h
	$(CHARMC) $(OPTS) -o ConfigList.o ConfigList.C

strlib.o: strlib.C strlib.h common.h
	$(CHARMC) $(OPTS) -o strlib.o strlib.C

common.o: common.C common.h InfoStream.h
	$(CHARMC) $(OPTS) -o common.o common.C

GromacsTopFile.o: GromacsTopFile.C common.h ResizeArray.h ResizeArrayRaw.h GromacsTopFile.h
	$(CHARMC) $(OPTS) -o GromacsTopFile.o GromacsTopFile.C

MStream.o: MStream.C Communicate.h MStream.h Vector.h common.h Debug.h
	$(CHARMC) $(OPTS) -o MStream.o MStream.C

Parameters.o: Parameters.C InfoStream.h Parameters.h parm.h common.h structures.h strlib.h MStream.h Vector.h SimParameters.h Lattice.h NamdTypes.h \
	      ResizeArray.h ResizeArrayRaw.h Tensor.h MGridforceParams.h GromacsTopFile.h Communicate.h ConfigList.h Debug.h
	$(CHARMC) $(OPTS) -o Parameters.o Parameters.C

Communicate.o: Communicate.C Communicate.h MStream.h Vector.h common.h
	$(CHARMC) $(OPTS) -o Communicate.o Communicate.C

ProcessorPrivate.o: ProcessorPrivate.C ProcessorPrivate.h BOCgroup.h Debug.h
	$(CHARMC) $(OPTS) -o ProcessorPrivate.o ProcessorPrivate.C

InfoStream.o: InfoStream.C InfoStream.h Vector.h common.h Tensor.h
	$(CHARMC) $(OPTS) -o InfoStream.o InfoStream.C

parm.o: parm.C strlib.h common.h InfoStream.h parm.h
	$(CHARMC) $(OPTS) -o parm.o parm.C

clean:
	rm -f *.decl.h *.def.h *.o
