export NAMD="$HOME/work/namd2"
export SRC="$NAMD/src"
export PLUGSRC="$NAMD/plugins/molfile_plugin/src"
export PLUGINC="$NAMD/plugins/include"

cp $SRC/BOCgroup.h .
cp $SRC/common.h .
cp $SRC/common.C .
cp $SRC/Communicate.C .
cp $SRC/Communicate.h .
cp $SRC/ConfigList.C .
cp $SRC/ConfigList.h .
cp $SRC/Debug.h .
cp $PLUGSRC/endianswap.h .
cp $PLUGSRC/fastio.h .
cp $PLUGSRC/fortread.h .
cp $SRC/GromacsTopFile.C .
cp $SRC/GromacsTopFile.h .
cp $PLUGSRC/hash.c .
cp $PLUGSRC/hash.h .
cp $SRC/InfoStream.C .
cp $SRC/InfoStream.h .
#cp $PLUGSRC/jsplugin.c .
cp $PLUGSRC/largefiles.h .
cp $PLUGINC/libmolfile_plugin.h .
cp $SRC/Lattice.h .
cp $SRC/MGridforceParams.h .
cp $PLUGINC/molfile_plugin.h .
cp $SRC/MStream.C .
cp $SRC/MStream.h .
cp $SRC/NamdTypes.h .
#cp $SRC/Parameters.C .
cp $SRC/Parameters.h .
cp $SRC/parm.h .
cp $SRC/parm.C .
cp $PLUGSRC/periodic_table.h .
cp $SRC/PluginIOMgr.C .
cp $SRC/PluginIOMgr.h .
cp $SRC/ProcessorPrivate.C .
cp $SRC/ProcessorPrivate.h .
cp $PLUGSRC/readpdb.h .
cp $SRC/ResizeArray.h .
cp $SRC/ResizeArrayRaw.h .
cp $SRC/strlib.C .
cp $SRC/strlib.h .
cp $SRC/structures.h .
cp $SRC/Tensor.h .
cp $SRC/Vector.h .
cp $PLUGINC/vmdplugin.h .

