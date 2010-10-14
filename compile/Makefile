include Makefile.config

#===============================================================================
all: full_code

#===============================================================================
full_code: charm namdsrc
	cd $(MOL3D)/compile; \
	$(CHARMC) -module CommonLBs -o $(BINARY) \
		-module CkMulticast -language charm++ \
		-L. -lCharm -lNamdsrc -lck

#===============================================================================
charm:
	(cd $(MAKE_MOL3D_CHARM); make OPTS="$(OPTS) -DUSE_SECTION_MULTICAST";)

namdsrc:
	(cd $(MAKE_MOL3D_NAMDSRC); make OPTS="$(OPTS)";)
#===============================================================================

clean: clean_charm clean_namdsrc
	rm -f lib*.a charmrun $(BINARY)

clean_charm:
	(cd $(MAKE_MOL3D_CHARM); make clean;)

clean_namdsrc:
	(cd $(MAKE_MOL3D_NAMDSRC); make clean;)
