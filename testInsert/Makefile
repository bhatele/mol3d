OPTS=-g
CHARMC=~/charm/bin/charmc $(OPTS)

default: testInsert

testInsert: testInsert.o 
	  $(CHARMC) -language charm++ -o testInsert testInsert.o

testInsert.o: testInsert.C testInsert.h testInsert.decl.h testInsert.def.h
	  $(CHARMC) -language charm++ -o testInsert.o testInsert.C

testInsert.decl.h testInsert.def.h: testInsert.ci
	  $(CHARMC) testInsert.ci

clean:
	rm -f testInsert charmrun testInsert.o testInsert.decl.h testInsert.def.h


