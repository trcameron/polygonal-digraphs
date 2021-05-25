#####################################################
#	Polygonal-Digraphs Make File					#
#####################################################
include make.inc

nrpolyn:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o nrpolyn c++/nrpolyn.cpp
	
nsearch9f:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o nsearch9f c++/nsearch9f.cpp
	
uninstall:
	@rm -f nrpolyn
	@rm -f nsearch9f