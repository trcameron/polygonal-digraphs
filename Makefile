#####################################################
#	Polygonal-Digraphs Make File					#
#####################################################
include make.inc

nrpolyn:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o nrpolyn c++/nrpolyn.cpp
	
nsearch9f:
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o nsearch9f c++/nsearch9f.cpp
	
run_nrpolyn:
	$(NAUTY)/geng $(N) | $(NAUTY)/directg | ./nrpolyn | tee d6files/polygonal$(N).d6
	
py_examples:
	$(PYTHON) python/examples.py
	
uninstall:
	@rm -f nrpolyn
	@rm -f nsearch9f