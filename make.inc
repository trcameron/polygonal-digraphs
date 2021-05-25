#####################################################
#   Polygonal-Digraphs Make Include File    		#
#####################################################
NAUTY = /usr/local/Cellar/nauty/27r1/bin
PYTHON = python3
N = 3
# ---------------------------------------------------------------------
# Compiler selection 
# --------------------------------------------------------------------- 
CCC = g++ -std=gnu++2a
# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
CCOPT = -Ofast -march=native -mtune=native -DMIN_VERTS=2 -DMAX_VERTS=8 -DNO_DYNAMIC_SIZE -DNDEBUG -Wno-deprecated-anon-enum-enum-conversion
# ---------------------------------------------------------------------
# Include options
# ---------------------------------------------------------------------
BOOSTINCDIR = /usr/local/Cellar/boost/1.75.0_3/include
EIGENINCDIR = /Users/thcameron/eigen-3.4-rc1
CCFLAGS = $(CCOPT) -I$(BOOSTINCDIR) -I$(EIGENINCDIR)