DIAGENIDSOLVERDIR =/home/tom/Documents/software/software/DiaGenIdSolver
UTILITIESDIR      = $(DIAGENIDSOLVERDIR)/Utilities
DIAGENDIR         = $(DIAGENIDSOLVERDIR)/DiaGen

CXX      = c++
CXXFLAGS = -std=c++11 -I$(DIAGENDIR) -I$(UTILITIESDIR) -O
LD       = c++
LDFLAGS  = -L$(DIAGENDIR) -L$(UTILITIESDIR)

LIBS = -lDiaGen -lUtilities

.PHONY: clean

diags: diags.o
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f diags *.o *.out *~
