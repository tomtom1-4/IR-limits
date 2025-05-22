# default optimization level, change for debugging
LEVEL = 0

CXX = c++
FC  = gfortran

LHADIR = $(shell lhapdf-config --prefix)

USE_QD   = true

ifeq ($(FC),ifort)
CXXFLAGS = -O$(LEVEL) -std=c++17 -I.. -I$(LHADIR)/include $(OPT) \
           -I /home/tom/Documents/software/software/Stripper/Stripper/
else
CXXFLAGS = -O$(LEVEL) -std=c++17 -I.. -I$(LHADIR)/include -DGFORTRAN $(OPT) \
           -I /home/tom/Documents/software/software/Stripper/Stripper/
endif
MAKE_SHR = true
ifeq ($(strip $(USE_OL)),true)
CXXFLAGS += -DUSE_OL
endif
ifeq ($(strip $(USE_QD)),true)
CXXFLAGS += -DUSE_QD
endif
ifeq ($(strip $(MAKE_SHR)),true)
CXXFLAGS += -fPIC -shared
endif
CXXFLAGS += -DGNU_NAME_MANGLING

ifeq ($(shell uname), Darwin)
  LD      = $(CXX)
  LDFLAGS =
  LIBS    = -L.. -lstripper -L../virtual/hardcoded -lhardcoded -L../AvH -lavh \
            -L$(LHADIR)/lib -lLHAPDF
  ifeq ($(FC),gfortran)
  LIBS += -lgfortran -lquadmath
  endif
else
  ifeq ($(FC),ifort)
  LDFLAGS = -nofor-main
  endif
  LD      = $(FC)
  LIBS    = -L.. -lstripper -L../virtual/hardcoded -lhardcoded -L../AvH -lavh \
            -L$(LHADIR)/lib -lLHAPDF -lstdc++ -ldl
endif

LIBS += -L/home/tom/Documents/software/software/recola_otter/recola2-2.2.3/ -lrecola

ifeq ($(strip $(USE_OL)),true)
LIBS += -lopenloops
endif
ifeq ($(strip $(USE_NJET)),true)
LIBS += -lnjet2
endif
ifeq ($(strip $(USE_CLN)),true)
LIBS += -lcln
endif
ifeq ($(strip $(USE_GNC)),true)
LIBS += -lginac
endif
ifeq ($(strip $(USE_FIVE)),true)
LIBS += -lPentagonFunctions -lLi2++ -lqd
endif
ifeq ($(strip $(USE_FEXT)),true)
LIBS += -lFivePointAmplitudes
endif
ifeq ($(strip $(USE_QD)),true)
LIBS += -lqd
endif
ifeq ($(strip $(USE_QCDL)),true)
LIBS += -lqcdloop
endif

SRCS = src/amplitude.cpp src/color.cpp src/eikonal.cpp src/tinyxml2.cpp src/splitting_functions.cpp src/utilities.cpp \
	/home/tom/Uni/phd/NnLO/phase_space_parameterization/src/Tree.cpp \
  /home/tom/Uni/phd/NnLO/phase_space_parameterization/src/PhaseSpace.cpp \
	/home/tom/Uni/phd/NnLO/phase_space_parameterization/src/Utilities.cpp \
  src/splitting_functions/P0_qqPqPbarg.cpp \
  src/splitting_functions/P0_qqqbarg.cpp \
  src/splitting_functions/P0_qggg.cpp \
  src/splitting_functions/P0_qqbarqPqPbar.cpp \
  src/splitting_functions/P0_qqbarqqbar.cpp \
  src/splitting_functions/P0_ggqqbar.cpp \
  src/splitting_functions/P0_gggg.cpp \
  src/splitting_functions/P1_qqPqPbar.cpp \
  src/splitting_functions/P1_qqqbar.cpp \
  src/splitting_functions/P1_qgg.cpp \
  src/splitting_functions/P1_gqqbar.cpp \
  src/splitting_functions/P1_ggg.cpp \
  src/splitting_functions/P2_qg.o \

OBJS = $(patsubst %.cpp,%.o,$(SRCS))

.PHONY: clean

.PRECIOUS: src/amplitude.o src/color.o src/eikonal.o src/tinyxml2.o src/splitting_functions.o src/utilities.o \
  src/splitting_functions/P0_qqPqPbarg.o \
  src/splitting_functions/P0_qqqbarg.o \
  src/splitting_functions/P0_qggg.o \
  src/splitting_functions/P0_qqbarqPqPbar.o \
  src/splitting_functions/P0_qqbarqqbar.o \
  src/splitting_functions/P0_ggqqbar.o \
  src/splitting_functions/P0_gggg.o \
  src/splitting_functions/P1_qqPqPbar.o \
  src/splitting_functions/P1_qqqbar.o \
  src/splitting_functions/P1_qgg.o \
  src/splitting_functions/P1_gqqbar.o \
  src/splitting_functions/P1_ggg.o \
  src/splitting_functions/P2_qg.o \
  /home/tom/Uni/phd/NnLO/phase_space_parameterization/src/Tree.o \
  /home/tom/Uni/phd/NnLO/phase_space_parameterization/src/PhaseSpace.o \
  /home/tom/Uni/phd/NnLO/phase_space_parameterization/src/Utilities.o

#target: $(OBJS)

%.exe: %.o $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) -DGNU_NAME_MANGLING -lpthread

# Compilation rule
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -DGNU_NAME_MANGLING $< -o $@

clean:
	rm -f *.o $(EXEC) $(OBJS) *.exe

