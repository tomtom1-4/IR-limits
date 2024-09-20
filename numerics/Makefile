# default optimization level, change for debugging
LEVEL = 3

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

SRCS = main.cpp src/amplitude.cpp src/color.cpp src/eikonal.cpp src/tinyxml2.cpp \
	/home/tom/Uni/phd/NnLO/phase_space_parameterization/src/Tree.cpp /home/tom/Uni/phd/NnLO/phase_space_parameterization/src/PhaseSpace.cpp \
	/home/tom/Uni/phd/NnLO/phase_space_parameterization/src/Utilities.cpp
OBJS = $(patsubst %.cpp,%.o,$(SRCS))
EXEC = $(patsubst %.cpp,%,$(SRCS))

.PHONY: clean

target: $(OBJS) $(EXEC)

%: %.o $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS) -DGNU_NAME_MANGLING -lpthread

# Compilation rule
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -DGNU_NAME_MANGLING $< -o $@

clean:
	rm -f *.o $(EXEC) $(OBJS)
