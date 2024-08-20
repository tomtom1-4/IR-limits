%.exe : %.o src/amplitude.o src/color.o src/eikonal.o src/tinyxml2.o /home/tom/Uni/phd/NnLO/phase_space_parameterization/src/Tree.o /home/tom/Uni/phd/NnLO/phase_space_parameterization/src/PhaseSpace.o /home/tom/Uni/phd/NnLO/phase_space_parameterization/src/Utilities.o
		g++ -std=c++17 -O3 -o $@ $^ -Wl,-rpath=/home/tom/Documents/software/software/recola_otter/recola2-2.2.3/ /home/tom/Documents/software/software/recola_otter/recola2-2.2.3/librecola.so -lpthread

%.o : %.cpp src/amplitude.cpp src/color.cpp src/eikonal.cpp src/tinyxml2.cpp
		g++ -std=c++17 -O3 -o $@ -c -I/home/tom/Documents/software/software/recola_otter/recola2-2.2.3/include -DGNU_NAME_MANGLING $<

clean:
		rm -f *.o *.exe

.PHONY: all clean
.PRECIOUS: src/amplitude.o src/color.o src/eikonal.o src/tinyxml2.o /home/tom/Uni/phd/NnLO/phase_space_parameterization/src/Tree.o /home/tom/Uni/phd/NnLO/phase_space_parameterization/src/PhaseSpace.o /home/tom/Uni/phd/NnLO/phase_space_parameterization/src/Utilities.o