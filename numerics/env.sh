export PHASE_SPACE_PATH=/home/tom/Uni/phd/NnLO/phase_space_parameterization
export STRIPPER_PATH=/home/tom/Documents/software/software/Stripper/Stripper
export RECOLA2_PATH=/home/tom/Documents/software/software/recola_otter_nf5/recola2-2.2.3
export OPENLOOPS_PATH=/home/tom/Documents/software/software/OpenLoops

export LD_LIBRARY_PATH="${STRIPPER_PATH}/libs":/usr/local/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH="${STRIPPER_PATH}/libs":$LIBRARY_PATH
export LD_LIBRARY_PATH=$RECOLA2_PATH:$LD_LIBRARY_PATH
export LIBRARY_PATH=$RECOLA2_PATH:$LIBRARY_PATH
export LD_LIBRARY_PATH=/home/tom/local/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH
export LIBRARY_PATH=/home/tom/local/lib/x86_64-linux-gnu/:$LIBRARY_PATH
export PKG_CONFIG_PATH=/home/tom/local/lib/x86_64-linux-gnu/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH="${OPENLOOPS_PATH}/lib":$LD_LIBRARY_PATH
export LIBRARY_PATH="${OPENLOOPS_PATH}/lib":$LIBRARY_PATH
export LD_LIBRARY_PATH=$PHASE_SPACE_PATH:$LD_LIBRARY_PATH