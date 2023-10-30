
export CXX="g++"
export CFLAGS=" -I/home/shilding/packages/nuSQuIDS/src/inlude -I/home/shilding/packages/SQuIDS/src/include -Wno-abi -I/home/shilding/packages/anaconda3/envs/deimos/include -I/home/shilding/packages/anaconda3/envs/deimos/include"
export CXXFLAGS=" -std=c++11"
export LDFLAGS=" -L/home/shilding/packages/nuSQuIDS/src/lib -Wl,-rpath -Wl,/home/shilding/packages/nuSQuIDS/src/lib -lnuSQuIDS -lpthread -L/home/shilding/packages/SQuIDS/src/lib -lSQuIDS -L/home/shilding/packages/anaconda3/envs/deimos/lib -lgsl -lgslcblas -lm -L/home/shilding/packages/anaconda3/envs/deimos/lib -lhdf5 -lhdf5_hl"

export LD_LIBRARY_PATH="/home/shilding/packages/nuSQuIDS/src/lib:/home/shilding/packages/SQuIDS/src/lib:/home/shilding/packages/anaconda3/envs/deimos/lib:/home/shilding/packages/anaconda3/envs/deimos/lib:/home/shilding/packages/SQuIDS/src/lib:/home/shilding/packages/nuSQuIDS/src/lib:/home/shilding/packages/anaconda3/envs/deimos/lib:"
