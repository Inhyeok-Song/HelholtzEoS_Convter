SOURCES=main.cpp helmholtzEoS.cpp read_helmholtz_eos.cpp set_bins.cpp write_eos_table.cpp
INCLUDES=helmholtz.h
HDF5DIR=/cluster/software/hdf5-parallel/1.8.21/gcc--8.3.0/openmpi--3.1.4
HDF5INCS=-I$(HDF5DIR)/include
MPIINCS=-I/cluster/software/openmpi/3.1.4/gcc--8.3.0/include
HDF5LIBS=-L$(HDF5DIR)/lib -lhdf5 -lhdf5_hl
OBJECTS=$(SOURCES:.cpp=.o )

CXX=mpicxx
CFLAGS=-g -O3 -std=c++03
EXTRALIBS=-lm

driver: $(OBJECTS) $(INCLUDES)
	$(CXX) $(CFLAGS) -o helm_eos_maker $(OBJECTS) $(HDF5LIBS) $(EXTRALIBS)

$(OBJECTS): %.o: %.cpp $(INCLUDES)
	$(CXX) $(CFLAGS) $(MPIINCS) $(HDF5INCS) -c $< -o $@

clean:
	rm -f *.o helm_eos_maker