# Make input file to compile and link the ramrun main program
# using the ramlib code used for SST.
# This is stand-alone, for testing; it does not need any other
# components of SST.
# The executable is functionally identical to the one

# Edit these for your local environment
#CXX = g++
#FC = gfortran
#FLIBS = -L/usr/local/gfortran/lib -lgfortran
#####FLIBS = -L/usr/local/lib -lgfortran
#CXXFLAGS = -g
#FFLAGS = -g

#ramrun_OBJECTS = RamMain.o RamRunner.o ramlib.o

#runrun: $(ramrun_OBJECTS)
#$(CXX) $(CXXFLAGS) -o ramrun $(OPTIONS) $(ramrun_OBJECTS) $(FLIBS)

#clean:
#-rm -f $(ramrun_OBJECTS) ramrun ramlib.mod


cmake_minimum_required (VERSION 3.8)
if(POLICY CMP0065)
        cmake_policy(SET CMP0065 NEW)
endif()

# set (FC gfortran CACHE STRING "Fortran Compiler short name")
# set (FC ifort CACHE STRING "Fortran Compiler short name")
find_program (CMAKE_Fortran_COMPILER "${FC}" DOC "Fortran Compiler")

# ----------------------------------------------------------------
# This defines what we are doing:
project (RAM
        VERSION 1.0
        LANGUAGES CXX C Fortran
        )

add_library (
        SSTram
        ramlib.f
        ramlib.hh
        RamRunner.cc
        RamRunner.hh
        RamTypes.hh
)

add_executable (ramrun RamMain.cc)
target_link_libraries (
        ramrun
        SSTram
        ${BOOST_FILESYSTEM_LIB}
        ${FLIBS}
)

install (
        TARGETS ramrun
        DESTINATION bin
)
