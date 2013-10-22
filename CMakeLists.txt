cmake_minimum_required(VERSION 2.8)

project(SMASH NONE)

# include variables for Rokko
include(${ROKKO_INSTALL_DIR}/config/UseRokko.cmake)
include(${ROKKO_INSTALL_DIR}/config/add_rokko_test.cmake)

enable_language(Fortran C)

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  # Almost always OpenMP flags are same both for C and for Fortran.
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
endif(OPENMP_FOUND)

# MPI library
find_package(MPI)

#set(CMAKE_EXE_LINKER_FLAGS ${MPI_C_LINK_FLAGS})

include_directories(${ROKKO_INCLUDE_DIRS})
include_directories(${ROKKO_INCLUDE_DIRS}/rokko_c-fortran)
message(STATUS "INNNNNNC = ${ROKKO_INCLUDE_DIRS}/rokko_c-fortran")
message(STATUS "LIBBBBBB = ${ROKKO_FORTRAN_LIBRARIES}")

set(SOURCES basis.F90 geom.F90 guess.F90 input.F90 int1.F90 machine.F90 math.F90 memory.F90 parallel.F90 scf.F90 scflib.F90 int2a.F90 int2ssss.F90 int2psss.F90 int2ppss.F90
int2psps.F90 int2ppps.F90 int2pppp.F90 int2dsss.F90 int2dpss.F90 int2dsps.F90 int2ddss.F90 int2dpps.F90 int2dsds.F90 int2dspp.F90 int2ddps.F90 int2dpds.F90 int2dppp.F90 int2ddds.F90 int2ddpp.F90 int2dpdp.F90 int2dddp.F90 int2dddd.F90 gradrhf.F90 rysquad.F90 int2rys.F90 dft.F90 dftfunc.F90 lebedev.F90 mp2.F90 grad1.F90 grad2.F90)
set(MODULES module.F90 modulefmt.F90 modulerys.F90)

set(SAMPLE main)
foreach(name ${SAMPLE})
  add_executable(${name} ${name}.F90 ${ROKKO_FORTRAN_LIBRARIES} ${SOURCES} ${MODULES})
  target_link_libraries(${name} ${ROKKO_FORTRAN_LIBRARIES} ${LAPACK_LIB})
endforeach(name)

#set(SAMPLE main)
#foreach(name ${SAMPLE})
#  add_executable(${name} ${name}.F90 ${SOURCES} ${MODULES})
#  target_link_libraries(${name} ${LAPACK_LIB})
#endforeach(name)

#install(TARGETS SMASH ARCHIVE DESTINATION lib LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
 
#copy files
set(PATTERN_FILES  "tests/small.inp" "tests/medium.inp" "tests/large.inp" "tests/tiny.inp")
foreach( pattern_file ${PATTERN_FILES} )
  message(STATUS "Display_CURRENT_BINARARY: ${CMAKE_CURRENT_BINARY_DIR}")
  message(STATUS "Display_CURRENT_SOURCE: ${CMAKE_CURRENT_SOURCE_DIR}")
  message(STATUS "Display_pattern: ${CMAKE_CURRENT_SOURCE_DIR}/${pattern_file}")
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${pattern_file} ${CMAKE_CURRENT_BINARY_DIR}/${pattern_file} COPYONLY)
endforeach( pattern_file )
