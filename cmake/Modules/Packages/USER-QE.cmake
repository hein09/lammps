if(PKG_USER-QE)
  enable_language(Fortran)
  enable_language(C)

  find_package(QE REQUIRED)
  include_directories(${QE_INCLUDE_DIRS})
  list(APPEND LAMMPS_LINK_LIBS ${QE_LIBRARIES} ${QE_EXT_LIBS} ${LAPACK_LIBRARIES})

  file(GLOB USER-QE_SOURCES ${LAMMPS_SOURCE_DIR}/USER-QE/[^.]*.f90)
  list(APPEND LIB_SOURCES ${USER-QE_SOURCES})
endif()
