if(PKG_USER-QE)
  enable_language(Fortran)
  enable_language(C)

  find_package(QE REQUIRED)
  include_directories(${QE_INCLUDE_DIRS})

  add_library(qepw ${LAMMPS_SOURCE_DIR}/USER-QE/fix_pw.f90)
  target_include_directories(qepw PRIVATE ${QE_PW_DIR} ${QE_INCLUDE_DIRS})
  target_link_libraries(qepw PRIVATE ${QE_PW_LIBRARY} ${QE_LIBRARIES})

  #add_library(qecp ${LAMMPS_SOURCE_DIR}/USER-QE/fix_cp.f90)
  #target_include_directories(qepw PRIVATE ${QE_CP_DIR} ${QE_INCLUDE_DIRS})
  #target_link_libraries(qepw PRIVATE ${QE_CP_LIBRARY} ${QE_LIBRARIES})

  list(APPEND LAMMPS_LINK_LIBS qepw ${QE_EXT_LIBS} ${LAPACK_LIBRARIES})
endif()
