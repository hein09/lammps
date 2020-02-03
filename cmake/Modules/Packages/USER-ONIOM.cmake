if(PKG_USER-ONIOM)
    # Oniom command with private fixes
    set(USER-ONIOM_SOURCES_DIR ${LAMMPS_SOURCE_DIR}/USER-ONIOM)
    set(USER-ONIOM_HEADERS ${USER-ONIOM_SOURCES_DIR}/fix_collect.h
                           ${USER-ONIOM_SOURCES_DIR}/fix_oniom.h
                           ${USER-ONIOM_SOURCES_DIR}/oniom.h)
    set(USER-ONIOM_SOURCES ${USER-ONIOM_SOURCES_DIR}/fix_collect.cpp
                           ${USER-ONIOM_SOURCES_DIR}/fix_oniom.cpp
                           ${USER-ONIOM_SOURCES_DIR}/oniom.cpp)
    AddStyleHeader(${USER-ONIOM_SOURCES_DIR}/fix_oniom.h FIX)
    AddStyleHeader(${USER-ONIOM_SOURCES_DIR}/oniom.h COMMAND)

    # Interface to QuantumEspresso
    find_package(qe QUIET)

    # PWScf
    option(BUILD_QE_PW "Build \"fix qe/pw\" that embeds PWScf in LAMMPS" ${qe_FOUND})
    if(BUILD_QE_PW)
        enable_language(Fortran)
        find_package(qe REQUIRED)
        find_package(MPI REQUIRED) #reimport to ensure get MPI::MPI_Fortran

        list(APPEND USER-ONIOM_SOURCES ${USER-ONIOM_SOURCES_DIR}/fix_pw.f90)
        list(APPEND USER-ONIOM_HEADERS ${USER-ONIOM_SOURCES_DIR}/fix_pw.h)
        list(APPEND USER-ONIOM_SOURCES ${USER-ONIOM_SOURCES_DIR}/fix_pw.cpp)
        AddStyleHeader(${USER-ONIOM_SOURCES_DIR}/fix_pw.h FIX)
        target_link_libraries(lammps PRIVATE qe::qe_pw)
    endif()

    # CP
    option(BUILD_QE_CP "Build \"fix qe/cp\" that embeds CP in LAMMPS" ${qe_FOUND})
    if(BUILD_QE_CP)
        enable_language(Fortran)
        find_package(qe REQUIRED)
        find_package(MPI REQUIRED) #reimport to ensure get MPI::MPI_Fortran

        list(APPEND USER-ONIOM_Sources ${USER-ONIOM_SOURCES_DIR}/fix_cp.f90)
        list(APPEND USER-ONIOM_HEADERS ${USER-ONIOM_SOURCES_DIR}/fix_cp.h)
        list(APPEND USER-ONIOM_SOURCES ${USER-ONIOM_SOURCES_DIR}/fix_cp.cpp)
        AddStyleHeader(${USER-ONIOM_SOURCES_DIR}/fix_cp.h FIX)
        target_link_libraries(lammps PRIVATE qe::qe_cpv)
    endif()

    # QE base class
    if(BUILD_QE_CP OR BUILD_QE_PW)
        list(APPEND USER-ONIOM_HEADERS ${USER-ONIOM_SOURCES_DIR}/fix_qe.h)
        list(APPEND USER-ONIOM_SOURCES ${USER-ONIOM_SOURCES_DIR}/fix_qe.cpp)
    endif()

    target_sources(lammps PRIVATE ${USER-ONIOM_SOURCES})
    target_include_directories(lammps PRIVATE ${USER-ONIOM_SOURCES_DIR})
endif()
