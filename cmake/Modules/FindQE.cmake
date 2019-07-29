# - Find quantum-espresso
# Find the native QE headers and libraries.
#
#  QE_INCLUDE_DIRS - where to find libqecouple.h, etc.
#  QE_LIBRARIES    - List of libraries when using quantum-espresso.
#  QE_FOUND        - True if quantum-espresso found.
#

find_path(QE_COUPLE_DIR libqecouple.h
    PATH_SUFFIXES COUPLE/include
    HINTS $ENV{QEROOT}
    )

find_library(QE_COUPLE_LIBRARY NAMES qecouple HINTS ${QE_COUPLE_DIR}/../src)

set(QE_PW_DIR ${QE_COUPLE_DIR}/../../PW/src)
find_library(QE_PW_LIBRARY NAMES pw HINTS ${QE_PW_DIR})

set(QE_MOD_DIR ${QE_COUPLE_DIR}/../../Modules)
find_library(QE_MOD_LIBRARY NAMES qemod HINTS ${QE_MOD_DIR})

set(QE_FFT_DIR ${QE_COUPLE_DIR}/../../FFTXlib)
find_library(QE_FFT_LIBRARY NAMES qefft HINTS ${QE_FFT_DIR})

set(QE_LAX_DIR ${QE_COUPLE_DIR}/../../LAXlib)
find_library(QE_LA_LIBRARY NAMES qela HINTS ${QE_LAX_DIR})

set(QE_CLIB_DIR ${QE_COUPLE_DIR}/../../clib)
find_library(QE_CLIB_LIBRARY NAMES clib.a HINTS ${QE_CLIB_DIR})

set(QE_IOTK_DIR ${QE_COUPLE_DIR}/../../S3DE/iotk/src)
find_library(QE_IOTK_LIBRARY NAMES iotk HINTS ${QE_IOTK_DIR})

set(QE_UTIL_DIR ${QE_COUPLE_DIR}/../../UtilXlib)
find_library(QE_UTIL_LIBRARY NAMES util PATHS ${QE_UTIL_DIR} NO_DEFAULT_PATH)

set(QE_FOX_DIR ${QE_COUPLE_DIR}/../../FoX/finclude)
find_library(QE_FOX_LIB_COMMON NAMES FoX_common PATHS ${QE_FOX_DIR}/../lib)
find_library(QE_FOX_LIB_DOM NAMES FoX_dom PATHS ${QE_FOX_DIR}/../lib)
find_library(QE_FOX_LIB_FSYS NAMES FoX_fsys PATHS ${QE_FOX_DIR}/../lib)
find_library(QE_FOX_LIB_SAX NAMES FoX_sax PATHS ${QE_FOX_DIR}/../lib)
find_library(QE_FOX_LIB_UTILS NAMES FoX_utils PATHS ${QE_FOX_DIR}/../lib)
find_library(QE_FOX_LIB_WCML NAMES FoX_wcml PATHS ${QE_FOX_DIR}/../lib)
find_library(QE_FOX_LIB_WKML NAMES FoX_wkml PATHS ${QE_FOX_DIR}/../lib)
find_library(QE_FOX_LIB_WXML NAMES FoX_wxml PATHS ${QE_FOX_DIR}/../lib)

set(QE_KS_DIR ${QE_COUPLE_DIR}/../../KS_Solvers)
find_library(QE_KS_LIBRARY NAMES ks_solvers PATHS ${QE_KS_DIR})
find_library(QE_KS_LIB_CG NAMES cg PATHS ${QE_KS_DIR}/CG)
find_library(QE_KS_LIB_DAVID NAMES david PATHS ${QE_KS_DIR}/Davidson)

set(QE_D3_DIR ${QE_COUPLE_DIR}/../../dft-d3)
find_library(QE_D3_LIBRARY NAMES dftd3qe PATHS ${QE_D3_DIR})

find_package(MPI REQUIRED)

set(QE_LIBRARIES ${QE_COUPLE_LIBRARY} ${QE_PW_LIBRARY} ${QE_MOD_LIBRARY} ${QE_FFT_LIBRARY} ${QE_IOTK_LIBRARY})
list(APPEND QE_LIBRARIES ${QE_FOX_LIB_DOM} ${QE_FOX_LIB_UTILS} ${QE_FOX_LIB_WCML} ${QE_FOX_LIB_WKML} ${QE_FOX_LIB_WXML} ${QE_FOX_LIB_SAX} ${QE_FOX_LIB_FSYS} ${QE_FOX_LIB_COMMON})
list(APPEND QE_LIBRARIES ${QE_KS_LIB_CG} ${QE_KS_LIB_DAVID} ${QE_KS_LIBRARY})
list(APPEND QE_LIBRARIES ${QE_D3_LIBRARY} ${QE_LA_LIBRARY} ${QE_UTIL_LIBRARY} ${QE_CLIB_LIBRARY} ${MPI_Fortran_LIBRARIES})
set(QE_INCLUDE_DIRS ${QE_COUPLE_DIR} ${QE_MOD_DIR} ${QE_PW_DIR} ${QE_FFT_DIR} ${QE_LAX_DIR} ${QE_CLIB_DIR} ${QE_IOTK_DIR} ${QE_UTIL_DIR} ${QE_MOD_DIR})

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set QE_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(QE DEFAULT_MSG QE_COUPLE_LIBRARY QE_PW_LIBRARY QE_MOD_LIBRARY QE_FFT_LIBRARY QE_LA_LIBRARY QE_CLIB_LIBRARY QE_IOTK_LIBRARY QE_UTIL_LIBRARY QE_COUPLE_DIR)

mark_as_advanced(QE_COUPLE_DIR QE_COUPLE_LIBRARY PW_LIBRARY QE_MOD_LIBRARY QE_FFT_LIBRARY QE_LA_LIBRARY QE_CLIB_LIBRARY QE_IOTK_LIBRARY)
