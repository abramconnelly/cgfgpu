cmake_minimum_required(VERSION 2.8)
set(PROJNAME cgf)
Project(${PROJNAME})
Message(STATUS "-------------------------------")
Message(STATUS "Processing Project ${PROJNAME}:")

####################################################################################
# Bootstrap
#
set( BASE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
find_path ( HELPERS "Helpers.cmake" HINTS 
    .
    ${CMAKE_MODULE_PATH}
    ${BASE_DIRECTORY}/utils
    ${BASE_DIRECTORY}/../utils
    ${BASE_DIRECTORY}/../source/utils
    ${BASE_DIRECTORY}/../../source/utils
)
if ( ${HELPERS} STREQUAL "HELPERS-NOTFOUND" )
    set ( CMAKE_MODULE_PATH "***FULL PATH TO***/utils/" CACHE PATH "Full path to /utils/" )
    message ( FATAL_ERROR "\n
    Please set the CMAKE_MODULE_PATH 
    to the full path of for /utils/ above
    and configure again." )
endif()
get_filename_component ( CMAKE_MODULE_PATH ${HELPERS} REALPATH )
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} CACHE PATH "Full path to /utils/" )
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "Executable path" )
include( ${CMAKE_MODULE_PATH}/Helpers.cmake )     
include( ${CMAKE_MODULE_PATH}/Packages.cmake )   

set ( USE_NVTX CACHE BOOL "NSight markers")
if ( USE_NVTX )
	add_definitions(-DUSE_NVTX)
endif()

####################################################################################
# Add CUDA
#
_add_package_CUDA()

####################################################################################
# Find SDSL
#
find_package(SDSL)

if ( SDSL_FOUND )
    message( STATUS "--> Using package SDSL (ver ${SDSL_VERSION})") 
	add_definitions(-DUSE_SDSL)	
	add_definitions(-DMSVC_COMPILER)	
    include_directories(${SDSL_INCLUDE_DIR})	
	LIST(APPEND LIBRARIES_OPTIMIZED ${SDSL_LIB_RELEASE} )
	LIST(APPEND LIBRARIES_DEBUG ${SDSL_LIB_DEBUG} )
	LIST(APPEND PACKAGE_SOURCE_FILES ${CUDA_TOOLKIT_INCLUDE} )    
else()
   message ( FATAL_ERROR "---> Unable to find package SDSL")
endif()

####################################################################################
# Compile PTX Files
#
get_filename_component ( _KERNEL_SEARCH_PATH "${BASE_DIRECTORY}/cgfx4" REALPATH )
file(GLOB CUDA_FILES "${_KERNEL_SEARCH_PATH}/*.cu" "${_KERNEL_SEARCH_PATH}/*.cuh")
message ( STATUS "Kernels Search Path: ${_KERNEL_SEARCH_PATH}" )
message ( STATUS "Kernels CUDA Build:  ${CUDA_FILES}" )

_COMPILEPTX ( SOURCES ${CUDA_FILES} TARGET_PATH ${EXECUTABLE_OUTPUT_PATH} GENERATED CUDA_PTX GENPATHS CUDA_PTX_PATHS INCLUDE "${CMAKE_CURRENT_SOURCE_DIR},${GVDB_INCLUDE_DIR}" OPTIONS -arch=compute_30 -code=sm_30 -Xptxas -v,-abi=no,-dlcm=cg --ptxas-options=-v -O3 --use_fast_math )

#####################################################################################
# Asset Path
#
if ( NOT DEFINED ASSET_PATH ) 
  get_filename_component ( _assets "${CMAKE_CURRENT_SOURCE_DIR}/data" REALPATH )
  set ( ASSET_PATH ${_assets} CACHE PATH "Full path to /data/" )  
endif()

###################
# Source files
#
file(GLOB SOURCE_FILES *.cpp *.hpp *.h cgf4/*.cpp cgf4/*.hpp cgf4/*.h cgfx4/*.cpp cgfx4/*.hpp cgfx4/*.h)
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/cgf4")    
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/cgfx4")    


###################################################################################
# Executable
#
_LIST_ALL_SOURCE ( FILES ${SOURCE_FILES} ${COMMON_SOURCE_FILES} ${PACKAGE_SOURCE_FILES} ${UTIL_SOURCE_FILES} ${UTIL_OPTIX_FILES} ${UTIL_OPTIX_KERNELS} ${CUDA_FILES} ${PTX_SOURCES} )
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")    
add_definitions(-DGVDB_IMPORTS -DGLEW_STATIC -DGLEW_NO_GLU)  
add_definitions(-DASSET_PATH="${ASSET_PATH}/")
add_executable (${PROJNAME} ${ALL_SOURCE_FILES} )
set_property ( TARGET ${PROJNAME} APPEND PROPERTY DEPENDS ${PTX_SOURCES} )

#####################################################################################
# Additional Libraries
#
_LINK ( PROJECT ${PROJNAME} OPT ${LIBRARIES_OPTIMIZED} DEBUG ${LIBRARIES_DEBUG} PLATFORM ${PLATFORM_LIBRARIES} )

#####################################################################################
# Windows specific
#
_MSVC_PROPERTIES()

#####################################################################################
# Install Binaries and PTX
#
_DEFAULT_INSTALL_PATH()
_INSTALL_PTX ( FILES ${CUDA_PTX_PATHS} DESTINATION ${EXECUTABLE_OUTPUT_PATH} OUTPUT INSTALL_LIST )

# Install to output location
install ( FILES ${INSTALL_LIST} DESTINATION ${BIN_INSTALL_PATH} )
install ( TARGETS ${PROJNAME} DESTINATION ${BIN_INSTALL_PATH} )
install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_INSTALL_PREFIX})

###########################
# Done
message ( STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}" )
message ( STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}" )
message ( STATUS "EXECUTABLE_OUTPUT_PATH: ${EXECUTABLE_OUTPUT_PATH}" )
message ( STATUS "ASSET_PATH: ${ASSET_PATH}")


