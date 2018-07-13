
#####################################################################################
# Find SDSL
#
unset(SDSL_FOUND CACHE)
unset(SDSL_INCLUDE_DIR CACHE)

if ( NOT DEFINED SDSL_ROOT_DIR )
  if (WIN32)
    get_filename_component ( BASEDIR "${CMAKE_MODULE_PATH}/shared_sdsl" REALPATH )
  else()
    get_filename_component ( BASEDIR "${CMAKE_MODULE_PATH}/shared_sdsl" REALPATH )
  endif()
  set ( SDSL_ROOT_DIR ${BASEDIR} CACHE PATH "Location of SDSL library" FORCE)
endif()
message ( STATUS "Searching for SDSL at.. ${SDSL_ROOT_DIR}")
set( SDSL_FOUND "YES" )

if ( SDSL_ROOT_DIR )

    #-- Paths to SDSL Library (cached so user can modify)
	set ( SDSL_INCLUDE_DIR "${SDSL_ROOT_DIR}/include" CACHE PATH "Path to include files" FORCE)
	set ( SDSL_LIB_DIR "${SDSL_ROOT_DIR}/lib" CACHE PATH "Path to libraries" FORCE)		

	#-------- Locate Header files
    set ( OK_H "0" )
	_FIND_FILE ( SDSL_HEADERS SDSL_INCLUDE_DIR "sdsl/sdsl_concepts.hpp" "sdsl/sdsl_concepts.hpp" OK_H )	
	if ( OK_H EQUAL 1 ) 
	    message ( STATUS "  Found. SDSL Header files. ${SDSL_INCLUDE_DIR}" )
	else()
	    message ( "  NOT FOUND. SDSL Header files" )
		set ( CUDPP_FOUND "NO" )
	endif ()

    #-------- Locate Library	
     set ( OK_LIB "0" )		
  	_FIND_FILE ( LIST_LIB SDSL_LIB_DIR "sdsl_${MSVC_YEAR}x64.lib" "libsdsl.a" OK_LIB )	        
	_FIND_FILE ( LIST_LIB SDSL_LIB_DIR "sdsl_${MSVC_YEAR}x64d.lib" "libsdsl.a" OK_LIB )	        	
	message ( STATUS "  Searching for SDSL LIBRARY: sdsl_${MSVC_YEAR}x64.lib or libsdsl.a" )
	if ( OK_LIB EQUAL 2 ) 
	   message ( STATUS "  Found. SDSL Library. ${SDSL_LIB_DIR}" )	   
	else()
	   message ( "SDSL Library not found. Please specify path to /shared_sdsl with SDSL library" )
	   message ( "Libraries found: ${LIST_LIB}" )	   
	endif()

endif()
 
if ( ${SDSL_FOUND} STREQUAL "NO" )
   message( FATAL_ERROR "
      Please set SDSL_ROOT_DIR to the root location 
      of the /shared_sdsl library containing /include and /lib.
      Not found at SDSL_ROOT_DIR: ${SDSL_ROOT_DIR}\n"
   )
endif()

set ( SDSL_LIB_DIR ${SDSL_LIB_DIR} CACHE INTERNAL "" FORCE)

if ( WIN32 ) 
	set ( SDSL_LIB_RELEASE "${SDSL_LIB_DIR}/sdsl_${MSVC_YEAR}x64.lib" CACHE INTERNAL "" FORCE)
	set ( SDSL_LIB_DEBUG "${SDSL_LIB_DIR}/sdsl_${MSVC_YEAR}x64d.lib" CACHE INTERNAL "" FORCE)
else()
	set ( SDSL_LIB_RELEASE "${SDSL_LIB_DIR}/libsdsl.a" CACHE INTERNAL "" FORCE)
	set ( SDSL_LIB_DEBUG "${SDSL_LIB_DIR}/libsdsl.a" CACHE INTERNAL "" FORCE)
endif()

#-- We do not want user to modified these vars, but helpful to show them
message ( STATUS "  SDSL_ROOT_DIR:    ${SDSL_ROOT_DIR}" )
message ( STATUS "  SDSL_INCLUDE_DIR: ${SDSL_INCLUDE_DIR}" )
message ( STATUS "  SDSL_LIB_DIR:     ${SDSL_LIB_DIR}" )

mark_as_advanced(SDSL_FOUND)






