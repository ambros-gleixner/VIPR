cmake_minimum_required(VERSION 3.3.0)

# set the project name and version
project(VIPR VERSION 1.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set function visibility default to hidden
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

# path to find modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# find GMP
find_package(GMP REQUIRED)
include_directories(${GMP_INCLUDE_DIRS})
set(libs ${libs} ${GMP_LIBRARIES})

# option to install viprcomp
option(VIPRCOMP "Use viprcomp" ON)

# add executables
add_executable(viprttn viprttn.cpp)
add_executable(vipr2html vipr2html.cpp)
add_executable(viprchk viprchk.cpp)

target_link_libraries(viprchk ${libs})

if(VIPRCOMP)
	# Only install viprcomp if working SoPlex is found
	find_package(ZLIB)
	if(ZLIB_FOUND)
		find_package(SOPLEX)
      find_package(Boost 1.77)
      if(Boost_FOUND)
         include_directories(${Boost_INCLUDE_DIRS})
      endif()
		if(SOPLEX_FOUND)
			# include ZLIB
			include_directories(${ZLIB_INCLUDE_DIRS})
			set(libs ${libs} ${ZLIB_LIBRARIES})

			# include SoPlex
			include_directories(${SOPLEX_INCLUDE_DIRS})
			set(libs ${libs} ${SOPLEX_LIBRARIES})

			# add viprcomp target and link
			add_executable(viprcomp viprcomp.cpp)
			add_definitions(-DSOPLEX_WITH_GMP)
			target_link_libraries(viprcomp ${libs})
         message(STATUS "Soplex found.")
		else()
			message(STATUS "viprcomp not installed, because SoPlex could not be found.")
		endif()
	else()
		message(STATUS "viprcomp not installed, because ZLIB could not be found (required by SoPlex).")
	endif()
endif()
