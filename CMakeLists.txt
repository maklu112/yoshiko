cmake_minimum_required( VERSION 3.9 )

#############################################

#GTesting

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                 ${CMAKE_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

#############################################

project( Yoshiko )

OPTION(BUILD_STATIC "Link libraries static" OFF)

set (Yoshiko_Version "2.2.0")

set( Yoshiko_Hdr

  src/ClusterEditingInstance.h
  src/ClusterEditingReduction.h
  src/ClusterEditingSolutions.h
  src/ClusterReductionInstance.h
  src/config.h
  src/config.h.in
  src/CoreAlgorithm.h
  src/CplexInformer.h
  src/Delta.h
  src/DynProgTable.h
  src/Globals.h
  src/Helpers.h
  src/ILPSolver.h
  src/InducedCostsHeuristic.h
  src/KClustifier.h
  src/ParameterizedInstance.h
  src/SolutionChecker.h
  src/WorkingCopyInstance.h

  src/input/ClusterEditingInput.h
  
  src/measurement/SilhouetteValue.h

  src/reductionRules/AlmostCliqueRule.h
  src/reductionRules/CliqueRule.h
  src/reductionRules/CriticalCliqueRule.h
  src/reductionRules/FPTReductionRule.h
  src/reductionRules/HeavyEdgeRule3in1.h
  src/reductionRules/MergingRule.h
  src/reductionRules/ParameterDependentReductionRule.h
  src/reductionRules/SimilarNeighborhoodRule.h

)

set( Yoshiko_Src
  src/ClusterEditingInstance.cpp
  src/ClusterEditingReduction.cpp
  src/ClusterEditingSolutions.cpp
  src/ClusterReductionInstance.cpp
  src/CoreAlgorithm.cpp
  src/Delta.cpp
  src/DynProgTable.cpp
  src/Globals.cpp
  src/Helpers.cpp
  src/ILPSolver.cpp
  src/InducedCostsHeuristic.cpp
  src/KClustifier.cpp
  src/ParameterizedInstance.cpp
  src/SolutionChecker.cpp
  src/WorkingCopyInstance.cpp

  src/input/ClusterEditingInput.cpp
  
  src/measurement/SilhouetteValue.cpp

  src/reductionRules/AlmostCliqueRule.cpp
  src/reductionRules/CliqueRule.cpp
  src/reductionRules/CriticalCliqueRule.cpp
  src/reductionRules/FPTReductionRule.cpp
  src/reductionRules/HeavyEdgeRule3in1.cpp
  src/reductionRules/MergingRule.cpp
  src/reductionRules/ParameterDependentReductionRule.cpp
  src/reductionRules/SimilarNeighborhoodRule.cpp

)


set (YoshikoApp_Src

  src/main.cpp

  src/input/JENAInput.cpp
  src/input/RowSimInput.cpp
  src/input/SIFInput.cpp

  src/output/ClusterEditingOutput.cpp
  src/output/CSVOutput.cpp
  src/output/GMLOutput.cpp
  src/output/MCLOutput.cpp
  src/output/XGMMLOutput.cpp
  src/output/PajekOutput.cpp
  src/output/TableOutput.cpp
  src/output/TableAppOutput.cpp
  src/output/TransClustOutput.cpp

  ${Yoshiko_Src}
)

set (YoshikoApp_Hdr

  src/input/JENAInput.h
  src/input/RowSimInput.h
  src/input/SIFInput.h
  src/input/StreamInput.h

  src/output/ClusterEditingOutput.h
  src/output/CSVOutput.h
  src/output/GMLOutput.h
  src/output/MCLOutput.h
  src/output/PajekOutput.h
  src/output/TableOutput.h
  src/output/TableAppOutput.h
  src/output/TransClustOutput.h
  src/output/XGMMLOutput.h

  ${Yoshiko_Hdr}
)

set (YoshikoLib_Src
	src/input/LibraryInput.cpp
	src/LibraryInterface.cpp
	src/LibraryInterface_wrap.cxx
	${Yoshiko_Src}
)

set (YoshikoLib_Hdr
	src/input/LibraryInput.h
	src/LibraryInterface.h
	${Yoshiko_Hdr}
)

set(YoshikoTest_Src
    src/test/runAllTests.cpp
    src/test/verifyExamples.cpp
    ${Yoshiko_Src}
    ${Yoshiko_Hdr}
    src/input/SIFInput.h
    src/input/JENAInput.h
    src/input/JENAInput.cpp
    src/input/SIFInput.cpp
)

set( GitCommand1 "symbolic-ref HEAD 2> /dev/null | cut -b 12-" )
execute_process( COMMAND git symbolic-ref HEAD  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                 COMMAND tail -c -7 COMMAND tr -d "\\n" OUTPUT_VARIABLE GitBranch )
set( GitCommand2 "log --pretty=format:\\\"%h\\\" -1" )
execute_process( COMMAND git log --pretty=format:%h WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                 COMMAND head -n 1  COMMAND tr -d "\\n" OUTPUT_VARIABLE GitHashtag )

set( YOSHIKO_VERSION \"${GitBranch}-${GitHashtag}\" )

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

find_package(JNI)

if (!JNI_FOUND)
    message("Could not find JNI, you will not be able to build the library!")
else()
    message("JNI found @: ${JNI_INCLUDE_DIRS}")
endif()
    
if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
	set( CPLEX_AUTODETECTION_PREFIXES "$ENV{HOME}/ibm/ILOG;/Applications;$ENV{HOME}/Applications/IBM/ILOG; $ENV{HOME}/ILOG;/opt/ibm/ILOG;$ENV{HOME}/Applications/IBM/ILOG" CACHE STRING "CPLEX automatic lookup directories" )
	foreach( CplexInstallPrefix ${CPLEX_AUTODETECTION_PREFIXES} )
		file( GLOB_RECURSE CplexLibraryPaths ${CplexInstallPrefix}/*/cplex/lib/*/libcplex.a )
		if( CplexLibraryPaths )
			list( LENGTH CplexLibraryPaths NumPath )
			if( NumPath GREATER 1 )
				message( WARNING "Found multiple CPLEX installations in ${CplexInstallPrefix}, automatically selecting the first one, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually if this is not the right selection" )
			else()
				message( "Found CPLEX's installation in ${CPlexInstallPrefix}" )
			endif()
			list( GET CplexLibraryPaths 0 CplexLibraryPath )

			set( MyRegex "${CplexInstallPrefix}/(.*)/cplex/?/lib/(.*)/libcplex.a" )
			string( REGEX MATCH "${MyRegex}" ValidatedCplexLibPath ${CplexLibraryPath} )
			if( NOT ValidatedCplexLibPath )
				message( WARNING "Unable to recognize ${CplexLibraryPath} as a valid CPLEX path" )
			else()
				string( REGEX REPLACE "${MyRegex}" "\\1" CplexHome ${ValidatedCplexLibPath} )
				string( REGEX REPLACE "${MyRegex}" "\\2" CplexArch ${ValidatedCplexLibPath} )

				set( CPLEX_INC_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/include"            CACHE PATH "CPLEX's include dir" )
				set( CPLEX_LIB_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/lib/${CplexArch}"   CACHE PATH "CPLEX's library dir" )
				set( CONCERT_INC_DIR "${CplexInstallPrefix}/${CplexHome}/concert/include"          CACHE PATH "Concert's include dir" )
        set( CONCERT_LIB_DIR "${CplexInstallPrefix}/${CplexHome}/concert/lib/${CplexArch}" CACHE PATH "Concert's library dir" )

				break()
			endif()
		else()
			message( "Unable to find CPLEX and Concert installation in ${CplexInstallPrefix}" )
		endif()
	endforeach()

	if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
		message( "Unable to autodeted CPLEX's and Concert's installation directories, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually" )
	else()
		message( "Autodetected CPLEX's include dir   : ${CPLEX_INC_DIR}" )
		message( "Autodetected CPLEX's library dir   : ${CPLEX_LIB_DIR}" )
		message( "Autodetected Concert's include dir : ${CONCERT_INC_DIR}" )
		message( "Autodetected Concert's library dir : ${CONCERT_LIB_DIR}" )
	endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wno-long-long -Wno-unused-parameter  -Wno-unknown-pragmas -g -ggdb -O0 -DDEBUG -DVERBOSE")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -Wno-long-long -Wno-unused-parameter  -Wno-unknown-pragmas -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_PROFILE "-Wall -Wextra -Wno-long-long -Wno-unused-parameter  -Wno-unknown-pragmas -O3 -g -ggdb -DNDEBUG")
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wno-ignored-attributes -fno-strict-aliasing -std=c++0x" )
add_definitions( -DIL_STD ) #???
configure_file( src/config.h.in ${PROJECT_SOURCE_DIR}/src/config.h )


#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -z noexecstack")

include_directories( ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} ${JNI_INCLUDE_DIRS} "${LIBLEMON_ROOT}/include" "src" "src/externlib" )
link_directories( ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR}  "${LIBLEMON_ROOT}/lib" )

set( CommonLibs
	m
	emon
	ilocplex
	cplex
	concert
	pthread
	dl
)

#############################

#GTest Target

enable_testing()
add_executable(unitTests ${YoshikoTest_Src} )
target_link_libraries(unitTests gtest_main ${CommonLibs})
add_test(unitTests unitTests)

#############################

if (JNI_FOUND)
    if (BUILD_STATIC)
            add_library(library STATIC ${YoshikoLib_Src} ${YoshikoLib_Hdr} )
            set_target_properties(library PROPERTIES OUTPUT_NAME "yoshiko.${Yoshiko_Version}")
            target_link_libraries( library ${CommonLibs} "-static")
    else()
            add_library(library SHARED ${YoshikoLib_Src} ${YoshikoLib_Hdr} )
            set_target_properties(library PROPERTIES OUTPUT_NAME "yoshiko.${Yoshiko_Version}")
            target_link_libraries( library ${CommonLibs})
    endif()
endif()

add_executable( executable ${YoshikoApp_Src} ${YoshikoApp_Hdr})
set_target_properties(executable PROPERTIES OUTPUT_NAME "yoshiko.${Yoshiko_Version}")

if (BUILD_STATIC)
	target_link_libraries( executable ${CommonLibs} "-static")
else()
	target_link_libraries( executable ${CommonLibs})
endif()


