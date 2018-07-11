cmake_minimum_required( VERSION 3.9 )

project(George)

set(George_src
	src/main.cpp
	src/GraphGenerator.cpp
	src/GraphGenerator.h
)

add_executable(executable ${George_src})
set_target_properties(executable PROPERTIES OUTPUT_NAME "george")


set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )
include_directories( "${LIBLEMON_ROOT}/include")

target_link_libraries(executable m emon)
