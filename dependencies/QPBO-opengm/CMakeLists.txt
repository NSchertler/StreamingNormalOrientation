cmake_minimum_required(VERSION 2.6)
project(QPBO-opengm)

# Allow specialization of templates in different namespaces
add_definitions(-fpermissive)

set(SOURCES
	QPBO.cpp
	QPBO_extra.cpp
	QPBO_maxflow.cpp
	QPBO_postprocessing.cpp)

add_library(QPBO-opengm ${SOURCES})

