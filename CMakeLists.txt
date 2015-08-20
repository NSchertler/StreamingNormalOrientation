cmake_minimum_required (VERSION 3.0)
project (Streaming_Normal_Orientation)

option (USE_CUDA "Enable normal estimation via CUDA" OFF)
option (USE_LSATR "Enable LSA-TR optimization" OFF) 
option (USE_HDF5 "Enable model serialization support" OFF)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(LIBRARIES)

set(QPBO_DIR "${CMAKE_CURRENT_LIST_DIR}/dependencies/QPBO-opengm")
add_subdirectory(${QPBO_DIR})
include_directories("${QPBO_DIR}")
set(LIBRARIES ${LIBRARIES} QPBO-opengm)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(Boost REQUIRED COMPONENTS thread system filesystem)
include_directories(${Boost_INCLUDE_DIRS})
set(LIBRARIES ${LIBRARIES} ${Boost_LIBRARIES})

if(USE_HDF5)
	add_definitions(-DWITH_HDF5)
	find_package(HDF5 REQUIRED)
	include_directories(${HDF5_INCLUDE_DIRS})
	set (LIBRARIES ${LIBRARIES} "${HDF5_CXX_LIBRARIES}")
endif(USE_HDF5)

find_package(OpenGM REQUIRED)
include_directories(${OPENGM_INCLUDE_DIR})

if(USE_LSATR)
	add_definitions(-DWITH_LSATR)
	find_package(MaxFlow REQUIRED)
	include_directories(${MAXFLOW_INCLUDE_DIR})
	set (LIBRARIES ${LIBRARIES} "${MAXFLOW_LIBRARIES}")

	find_package(Vigra REQUIRED)
	include_directories(${VIGRA_INCLUDE_DIR})
endif(USE_LSATR)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set (SOURCES  	StreamingNormalOrientation StreamingNormalOrientation/basicProcessing.h
				StreamingNormalOrientation/BlockTimer.h
				StreamingNormalOrientation/common.h
				StreamingNormalOrientation/exceptions.h
				StreamingNormalOrientation/kdTree.h
				StreamingNormalOrientation/normalEstimation.h
				StreamingNormalOrientation/normalOrientation.h
				StreamingNormalOrientation/normalOrientationHooks.h
				StreamingNormalOrientation/OpenPointCloudStreamBinary.h
				StreamingNormalOrientation/OrientationProblem.h
				StreamingNormalOrientation/PointCloudStreamAscii.h
				StreamingNormalOrientation/PointCloudStreamBinary.h
				StreamingNormalOrientation/pointCloudWriter.h
				StreamingNormalOrientation/SignedUnionFind.h
				StreamingNormalOrientation/utils.h
				StreamingNormalOrientation/BlockTimer.cpp
				StreamingNormalOrientation/common.cpp
				StreamingNormalOrientation/kdTree.cpp
				StreamingNormalOrientation/main.cpp
				StreamingNormalOrientation/normalOrientationHooks.cpp
				StreamingNormalOrientation/OrientationProblem.cpp
				StreamingNormalOrientation/utils.cpp )

if(USE_CUDA)
	add_definitions(-DWITH_CUDA)
	
	FIND_PACKAGE(CUDA REQUIRED)

	set(CUDA_SEPARABLE_COMPILATION ON)
	set(LIBRARIES ${LIBRARIES} ${CUDA_LIBRARIES})

	set(SOURCES 	${SOURCES}
					StreamingNormalOrientation/nih/priority_queue.h
					StreamingNormalOrientation/nih/priority_queue_inline.h
					StreamingNormalOrientation/nih/types.h
					StreamingNormalOrientation/SVD/Givens_QR_Factorization_Kernel.cuh
					StreamingNormalOrientation/SVD/Jacobi_Conjugation_Kernel.cuh
					StreamingNormalOrientation/SVD/Kernel_Declarations.cuh
					StreamingNormalOrientation/SVD/Main_Kernel_Body.cuh
					StreamingNormalOrientation/SVD/SVD.cuh
					StreamingNormalOrientation/SVD/SVD.cu
					StreamingNormalOrientation/MortonUtils.cuh
					StreamingNormalOrientation/normalEstimation.cuh
					StreamingNormalOrientation/MortonUtils.cu
					StreamingNormalOrientation/normalEstimation.cu)	
endif(USE_CUDA)

add_definitions(-DGATHER_STATISTICS)
add_definitions(-DIL_STD)
add_definitions(-DWITH_BOOST_GRAPH)
add_definitions(-DWITH_BOOST)
add_definitions(-DEIGEN_DONT_PARALLELIZE)
add_definitions(-std=c++14)
# Allow specialization of templates in different namespaces
add_definitions(-fpermissive)

if(USE_CUDA)
	cuda_add_executable( StreamingNormalOrientation ${SOURCES})
else(USE_CUDA)
	add_executable( StreamingNormalOrientation ${SOURCES})
endif(USE_CUDA)
if(MSVC)
	SET_TARGET_PROPERTIES( StreamingNormalOrientation PROPERTIES LINK_FLAGS "/LIBPATH:${Boost_LIBRARY_DIRS} /LIBPATH:${MAXFLOW_LIBRARY_DIR}" )
endif()

target_link_libraries(StreamingNormalOrientation ${LIBRARIES})
