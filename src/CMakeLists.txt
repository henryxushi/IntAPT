cmake_minimum_required(VERSION 2.6)
project(IntAPT)

set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(INTAPT_DIR /home/cbil/Henry/IntAPT_v1.1_code)
set(BOOST_LIB_DIR /usr/lib/x86_64-linux-gnu)
set(BOOST_INCLUDE_DIR /usr/include)
set(BAMTOOLS_LIB_DIR /home/cbil/Henry/bamtools-2.4.1/lib)
set(BAMTOOLS_INCLUDE_DIR /home/cbil/Henry/bamtools-2.4.1/include)
set(EIGEN_INCLUDE_DIR /home/cbil/Henry/eigen-eigen-5a0156e40feb)

include_directories(. ${BOOST_INCLUDE_DIR} ${BAMTOOLS_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR} ${INTAPT_DIR})
link_directories(${BOOST_LIB_DIR} ${BAMTOOLS_LIB_DIR})

add_executable(${PROJECT_NAME} IntAPT.cpp process_junc.cpp mergerange.cpp filter_bam.cpp merge_instances.cpp readinstance.cpp utility.cpp options.cpp Info.cpp bedio.cpp rangeset.cpp)

target_link_libraries(${PROJECT_NAME} pthread boost_math_c99 bamtools boost_program_options boost_regex boost_graph boost_thread boost_system boost_filesystem)
