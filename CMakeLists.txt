PROJECT( doc_rotator )
cmake_minimum_required(VERSION 2.8)
# Debug, Release, RelWithDebInfo and MinSizeRe
set(CMAKE_BUILD_TYPE RelWithDebInfo)
SET(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra") # add extra warnings

FIND_PACKAGE( OpenCV REQUIRED )
ADD_EXECUTABLE( doc_rotator doc_rotator.cpp )
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TARGET_LINK_LIBRARIES( doc_rotator ${OpenCV_LIBS} )
