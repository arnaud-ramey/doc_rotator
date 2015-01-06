PROJECT( doc_rotator )
cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE Release)

FIND_PACKAGE( OpenCV REQUIRED )
ADD_EXECUTABLE( doc_rotator doc_rotator.cpp )
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TARGET_LINK_LIBRARIES( doc_rotator ${OpenCV_LIBS} )
