CMAKE_MINIMUM_REQUIRED(VERSION 3.28.3)

ADD_LIBRARY(Rupert STATIC matrix.h matrix.cpp
                          rotate.h rotate.cpp
                          rupert.h rupert.cpp
           )

PROJECT(test01)
ADD_EXECUTABLE(test01 test01.cpp)
TARGET_LINK_LIBRARIES(test01 PUBLIC Rupert)

PROJECT(line)
ADD_EXECUTABLE(line line.cpp)
TARGET_LINK_LIBRARIES(line PUBLIC Rupert)
