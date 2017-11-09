cmake_minimum_required(VERSION 3.8)
project(PasiveLocation)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp src/coordinates_system.h src/Pasive_Location.cpp src/Pasive_Location.h)
add_executable(PasiveLocation ${SOURCE_FILES}  src/coordinates_system.h)