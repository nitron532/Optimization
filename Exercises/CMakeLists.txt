cmake_minimum_required(VERSION 3.16)

set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")

project(Exercises)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
find_package(Ceres REQUIRED)

file(GLOB EXERCISE_SOURCES "*.cpp")

foreach(SOURCE ${EXERCISE_SOURCES})

    get_filename_component(TARGET_NAME ${SOURCE} NAME_WE)

    add_executable(${TARGET_NAME} ${SOURCE})

    target_link_libraries(${TARGET_NAME}
        PRIVATE
            Eigen3::Eigen
            Ceres::ceres
    )
endforeach()