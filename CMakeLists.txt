cmake_minimum_required(VERSION 3.20.0)
project(Mod-InterleaveAVX_MultiThreads)
include(CTest)

set(SOURCES
    CTool.cpp
    main.cpp
    CSimulate.cpp
    CModulate.cpp
    CLDPC.cpp
    CDecoder_OMS.cpp
    CDecoder_FAID.cpp
    CDecoder_OMSBF.cpp
    CDecoder_OMS_DTBF.cpp
    CChannel.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

# set(CMAKE_CXX_FLAGS
#     "${CMAKE_CXX_FLAGS} -march=native -fast -ansi-alias -no-prec-div -qopt-prefetch -unroll-aggressive -auto-ilp32 -fopenmp -finline -funroll-loops -m64 -fma -mkl"
# )

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -march=native -fast -ansi-alias -fopenmp -finline -funroll-loops -m64 -fma -qmkl"
)

# # for debug 
# set(CMAKE_CXX_FLAGS "-O0 -march=native -fopenmp -qmkl")

target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
