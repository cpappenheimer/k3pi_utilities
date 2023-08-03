cmake_minimum_required(VERSION 3.24)

### project
project(K3PiStudiesUtils VERSION 1.0
        LANGUAGES CXX)

### find packages
find_package(ROOT 6.16 CONFIG REQUIRED)
find_package(Boost 1.50 REQUIRED)

### get path to include dirs
cmake_path(GET CMAKE_CURRENT_SOURCE_DIR PARENT_PATH K3PISTUDIESUTILS_ROOT_DIR)
set(K3PISTUDIESUTILS_INC_DIR "${K3PISTUDIESUTILS_ROOT_DIR}/include")

### add library
add_library(K3PiStudiesUtils SHARED K3PiStudiesUtils.cpp "${K3PISTUDIESUTILS_INC_DIR}")
set_target_properties(K3PiStudiesUtils PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
### include dirs
target_include_directories(K3PiStudiesUtils 
                            PUBLIC "${K3PISTUDIESUTILS_INC_DIR}")
### add required ext libs
target_link_libraries(K3PiStudiesUtils PUBLIC 
                        ROOT::Core 
                        ROOT::MathCore
                        ROOT::Physics
                        ROOT::ROOTVecOps)