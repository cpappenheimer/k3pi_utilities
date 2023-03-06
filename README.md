# K3Pi Utilities

## Building the code
```
mkdir build
cmake -S <path to K3PiStudies root> -B <path to build dir> 
```
(add `-DCMAKE_BUILD_TYPE=Debug` flag for debug build)
```
cmake --build <path to build dir>  
```
(`-v` for verbose builds, `-j N` for parallel builds on `N` cores)