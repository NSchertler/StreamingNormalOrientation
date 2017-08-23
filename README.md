# Streaming Normal Orientation
This repository contains the code for the paper "Towards Globally Optimal Normal Orientations for Large Point Clouds"

# Build Instructions
The repository contains Visual Studio 2013 solutions as well as CMake files. As soon as the dependencies are fulfilled, the code should build out of the box. The CMake version allows you to choose the included components. The Visual Studio solution includes a project with minimal requirements and a full project. Make sure to compile as x64.

## Dependencies
### Eigen3
Eigen3 is a mandatory dependency. The Visual Studio solution uses the environment variable `EIGEN3_ROOT` to find the library. The CMake version additionally looks up some common directories.

### Boost
The following Boost components are required: `thread`, `system`, `filesystem`. The environment variable `BOOST_ROOT` should be set. Under Windows, the libraries must be built in `%BOOST_ROOT%/lib64-msvc-14.0` (this is the directory structure used by the pre-built binary installer). The Visual Studio solution uses statically linked libraries, whereas the CMake version uses dynamically linked libraries. In this case, the shared libraries must be in the PATH.

### OpenGM
OpenGM is used for graph-based optimization. The necessary files are included as submodules of the repository. Make sure to check them out.

### HDF5 (optional)
The HDF5 library can be used to serialize optimization models. The environment variable `HDF5_ROOT` should be set.

### MaxFlowLib and VIGRA (optional)
MaxFlowLib and VIGRA are used by OpenGM for the LSA-TR solver. They require the environment variables `MAXFLOW_ROOT` and `VIGRA_ROOT`. MaxFlowLib must be compiled and the according CMake find script might not work under Linux.

### CUDA (optional)
CUDA can be used to enable a simple PCA-based normal estimation step.

### OpenMP (optional)

# Usage
The program is a command line application with the following parameters:

```
Parameters: filename [options]
	filename        first argument must be the input file name (PLY OR BIN file without extension)
options:
	-stream         streaming orientation.
	-prepare        convert the input PLY file to the streaming BIN format.
	-randomize      randomize normal orientations before the optimization.
	-estimateNormals    run a simple PCA-based normal estimation step.
	-calculateEnergies  calculate initial and final energies if file size is less than 50 MB.
	-k %i           number of neighbors.
	-r %f           neighborhood search radius.
	-thetaS %f      theta_S from paper.
	-thetaAcc %f    theta_acc from paper.
	-flipcrit [hoppe|xie]  flip criterion.
```

The input PLY file must be in ASCII format with the following properties where the color and normal are optional (but normal must be present if color is present). The header may also be absent:

```
property float x
property float y
property float z
property float nx
property float ny
property float nz
property uchar red
property uchar green
property uchar blue
```

During execution, some intermediary files are created and statistics are printed to the console. The final result is saved as an ASCII PLY file.