# ROCKCT-Imaging
This repository aims to reconstruct and segment 3D digital rock volumes from binary 2D slices. The resulting 3D rock structures can be visualized using [Paraview](https://www.paraview.org/) for visualization.

## Setup
This project requires the installation of the FEM environment [NeoPZ](https://github.com/labmec/neopz).
The class `TPZVTKGeoMesh` is currently required for VTK output.

Additionally, CMake version 3.11.0 or higher is required.

### CMake Configuration
```markdown
cmake_minimum_required(VERSION 3.11.0)

project(MyProject)

find_package(
  NeoPZ REQUIRED
  HINTS
    ${CMAKE_SOURCE_DIR}/../neopz_install/
    ${CMAKE_SOURCE_DIR}/neopz_install/
)
```
## Reproducibility
To reproduce the main experiment:

1. Place the binary 2D CT slice images in the following directory:
   DADOS_TIF_10Layers/
2. Compile and execute the main program located in:
   RAW3dVTK/main.cpp
3. The following VTK files will be generated:
- `RawTest.vtk`: 3D reconstructed digital rock volume.
- `RawTestoutput.vtk`: Segmented volume with a unique label assigned to each connected component.
- `RawTestordered.vtk`: Segmented volume with connected components ordered by size.


## Main Contributions
This work provides:
- A framework for reconstructing 3D rock volumes from 2D CT slices.
- A segmentation pipeline tailored to heterogeneous rock structures.
- A reproducible experimental setup for quantitative evaluation.
