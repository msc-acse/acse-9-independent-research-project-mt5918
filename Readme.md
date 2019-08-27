# Independent Research Project

- Author: Michael Trapp
- CID: 01627245
- Title: 

# Architecture for this Repository

This repository contains three separate folders. The first folder "5cf633a9f7f8694a2b43950a" (overleaf key) contains the documentation. The second folder, "setup2d", contains data and original program files to produce FEMDEM simulations. The third folder, "Y2D", contains modified code and program files which produces VTK output files.

# Folder "setup2d"
## System requirements and Installation

Tests performed using the original code produced unrealistic simulations. Therefore, the system requirements are not clear. The following system requirements are given with the best intention of the author, but may not be final.

- Operating System: Linux Ubuntu 14.04 LTS. For newer versions, errors are likely to occur.
- Optional: Access to High Performance Computing HPC system at Imperial College. https://www.imperial.ac.uk/computational-methods/hpc/

- Software:
    - Visualization Toolkit VTK version 5.8.0. Newer versions fail. Version 5.8.0 cannot be downloaded from https://vtk.org/download/ . 
    - Paraview version 5.4.1. Paraview version 5.6.1 was tested as well. Other versions may work as well. https://www.paraview.org/download/
    - The personal pre and post processor GiD version 10.0.9. Newer versions are incompatible, as they do not support access to the B2D problem set. https://www.gidhome.com/download/

Access to the HPC system is optional. It is advantageous primarily due to the module support for VTK version 5.8.0. 

## 


# Folder "Y2D"
## System requirements and Installation

- Operating System: Linux Ubuntu 18.04.2 LTS. The author used Windows subsystem for Linux Ubuntu 18.04 LTS https://docs.microsoft.com/en-us/windows/wsl/install-win10

- Software:
    - OpenGL libaries libgl1-mesa-dev and libglu1-mesa-dev, installationable via "sudo apt-get install -y libgl1-mesa-dev libglu1-mesa-dev".


