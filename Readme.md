# Independent Research Project

- Author: Michael Trapp
- CID: 01627245
- Title: Elementary Implementation of *VTK* files in C

# Architecture for this Repository

This repository contains three separate folders. The first folder *Documentation* contains the project documentation. The second folder, *setup2d*, contains data and original program files to produce FEMDEM simulations. 

The code for this project is stored in a third folder, *Y*. It contains the original code, and the modified code and program files which produce *VTK* output files.

*Travis* https://travis-ci.com/ is configured in file .travis.yml.

# Folder *Y*

## System requirements and Installation

- Operating System: Linux Ubuntu 18.04.2 LTS. The author used Windows subsystem for Linux Ubuntu 18.04 LTS https://docs.microsoft.com/en-us/Windows/wsl/install-win10 due to lack of alternatives.

- Software:
    - OpenGL libaries libgl1-mesa-dev and libglu1-mesa-dev, installationable via 
    
```bash
    sudo apt-get install -y libgl1-mesa-dev libglu1-mesa-dev
```

## Usage

The folder contains two subfolders, *originalcode* and *finalcode*. *originalcode* contains the initially provided code files. Excessive debugging is neccessary to even get the code running. This debugging is achieved in the other, completely separate subfolder finalcode. The subfolder also contains the algorithm to produce *VTK* output files.

As both subfolders have very similar structure, only the structure of the subfolder *finalcode* is explained. 

*finalcode* contains two folders, *test* and *Y2D*.

*Y2D*, contains the code files of the *Y* *FEMDEM* code. *test* contains a bash script script.sh to produce the output files. The script is run via the command:

```bash
bash script.sh
```

The script remakes the executables and runs the *Yf* program. The subfolder *results* contains valid vtu output files. The files can be opened as a time series using paraview.

The primary code file for this project is file *Yod.c*, which contains the code to project the *VTK* output files. 

# Folder *setup2d*

## System requirements and Installation

***
Tests performed using the original code produced unrealistic simulations. Therefore, the system requirements are not clear. The following system requirements are given with the best intention of the author, but are likely not final.
***

- Operating System: Linux Ubuntu 14.04 LTS. The author used a Windows subsystem for Linux Ubuntu 18.04 LTS https://docs.microsoft.com/en-us/Windows/wsl/install-win10 due to lack of alternatives.

- Software:
    - Visualization Toolkit *VTK* version 5.8.0. Newer versions fail. Version 5.8.0 cannot be downloaded from https://VTK.org/download/ .
    - Paraview version 5.4.1. Paraview version 5.6.1 was tested as well. Other versions may work as well. https://www.paraview.org/download/
    - The personal pre and post processor *GID* version 10.0.9 https://www.GIDhome.com/download/ . Newer versions are incompatible, as they do not support access to the *B2D* problem set.

Tests were perfomed on a local Windows 10 operating system using Windows subsystem for Linux Ubuntu 18.04 LTS. Resorting to remote Linux Ubuntu 14.04 LTS university server ese-pollux.ese.ic.ac.uk introduced extreme performance issues. For testing and small calculations, access to a remote server like the High Performance Computing *HPC* system at Imperial College https://www.imperial.ac.uk/computational-methods/HPC/ is possible, but inconvenient. 
The only reason to consider these remote servers is due to the module support for *VTK* version 5.8.0. Small test codes should ideally be able to be run locally instead if the dependencies were available to maximize code development.

## Usage 

1. Generate a y input file using *GID*. 

Starting *GID* invokes a password window. Entering a false address or ignoring this window appears to have no effect. The specifics and purpose of this window are not entirely clear to the author. Inquiries pertaining to this password should thus be directed towards the *AMCG* group at Imperial College.

For now, no other alternatives than _GID_ are available to create a *y* input file. Generating this file requires the following tasks, details can be found in the author's documentation:

    - Opening GiD. The author only managed to do so on a Windows operating system.
    - Setting the problem type to 'B2D'. If the problem type does not exist, either the *GID* version is wrong or there is another problem.
    - Creating geometries via points, lines and surfaces (the case of creating 3D volumes is not considered here).
    - Generating a mesh from the geometry.
    - Specifying problem data. 
    - Specifying the materials.
    - Specifying initial conditions.
    - Generating the _y_ file. This is achieved via clicking _calculate_ in the calculate menu.

To the knowledge of the author, any command line options whatsoever are not available for _GID_. Within the program, a command line exists to specify numerical values such as coordinates or indices. The commandline is linked with the UI, which makes any automation attempts impossible. Saving a _GID_ project creates a folder with several files. Besides the _y_ file, any other _GID_ files from this folder are irrelevant for further consideration here.

***
Warning: _GID_ is extremely unreliable, program code for _GID_ is inaccessible and the _GID_ documentation is extremely poor. Even when using Linux 14.04 operating system, test files generated by _GID_ almost always resulted in errors when the code was compiled. 
***

For the simulations here, *y* files were generated on specific machines where *GID* seems to generate valid *y* files. An attempt was made to replace _GID_ using _Gmesh_ http://gmsh.info/ . Due to several reasons, the code is left unfinished:

    - Problems implementing the Gmesh API for C++, with the poor alternative to write two C++ files and execute Gmesh in between files
    - Receiving the Y2D code, written in C (instead of C++)
    - Complexity of the data - input of arbitrary 2D geometries via a text file
    - The appearent lack of "scientific relevance" (which, in the author's opinion, is only at most partially justified)

Compatibility issues between *GID* and the operating system may have generated an invalid *y* file. This may be one potential reason why the Travis tests are failing in a clean environment.

2. Configure the following items in the user config file *config.txt*:

- The test file name, which is the y input file generated by GID. 
- The bash path of the SSH Login for Imperial College HPC High Performance Computing System
- The HPC project directory, where the files are stored
- HPC qsubdir, with the qsub and qstat commands 
- A flag to run y2d code locally instead of on HPC (this requires VTK 5.8.0)
- The bash path to paraview

After configuration, one or more bash scripts (depending on whether *HPC* is used) are ready to be executed to generate the output files.

3. Run the following bash script(s):

    - test.sh runs the following programs subsequently, using the *GID*-generated input y file as a command argument 
        - Yf
        - m2vtu and
        - m2vtu_crack
    - script.sh is a submittable HPC script in case the HPC system is used.
    - qstat.sh to check whether the HPC job has finished in case the HPC system was used.
    - getresults.sh fetches the results from the HPC server in case the HPC system was used.

All scripts can be run in a Linux shell via the command:

```bash
    bash myscript.sh
```

where myscript.sh is the name of the script. Ytmp appears to be an additional temporar *y* file but must be present to prevent an error. The exact purpose of this file is not clear to the author.

4. The results are available in folder results in a subfolder named as specified in *config.txt* (default name: *test*) . The files can be opened as a time series using paraview.

***
Programs *Yf*, *m2vtu* and *m2vtu_crack* contain an excessive amount of bugs, as detailed in the project documentation. This may be another reason as to why *Travis* is failing and errors are generated. The original code thus does not run in a clean environment, and potentially does not run using the specified Linux version depending on whether the specific *GID*-generated *y* file is compatible.
***


 


