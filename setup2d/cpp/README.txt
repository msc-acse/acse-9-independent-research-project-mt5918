this is directory {localdirectory}/cpp

this directory contains c++ files

ygen.cpp: cpp file
- file to generate the y code
- for this purpose, the finite element meshing program gmsh is used
- the cpp file is automatically called when invoking the {localdirectory}/scripts/testscript.sh with bash
- the cpp output file is a {localdirectory}/gid/{testfile}/{testfile}.y file which is used further (Y2D code is applied on it)