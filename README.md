# QUDA File Utilities

### Package for manipulating QUDA output files

This package contains a set of programs and scripts,
particularly useful to manipulate ASCII and HDF5 files.
The repository supports ASCII/HDF5 files containing disconnected
loops, as well as two- and three-point correlation functions.
These ASCII/HDF5 files are expected to be output of
[quda-QKXTM-Multigrid](https://github.com/ETMC-QUDA/quda-QKXTM-Multigrid).

## Package contents

* **extract:**
This folder contains C/C++ programs to extract HDF5 files
into their corresponding ASCII format.

* **convert:**
This folder contains C/C++ programs to convert ASCII files
into the HDF5 format.

* **scripts:**
This folder contains scripts which:
 * run the programs that extract HDF5 files.
 * reshape files already in ASCII format into a different,
   predefined ASCII form to ensure compatibility with the
   extracted HDF5 files.
 * compare two- and three-point functions and disconnected
   quark loops.

## Dependencies
The package requires linking to the HDF5 and FFTW libraries.

## Compilation/Installation
The user needs to modify make.inc to contain the appropriate
directories of the HDF5 and FFTW libraries.

## Author
**Christos Kallidonis**.

E-mail: c.kallidonis@cyi.ac.cy

Web: christoskallidonis.com
