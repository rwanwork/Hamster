HAMSTER
=======


Introduction
------------

HAMSTER (Helpful Abstraction using Minimum Spanning Trees for Expression Relations) is a system for depicting the relationships between experiments in a microarray data set as a series of minimum spanning trees (MSTs). 

This document accompanies the archive, which also includes:
    * source code in C++,
    * other important documentation and license information, and
    * a very small data file for testing.

This archive does not include any binaries. This program is described in the
paper:

     R. Wan, L. Kiseleva, H. Harada, H. Mamitsuka and P. Horton. HAMSTER: visualizing microarray experiments as a set of minimum spanning trees. Source Code for Biology and Medicine, 4(8), 2009.
     
which we refer to as "the paper" throughout this document.  The software has been updated in 2020 for current compilers.


Requirements
------------

| Software        | Minimum version | Tested version | Required? | Web site                              |
| --------------- | :-------------: | :------------: | :-------: | ------------------------------------- |
| g++             | 4.3.2           | 9.3.0          | Yes       | http://gcc.gnu.org/                   |
| CMake           | 2.8.2           | 3.16.3         | Yes       | http://www.cmake.org/                 |
| Boost library   | 1.39.0          | 1.71.0         | Yes       | http://www.boost.org/                 |
| Graphviz        | 2.8-2.4         |                | No        | http://www.graphviz.org               |
| MPI (Open MPI)  | 1.1-2.3         |                | No        | http://www.open-mpi.org/              |
| Doxygen         | 1.8.9.1         | 1.8.17         | No        | http://www.stack.nl/~dimitri/doxygen/ |


Experiments in the paper using Hamster was executed on Linux systems running Debian 6.0 (squeeze) or CentOS 5.4.  Currently, it is being maintained on an Ubuntu 20.04.1 system (i.e., it's been tested recently on such a system).

Both optional and required tools for compiling or using Hamster is listed in the table above.  The column "Minimum version" refers to the software versions used during software development and when running the experiments in the paper.  They do not represent the minimum requirements; it is possible that lower versions can be used.  The column "Tested version" refers to the versions used for the most recent tests on Ubuntu 20.04.1.

The main tools used by HAMSTER are:
  * The Boost Library is required for compilation and must be both installed and compiled (usually, installation of the Boost Library does not require compilation). The compiled components of Boost that are required are:  regex, graph, mpi, program_options, and serialization.
  * In order to create images, the Graphviz system needs to be available to HAMSTER. It is assumed that the executable neato is available in the path. If it is not found, then Graphviz source files are generated, but not their corresponding images.
  * HAMSTER's workload can be distributed if the Message Passing Interface libraries (MPI) are installed for layout-mst. While any version of MPI that conforms to the standard should work, HAMSTER has been tested only with Open MPI.
  * Doxygen is a documentation system to extract comments that have been placed inline in the source code. See the section below entitled "Software_Documentation" for more information.

Under some Linux distributions, the above libraries can be installed using its associated package manager (such as "aptitude" for Debian and Ubuntu). Consult your Linux distribution's documentation for further information.


Installation Tips
-----------------

We refer the user to the web sites of the above programs and libraries for up-
to-date instructions on how to install them. Some brief notes are provided next
which are specific to HAMSTER and are meant to supplement their respective
documentation. Generally, we suggest that MPI be installed before Boost.
  1. MPI should be installed according to its documentation. The Open MPI web site (for example) also includes small examples that should be tested to ensure that the installation was successful.
  2. Boost library - For version 1.39 and later, unarchive Boost and perform the following steps in `/usr/local/boost_1_39_0/`:
    a. `./bootstrap.sh --with-libraries=regex,graph,mpi,program_options,serialization --libdir=/usr/local/boost_1_39_0/lib --includedir=/usr/local/boost_1_39_0/include`
    b. Edit `/usr/local/boost_1_39_0/project-config.jam` and add to the end of the file:

                 using mpi ;
                 
    c. `./bjam --layout=system --libdir=/usr/local/boost_1_39_0/lib --includedir=/usr/local/boost_1_39_0/include install`

    The above instructions assume the installation directory is in /usr/local/ and, as such, system administration privileges are required. Also, for Boost versions prior to 1.39 and more recent (i.e., in 2020), the instructions may differ slightly; please consult your Boost documentation.
  3. Graphviz should be installed according to its documentation. Afterwards, ensure that the program neato can be found in the path.


Compiling
---------

The HAMSTER software is written in C++ and has been compiled using v4.3.2 of g++. The system has been tested on both 32-bit and 64-bit systems, but it does not make use of any features from 64-bit architectures.

CMake (at least version 2.8) is used to compile the software and it is recommended that an "out-of-source" build is performed so as not to clutter the original source directories. We give some brief instructions below on how to do this:
  1. After having installed Boost, set the variable BOOST_ROOT to the location of Boost if it has not already been set:

           export BOOST_ROOT=/usr/local/boost_1_39_0/
           
  2. Expand the HAMSTER archive in a temporary directory [i.e., `~/tmp/`]. As a result, another directory will be created inside [i.e., `~/tmp/hamster-1.3.0-src/`].
  3. Within ~/tmp/hamster-1.3.0-src/, create a build/ subdirectory and then enter it (Actually, build/ can be anywhere since it will be deleted later; this is just an example.). Then run

           cmake ..
           
  where ".." represents the location of the top-level CMakeLists.txt. By default, this will set up a Makefile to install the program into /usr/local/, which would require system administrator access. To use another directory, type this:

           cmake .. -DCMAKE_INSTALL_PREFIX=~/tmp
      replacing the installation prefix with whatever you prefer.
  4. If your machine has MPI installed but you do not want to have it enabled and running HAMSTER as `mpirun -np 1 ...` is not satisfactory, then you can compile HAMSTER without MPI linked in. To do this, after the previous step with cmake, edit the `LayoutMSTConfig.hpp` in the current directory.  Change the value of `HAVE_MPI` from 1 to 0. Then proceed to the next step.
  5. Type make to compile the C++ source code of HAMSTER. If this succeeds, then the executables will be in the build subdirectory as `build-mst/build-mst` and `layout-mst/layout-mst`.
  6. Finally, type make install to install the software. This copies the important files from the archive to the installation prefix specified in the cmake line above (see "Files_and_Directories" for information about the structure) . The `~/tmp/hamster-1.3.0-src/` directory, including the `build/` directory, can now be deleted, unless you are interested in viewing the source code.


Software Documentation
----------------------

HAMSTER was developed with inline comments that can be extracted using the Doxygen documentation system. They can be created using make as well. Here, we show how to generate the documentation for Built-MST.
  1. In the directory where you ran make to compile the executable, type `make build-mst-doc`. If Doxygen can be found, then the documentation will be generated in the `doc/` subdirectory.
  2. You may choose to do a make install which would put this documentation in the `doc/` subdirectory under the installation prefix. See "Files_and Directories" for more information.

The following HTML file can be used if the documentation has been installed:  documentation.html. Of course, the documentation has to be created first according to step 1.


Files and Directories
---------------------

After installation, the following directory structure should result:

    .                            Directory you specified when you ran cmake
    ├── bin                      Binary files
    │   ├── build-mst            build-mst executable
    │   └── layout-mst           layout-mst executable
    ├── data
    │   ├── build-mst.cfg        Sample configuration file for build-mst
    │   ├── layout-mst.cfg       Sample configuration file for layout-mst
    │   ├── sample.attr          Attributes for sample data
    │   └── sample.data          Sample data file
    ├── doc
    │   ├── AUTHORS              Authors of the software
    │   ├── ChangeLog            History of changes (pre-GitHub)
    │   ├── documentation.html   Top-level HTML file for viewing Doxygen-generated documents
    │   └── VERSION              How to determine the software version
    ├── CMakeLists.txt           Top-level CMakeLists.txt
    ├── LICENSE                  Copy of GNU GPL license v3
    └── README.md                This file in markdown

The HAMSTER system is comprised of two executables: `build-mst` and `layout-mst`.


Running HAMSTER
---------------

### Configuration

Running either executable (`build-mst` or `layout-mst`) without any arguments displays a list of available options. Options can also be placed in configuration files called `build-mst.cfg` and `layout-mst.cfg`, respectively. In the configuration file, the options take the form of:

     <key> = value

where <key> is the name of the command-line option. Options on the command-line take priority over those in the configuration files.

To see a list of available options, use the `--help` option flag. Verbose and debugging outputs are available using `--verbose` and `--debug`, respectively for both programs.

### Experiment Attributes

Every experiment is assigned default attributes which include a color and a shape. Since Graphviz is used to layout the MSTs, the names of both are based on the terms used by Graphviz. Colors are given in hexadecimal format (#C0C0C0 represents gray). Consult the documentation for Graphviz for more information.  In particular, take a look at the documents that talk about "Colors" and "Node Shapes". Some example colors and shapes are as follows:

|-------|-------|-------|-------|-------|
|#FFFFD5|#FFFF80|#FFFF2A|#FFFF00|#FFE300|
|#FFC600|#FFAA00|#FF8E00|#FF7100|#FF5500|
|#FF3900|#FF1C00|#FF0000|#B6FF00|#49FF00|
|#00FF24|#00FF92|#00FFFF|#0092FF|#0024FF|
|#4900FF|#B600FF|#FF00DB|#FF006D|

|-------|-------|-------|-------|-------|
|ellipse|house  |box    |diamond|hexagon|


### Sample Run

A sample data file of four experiments and two probes is provided in the `data/` directory to check if the system is working correctly. There are two files in total for this sample:
   1. sample.data
   2. sample.attr

The file sample.data is the main gene expression data file and is required. It is tab-separated with an experiment on each line. The first row and first column are used as headers and may contain any text. Expression levels are floating point values OR the string NULL to indicate a missing value.

The optional file sample.attr lists the attributes of each experiment. The file is a tab-separated file with an experiment on each line and three fields per line. The three fields are: experiment name, color, and shape. Colors and shapes are specified as given in the previous section.

Running HAMSTER effectively means running first build-mst to construct the MSTs and then layout-mst to generate the images. So, with the sample data in the current directory, the following command would apply Euclidean distance and single linkage to our sample data set with verbose output:

`build-mst sample-data --attrs sample.attr --distance euclidean --linkage single --verbose`

The output should look like the following:
==  Distance method:                Euclidean distance
==  Linkage method:                 Single linkage
==  Scoring method:                 Gap-based
==  Microarray filename:            sample.data
==  Attribute filename:             sample.attr
==  Output path:                    ./
==  Microarray dimensions:          4 by 2
==  Experiment attributes read in:  4
==  Number of pairs calculated:     6

A summary of the merging process is provided in the file summary.txt:
4
0 4294967295  4294967295  0 0 0
1 0 1 1 1.41421 69.9261
2 4 2 2.23607 2.82843 100
3 5 3 0 0 0

Except for the first line which gives the number of experiments, the remaining lines indicate the merging process and the score assigned. This file is used by the next program to generate the images:

layout-mst summary.txt

If Open MPI has been installed, then layout-mst has to be run as follows:

mpirun -np <proc> layout-mst

where <proc> is the number of processors to use.

Each MST corresponds to a Graphviz file in the format <number>.graphviz. If Graphviz was found, then a set of PNG images would also have been created. The Graphviz source files can be edited by hand and re-processed manually using neato:

`neato -Tpng <0.graphviz >0.png`


HAMSTER+
--------

Around 2009, HAMSTER was available as a part of a web service at [http://hamster.cbrc.jp/](http://hamster.cbrc.jp/).  

Dubbed HAMSTER+, it was a web server based on HAMSTER which adds a wrapper specific to NCBI Gene Expression Omnibus (GEO) data. Some support for tab-separated data is also available with one significant difference -- a data file in the format of sample.data is accepted by transposed so that each row is a probe and each column is an experiment. This coincides with the format used by the Simple Omnibus Format in Text (SOFT) used by GEO.

Unfortunately, that web site is no longer available.  The web site was not the focus of the above paper (it was only mentioned briefly) and was not the focus of any publication.


About HAMSTER 
-------------

This software was implemented while I was at Kyoto University and the Computational Biology Research Centre (Tokyo, Japan), around 2008-2009).  My contact details:

     E-mail:  rwan.work@gmail.com 

My homepage is [here](http://www.rwanwork.info/).

The latest version of HAMSTER can be downloaded from [GitHub](https://github.com/rwanwork/Hamster).

If you have any information about bugs, suggestions for the documentation or just have some general comments, feel free to contact me via e-mail or GitHub.

As noted above, many tools were used used for the development of HAMSTER. In addition, the text version of this document was generated using `html2text`.


Copyright and License
---------------------

     HAMSTER (Helpful Abstraction using Minimum Spanning Trees for Expression Relations)
     Copyright (C) 2009-2020 by Raymond Wan

HAMSTER is distributed under the terms of the GNU General Public License (GPL, version 3 or later) -- see the file LICENSE for details.

Permission is granted to copy, distribute and/or modify this document under the terms of the GNU Free Documentation License, Version 1.3 or any later version published by the Free Software Foundation; with no Invariant Sections, no Front-Cover Texts and no Back-Cover Texts. A copy of the license is included with the archive as LICENSE.
