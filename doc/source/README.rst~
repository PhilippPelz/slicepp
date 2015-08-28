.. Multislice Simulation documentation master file, created by
   sphinx-quickstart on Thu Aug 20 11:30:21 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Slicepp's documentation!
=================================================

About
=====
Slicepp is a multislice simulation software.It simulates an incoming electron wave striking and propagating through a user-defined sample, eventually picked up by a detector. This software is accelerated by the NVIDIA CUDA platform and ArrayFire.

Installation
============

**Libraries Required**

=============================================================  ================= ================
Library							       Version		 Comments
=============================================================  ================= ================
`Armadillo <http://arma.sourceforge.net/>`_ 			>= 5.0		 Required
`ArrayFire <http://www.arrayfire.com/>`_ 			>= 3.0		 Required
`Boost <http://www.boost.org/>`_ 				>= 1.5		 Required
`CUDA Toolkit <https://developer.nvidia.com/cuda-downloads>`_ 	>= 7.0		 Required
`FFTW3 <http://www.fftw.org/>`_ 				>= 3.3		 Required
`HDF5 <https://www.hdfgroup.org/HDF5/>`_ 			>= 1.8		 Required
`NLopt <http://ab-initio.mit.edu/wiki/index.php/NLopt>`_ 	>= 2.4		 Required
`OpenBabel <http://openbabel.org/wiki/Main_Page>`_ 		>= 2.3		 Required
`OpenMP <http://openmp.org/wp/>`_ 				>= 4.0		 Required
`Python2 <https://www.python.org/>`_				2.7		 GUI
=============================================================  ================= ================

**Compilation**

*Debug*
::

	% cd /slicepp/build/debug
	% ./remake.sh
	% make -j2

*Release*
::

	% cd /slicepp/build/release
	% ./remake.sh:

Execution
=========

*From terminal*

Run the following commands in terminal: 
::
	% cd /slicepp/build/release/bin
	%./stem3 [PATH_TO_CONFIG_FILE]

*From GUI*

Run the following commands in terminal: 
::
	% cd /slicepp/GUI_new
	% python window.py

Parameters
==========

mode

*Experiment type: CBED = 1, NBED = 2, STEM = 3, TEM = 4, PTYC = 5*

nthreads

**Scan**

**Wave**

**Beam**

**Output**

**Detector**

**Structure**

**Model**

**Potential**


Programming Guide
=================

**Source Code Structure**

Source code for the main project are in the **/libs** folder. Families of hierarchial code are grouped into subdirectories (e.g. wavefunctions, potentials, detectors). In each of these subdirectories, there is a header file, usually called **[DIRECTORYNAME]_interface.hpp**, that serves as a template for other classes in the same directory.To create a new class, simply inherit from the template and override methods where appropriate. For instance, if one wants to create a new way of calculate potential, one can make a subclass of potential_interface called NewPotential. Then, implement pure virtual and other additional methods in the NewPotential class. When creating a new subdirectory of **/libs**, make sure also to create a CMakeList.txt to change the scopes of the files (examples can be found in existing subdirectories) and modify the CMakeLists in **/libs** to include the newly created subdirectory. Project scope classes should be added in the **/libs** folder.

Important definitions (variables types, constants, precision) are in the file called **stemtypes_fftw3.hpp**, also in the **/libs** folder. New definitions should be appended to this file.

**Parallelization with CUDA/ArrayFire**

ArrayFire will be imported automatically when building the project with cmake. 

When creating CUDA source files (files with .cu extensions), make sure it stands alone and does not include any non-native C++ libraries. The NVCC compiler is known to have compilation issues with fftw3, Boost, etc. Additionally, modify the CMakeList.txt in the same directory to explicitly compile .cu files. An example of implementing raw CUDA files is **CUDA2DPotential.cu** in **/libs/potentials**. Using ArrayFire along side CUDA memory allocation (cudaMalloc()) and libraries (cuBLAS, cuFFT, etc.) sometimes produce wrong results or even memory violations. Therefore, implement CUDA kernels only when absolutely needed as most operations have equivalent ArrayFire counterparts.

**Main Executable**

The main **stem3.cpp** is in the **/stem3** folder. The **Bootstrapper** class handles creation of wavefunctions, potentials, etc. Append to the appropriate **Register[TYPENAME]Types()** functions when implementing new features.

**New Config Parameters**

To add new external parameters, modify appropriate sections of the **read_qsc** class in **/libs/Config_IO**.

API
===
.. toctree::
   :maxdepth: 1

   Functions
   

Support
=======

License
=======

