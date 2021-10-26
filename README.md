# EDIpack: A *parallel* Exact Diagonalization solver for Quantum Impurity problems

A Lanczos based solver for generic quantum impurity models exploiting distributed memory MPI parallelisation. This software focuses on the *normal* case (as opposed to superconducting or spin non-conserving cases) including long range magnetic ordering and arbitrary unit cells. 
See [arxiv.org/2105.06806](https://arxiv.org/abs/2105.06806). 

### Dependencies

The code is written around the SciFortran library. Dependencies are:   

* gfortran > 9.0.0 **OR** ifort  > 13.0
* cmake > 3.0.0    
* [MPI](https://github.com/open-mpi/ompi)
* [SciFortran](https://github.com/QcmPlab/SciFortran)


### Installation

Installation is  available using CMake.    

Clone the repo:

`git clone https://github.com/QcmPlab/EDIpack EDIpack`

and from the just created directory make a standard out-of-source CMake compilation:

`mkdir build`  
 `cd build`  
`cmake ..`     
`make`     
`make install`   

The `CMake` compilation can be controlled using the following additional variables, default values between `< >`:   

* `-DPREFIX=prefix directory <~/opt/edipack/PLAT/[VERSION]>` 
* `-DUSE_MPI=<yes>/no`  
* `-DVERBOSE=yes/<no> `  
* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  

### System loading: 

The library can be loaded into the operative system using one of the following, automatically generated, methods:    

* environment module file `~/.modules.d/edipack/<PLAT>`  
* homebrew `bash` script `<PREFIX>/bin/edipack_configvars.sh`
* pkg-config file in `~/.pkg-config.d/edipack.pc`

### Python binding

Python binding (API) through module `edipy` can be  installed, once the library is successfully loaded in the OS, using the conventional toolchain:

`export F90=mpif90` (required if library has been compiled and installed with MPI support)  

1. `python setup.py install`
2. `pip install .`

Method 2. has the advantage of making `uninstall` operation feasible. 

### Uninstall

The library is removed with the command:

`make uninstall`

from the same building directory as for the installation part. 



For any information contact the author as:  
adriano DOT amaricci @ gmail DOT com

--

***COPYRIGHT & LICENSING***  
Copyright  (c), Adriano Amaricci, Lorenzo Crippa, Alberto Scazzola, Francesco Petocchi, Giacomo Mazza, Luca de Medici, Massimo Capone.  
All rights reserved. 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. 
Using this code please consider to cite the corresponding article [arxiv.org/2105.06806](https://arxiv.org/abs/2105.06806)

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

--



