# A *massively parallel* Exact Diagonalization solver for Quantum Impurity problems

A Lanczos based solver for generic quantum impurity models exploiting distributed memory MPI parallelisation. This software focuses on the *normal* case (as opposed to superconducting or spin non-conserving cases) including long range magnetic ordering and arbitrary unit cells. 

### Dependencies

The code is written around the SciFortran library, which can be found here with installatio notes.   

* SciFortran [https://github.com/aamaricci/SciFortran](https://github.com/aamaricci/SciFortran)  

### Installation

Installation is  available using CMake. In the current v0.0.1 API are only provided in Fortran.   

Clone the repo:

`git clone https://github.com/aamaricci/lib_dmft_ed scifor`

And from the repository directory (`cd lib_dmft_ed`) make a standard out-of-source CMake compilation:

`mkdir build`
`cd build`
`cmake ..`     
`make`     
`make install`   
`make post-install`    

Please follow the instructions on the screen to complete installation on your environment.  
The library can be loaded using one of the following, automatically generated, files :  

* pkg-config file in `~/.pkg-config.d/dmft_ed.pc`  
* environment module file `~/.modules.d/dmft_ed/<PLAT>`  
* homebrew `bash` script `<PREFIX>/bin/configvars.sh`


The `CMake` compilation can be controlled using the following additional variables, default values between `< >`:   

* `-DPREFIX=prefix directory <~/opt/dmft_ed/VERSION/PLAT/[GIT_BRANCH]>` 

* `-DUSE_MPI=<yes>/no`  

* `-DVERBOSE=yes/<no> `  

* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  


For any information contact the author as:  
adriano DOT amaricci @ gmail DOT com

--

***COPYRIGHT & LICENSING***  
Copyright  (c), Adriano Amaricci, Lorenzo Crippa, Alberto Scazzola, Luca de Medici, Massimo Capone.  
All rights reserved. 

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

<!--This program is free software: you can redistribute it and/or modify-->
<!--it under the terms of the GNU General Public License as published by-->
<!--the Free Software Foundation, either version 3 of the License, or-->
<!--(at your option) any later version.-->

<!--You should have received a copy of the GNU General Public License-->
<!--along with this program.  If not, see <http://www.gnu.org/licenses/>.-->

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.


--



