# beam_down
Solar flare electron beam propagation downwards

The code to simulatione downward electron transport in collisional plasma of solar corona;
Includes self-consistently Coulumb collisions, Langmuir wave generation/absorption, and non-uniform plasma effects

#### If not already installed on your system, download and install git for command line from https://git-scm.com/
#### Launch the Git Bash terminal
```bash
git clone https://github.com/edkontar/wkt wkt
cd wkt
git submodule update --init --recursive --remote
```
#### To compile with f95 or Intel Fortran compiler and build a binary file 
```bash
make
```
#### To plot results and/or make animations using IDL use fw.pro 
```bash
IDL> fw
```
