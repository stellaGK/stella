# stella compilation: Load modules

## Download stella

First, download stella from the GitHub repository:
>> git clone https://github.com/stellaGK/stella.git
>> git submodule update --init --recursive

## Load the system dependent modules

Each system has its own modules installed, therefore a load_modules.SYSTEM file
has been written for each specific system. Source the load_modules.SYSTEM file:

>> source COMPILATION/LoadModules/load_modules.marenostrum
>> source COMPILATION/LoadModules/load_modules.xula

## Compile with CMake

>> cmake . -B COMPILATION/build_cmake
>> cmake --build COMPILATION/build_cmake


-- Found MPI_Fortran: /mnt/lustre/home/spack/spack/opt/spack/linux-rocky8-icelake/oneapi-2023.1.0/intel-oneapi-mpi-2021.9.0-alz6lstmpi2dfpkk3trayc2xqpj54hi6/lib/libmpifort.so (found version "3.1") 




# Set-up system dependent modules

## Module availability

To check the available modules:
>> module avail intel
>> module avail mpi
>> module avail netcdf
>> module avail hdf5

## Module location 

To check the location of a loaded module:
>> module show intel
>> module show mpi
>> module show netcdf
>> module show hdf5

## Local directory without GitHub link

If the stella directory is not linked to GitHub, the git_version library will
most likely refuse to compile, to resolve this issue do:
>> git init


## Check if the correct modules are loaded with CMake

When setting up the modules, it is the easiest to use CMake to check the compilation.
>> cmake . -B COMPILATION/build_cmake




