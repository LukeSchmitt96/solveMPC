## Need to install osqp-eigen and some of its dependencies. Find information [here](https://github.com/robotology/osqp-eigen).
The required dependencies are:
- [cmake](https://cmake.org/install/)
- [osqp](http://osqp.readthedocs.io/en/latest/index.html)
- Eigen3 is already included in this repo and since it only includes header files, we don't need to build it.

## Get cmake:
- Download files from [here](https://github.com/Kitware/CMake/releases/download/v3.19.1/cmake-3.19.1.tar.gz).
- Extract the files, enter the directory in the command line and run:
```bash
./bootstrap
make
make install
```


## Get osqp:
```bash
git clone --recursive https://github.com/oxfordcontrol/osqp
cd osqp
mkdir build
cd build
cmake -G "Unix Makefiles" ..
cmake --build .
cmake --build . --target install
```

## Get osqp-eigen:
```bash
cd ~
git clone https://github.com/robotology/osqp-eigen.git
cd osqp-eigen
mkdir build && cd build
cmake ../
make
sudo make install
```

## Get solveMPC
Once all dependencies are installed, clone the repo using:
```bash
git clone https://github.com/LukeSchmitt96/solveMPC.git
```

## Build solveMPC
To build the project for the first time:
```bash
cd [project-directory]
cmake .
cmake --build ./
```

For each subsequent build, just run:
```bash
cmake --build ./
```

To start solver:
```bash
./solveMPC.cpp
```

Note that you may need to run it once or twice before it 'takes'.

Note that you can set a verbose output by adding a true or false flag to the command line statement. It is set to `false` by default.
```bash
./solveMPC true
```
