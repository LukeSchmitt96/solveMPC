## Need to install osqp-eigen and some of its dependencies. Find information [here](https://github.com/robotology/osqp-eigen).
The required dependencies are:
- [cmake](https://cmake.org/install/)
- [osqp](http://osqp.readthedocs.io/en/latest/index.html)
- [Eigen3](http://eigen.tuxfamily.org/index.php?title=Main_Page)

## Get cmake:
```bash
sudo apt-get install cmake -y
```

## Get Eigen:
```bash
sudo apt-get install libeigen3-dev -y
```

## Get osqp:
```bash
cd ~
git clone --recursive https://github.com/oxfordcontrol/osqp
cd osqp
mkdir build
cd build
cmake -G "Unix Makefiles" ..
cmake --build .
sudo cmake --build . --target install
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
cd ~
git clone https://github.com/LukeSchmitt96/solveMPC.git
```

## Building the solver
To build the project for the first time:
```bash
cd solveMPC
cmake .
cmake --build ./
```

For each subsequent build, run the following command in the project directory:
```bash
cmake --build ./
```

To start solver, run the following command in the project directory:
```bash
./solver.cpp
```

Note that you may need to run it once or twice before it 'takes'.

Note that you can set a verbose output by setting the command flag `-v`. Verbose output is turned off by default. A verbose output will output information like the QP problem matrices, the messages read from serial, the solution to the QP problem, etc.
```bash
./solver.cpp -v
```

## Configuring the solver

To edit settings like position reference and serial port name, check the configuration files.

To check which serial port the Pi connects to Arduino through, run the following command and update the Serial Port config file accordingly.
```bash
dmesg | grep tty
```
