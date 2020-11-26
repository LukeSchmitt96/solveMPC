Need to install dependencies of osqp-eigen. Find information [here](https://github.com/robotology/osqp-eigen).

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