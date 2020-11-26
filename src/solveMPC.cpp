/**
 * @file solveMPC.cpp
 * @author Luke Schmitt
 * @date 2020
 */

#define _GLIBCXX_USE_CXX11_ABI 0

#include <chrono>

#include "OsqpEigen/OsqpEigen.h"

#include "ModelPredictiveControlAPI.h"
#include "SerialPort.h"

using namespace std::chrono;

int main(int argc, char** argv)
{
    // instantiate mpc api object
    ModelPredictiveControlAPI mpc;

    // instantiate serial port object
    SerialPort sp("/dev/ttyUSB0");

    // instantiate the solver
    OsqpEigen::Solver solver;

    // settings
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // set up the QP solver
    solver.data()->setNumberOfVariables(N_S*mpcWindow);
    solver.data()->setNumberOfConstraints(2*mpcWindow);
    if(!solver.data()->setHessianMatrix(mpc.H)) return 1;
    if(!solver.data()->setGradient(mpc.f)) return 1;
    if(!solver.data()->setLinearConstraintsMatrix(mpc.Gbar)) return 1;
    if(!solver.data()->setLowerBound(mpc.lb)) return 1;
    if(!solver.data()->setUpperBound(mpc.ub)) return 1;

    // instantiate the solver
    if(!solver.initSolver()) return 1;

    // controller input and QPSolution vector
    Eigen::Vector4d U;
    Eigen::VectorXd QPSolution;

    // set timinig variables
    auto start = high_resolution_clock::now(), stop = high_resolution_clock::now();
    
    // set serial read variables
    int num_bytes; 
    char read_buf [32];


    int count = 0;

    // enter loop
    std::cout << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "-------------- Entering control loop. --------------" << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << std::endl;

    while(true)
    {

        // get state by reading from serial
        sp.readPort(mpc.dt, mpc.X);

        mpc.t0 += mpc.dt;

        // update reference trajectory
        mpc.t = mpc.linspace(mpc.t0, mpc.t0+mpc.dt/1000*(mpcWindow-1), mpcWindow);
        
        mpc.updateRef(mpc.ref, mpc.t);

        // update gradient and constraints
        mpc.setF(mpc.f, mpc.Fu, mpc.Fr, mpc.Fx, mpc.X, mpc.ref);

        if(!solver.updateGradient(mpc.f)) return 1;
        if(!solver.updateUpperBound(mpc.W+mpc.Sbar*mpc.X)) return 1;

        // solve the QP problem
        if(!solver.solve()) return 1;
        
        // get the controller input
        QPSolution = solver.getSolution();

        sp.writePort(QPSolution.block<N_C, 1>(0, 0));

        count++;
    }

    sp.~SerialPort();

    return 0;
}
