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
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "[solveMPC]\tStarting MPC solver." << std::endl;
    std::cout << std::endl;;
    std::cout << std::endl;

    // instantiate mpc api object
    ModelPredictiveControlAPI mpc;

    // instantiate serial port object
    SerialPort sp("/dev/ttyUSB0");

    // instantiate the solver
    OsqpEigen::Solver solver;

    // solver settings
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // set up the QP solver
    solver.data()->setNumberOfVariables(mpc.n_variables);
    solver.data()->setNumberOfConstraints(mpc.n_constraints);
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
    char read_buf [64];

    // while(true){
    //     sp.sendInit();
    // }

    int count = 0;

    // enter loop
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "-------------- Entering control loop. --------------" << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << std::endl;


    while(true)
    {

        // get state by reading from serial
        if(sp.readPort(mpc.dt, mpc.X))  // check for valid message
        {
            // start = high_resolution_clock::now();

            mpc.t0 += mpc.dt;

            // update reference trajectory
            mpc.t = mpc.linspace(mpc.t0, mpc.t0+mpc.dt/1000.0*(mpcWindow-1), mpcWindow);
            mpc.updateRef(mpc.ref, mpc.t);

            // update gradient and constraints
            mpc.setF(mpc.f, mpc.Fu, mpc.Fr, mpc.Fx, mpc.X, mpc.ref);
            if(!solver.updateGradient(mpc.f)) return 1;
            if(!solver.updateUpperBound(mpc.W+mpc.Sbar*mpc.X)) return 1;

            // solve the QP problem
            if(!solver.solve()) return 1;
            
            // get the controller input
            QPSolution = solver.getSolution();

            // TODO: add gain to the U and the X term
            // TODO: add ff term, just control signal propogated through the system

            U = QPSolution.block<N_C, 1>(0, 0);

            // Xact = mpc.Ad*mpc.X + mpc.Bd*U

            std::cout << "[solveMPC]\t" << "Control output: " << U.transpose() << std::endl;

            sp.writePort(U);

            count++;

            // stop = high_resolution_clock::now();
            
            // std::cout << "[solveMPC]\t" << "Cycle time: " << duration_cast<milliseconds>(stop - start).count() << "ms" << std::endl;
        }
        else
        {
            sp.writePort(U);
        }
        
    }

    sp.~SerialPort();

    return 0;
}
