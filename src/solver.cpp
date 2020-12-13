/**
 * @file solveMPC.cpp
 * @author Luke Schmitt
 * @date October 2020
 */

#define _GLIBCXX_USE_CXX11_ABI 0

#include "ModelPredictiveControlAPI.h"
#include "SerialPort.h"
#include "solver.h"

int main(int argc, char** argv)
{
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "[solveMPC]\tStarting MPC solver." << std::endl;
    std::cout << std::endl;;
    std::cout << std::endl;

    if(cmdOptionExists(argv, argv+argc, "-v"))
    {
        verbose = true;
        std::cout << "[solveMPC]\tVerbose output on." << std::endl;
    }

    // instantiate mpc api object
    ModelPredictiveControlAPI mpc(verbose);

    if(!mpc.solverFlag){return 1;}

    // instantiate serial port object
    SerialPort sp(verbose);

    // enter solver loop
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << "-------------- Entering control loop. --------------" << std::endl;
    std::cout << "----------------------------------------------------" << std::endl;
    std::cout << std::endl;


    while(true)
    {
        if(sp.readPort(mpc.dt, mpc.X))  // check for valid message
        {
            // start = high_resolution_clock::now();

            // update reference, QP problem, and update solution
            if(!mpc.controllerStep()) return 1;

            // print out state and output
            if(mpc.verbose)
            {
                std::cout << "[solveMPC]\t" << "Current state: " << mpc.X.transpose() << std::endl;
                std::cout << "[solveMPC]\t" << "Control output: " << -mpc.U.transpose() << std::endl;
            }
            
            // send control signal to controller
            sp.writePort(-mpc.U);

            // iterate step count
            count++;

            // stop = high_resolution_clock::now();
            
            // std::cout << "[solveMPC]\t" << "Cycle time: " << duration_cast<milliseconds>(stop - start).count() << "ms" << std::endl;
        }
        else
        {
            // send last control signal to controller if we got a bad reading
            sp.writePort(-mpc.U);
        }
    }

    mpc.~ModelPredictiveControlAPI();
    sp.~SerialPort();

    return 0;
}


char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}


bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}