/**
 * @file solveMPC.cpp
 * @author Luke Schmitt
 * @date October 2020
 */

#define _GLIBCXX_USE_CXX11_ABI 0

#include "ModelPredictiveControlAPI.h"
#include "SerialPort.h"
#include "solveMPC.h"

int main(int argc, char** argv)
{
    ss << argv[1];

    if(!(ss >> std::boolalpha >> verbose)) {
        // Parsing error.
    }

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "[solveMPC]\tStarting MPC solver." << std::endl;
    std::cout << std::endl;;
    std::cout << std::endl;

    // instantiate mpc api object
    ModelPredictiveControlAPI mpc(verbose);

    if(!mpc.solverFlag){return 1;}

    // instantiate serial port object
    SerialPort sp("/dev/ttyUSB0");


    // enter loop
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

            if(!mpc.controllerStep()) return 1;

            if(mpc.verbose)
            {
                std::cout << "[solveMPC]\t" << "Current state: " << mpc.X.transpose() << std::endl;
                std::cout << "[solveMPC]\t" << "Control output: " << mpc.U.transpose() << std::endl;
            }
            
            sp.writePort(-mpc.U);

            count++;

            // stop = high_resolution_clock::now();
            
            // std::cout << "[solveMPC]\t" << "Cycle time: " << duration_cast<milliseconds>(stop - start).count() << "ms" << std::endl;

            printf("\n");
        }
        else
        {
            sp.writePort(mpc.U);
        }
        
    }

    sp.~SerialPort();

    return 0;
}
