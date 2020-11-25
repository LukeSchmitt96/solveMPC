/**
 * @file solveMPC.cpp
 * @author Luke Schmitt
 * @date 2020
 */

#define _GLIBCXX_USE_CXX11_ABI 0

#include "solveMPC.h"


void parseArduinoData(char *arr,
                     Eigen::VectorXd &arduinoOutput)
{
    std::cout << "Parsing data from Arduino..." << std::endl;

    std::string str;
    std::cout << "str created" << std::endl;

    // strcpy(str, arr);
    // str = arr;
    // std::cout << "arr set to str" << std::endl;

    // float temp = getDataFromSerial(arr);
    // std::vector<double> split_str = splitToDouble(arr, delimiter);
    // std::vector<double> split_str = splitToDouble(str, delimiter);
    std::cout << "returned from splitToDouble" << std::endl;

    // arduinoOutput(0) = temp[0];
    // arduinoOutput(1) = temp[1];
    // arduinoOutput(2) = temp[2];
    // arduinoOutput(3) = temp[3];
}


int main(int argc, char** argv)
{
    // Arduino output vector in the format [time (ms), x (m), dx (m/s), theta (rad), dtheta (rad/s)]
    Eigen::Matrix<double, 5, 1> arduino_output;

    // The following code implements Linux Serial Port communications detailed in the article:
    // https://blog.mbedded.ninja/programming/operating-systems/linux/linux-serial-ports-using-c-cpp/#overview

    // Open serial port, return error if it fails
    int serial_port = open("/dev/ttyUSB0", O_RDWR);
    if (serial_port < 0) 
    {
        printf("Error %i from open: %s\n", errno, strerror(errno)); return 1;
    }

    // Create new termios struc, we call it 'tty' for convention
    // No need for "= {0}" at the end as we'll immediately write the existing
    // config to this struct
    struct termios tty;

    // Read in existing settings, and handle any error
    // NOTE: This is important! POSIX states that the struct passed to tcsetattr()
    // must have been initialized with a call to tcgetattr() overwise behaviour
    // is undefined
    if(tcgetattr(serial_port, &tty) != 0)
    {
        printf("Error %i from tcgetattr: %s\n", errno, strerror(errno));
    }

    tty.c_cflag &= ~PARENB; // Clear parity bit, disabling parity (most common)
    tty.c_cflag &= ~CSTOPB; // Clear stop field, only one stop bit used in communication (most common)
    tty.c_cflag &= ~CSIZE; // Clear all the size bits, then use one of the statements below
    tty.c_cflag |= CS8; // 8 bits per byte (most common)
    tty.c_cflag &= ~CRTSCTS; // Disable RTS/CTS hardware flow control (most common)
    tty.c_cflag |= CREAD | CLOCAL; // Turn on READ & ignore ctrl lines (CLOCAL = 1)

    tty.c_lflag &= ~ECHO; // Disable echo
    tty.c_lflag &= ~ECHOE; // Disable erasure
    tty.c_lflag &= ~ECHONL; // Disable new-line echo
    tty.c_lflag &= ~ISIG; // Disable interpretation of INTR, QUIT and SUSP

    tty.c_iflag &= ~(IXON | IXOFF | IXANY); // Turn off s/w flow ctrl
    tty.c_iflag &= ~(IGNBRK|BRKINT|PARMRK|ISTRIP|INLCR|IGNCR|ICRNL); // Disable any special handling of received bytes
    tty.c_oflag &= ~OPOST; // Prevent special interpretation of output bytes (e.g. newline chars)
    tty.c_oflag &= ~ONLCR; // Prevent conversion of newline to carriage return/line feed
    // tty.c_oflag &= ~OXTABS; // Prevent conversion of tabs to spaces (NOT PRESENT IN LINUX)
    // tty.c_oflag &= ~ONOEOT; // Prevent removal of C-d chars (0x004) in output (NOT PRESENT IN LINUX)

    tty.c_cc[VTIME] = 1;    // Wait for up to 1s (10 deciseconds), returning as soon as any data is received.
    tty.c_cc[VMIN] = 0;

    // Set in/out baud rate to be 9600
    cfsetispeed(&tty, B9600);
    cfsetospeed(&tty, B9600);

    // Save tty settings, also checking for error
    if (tcsetattr(serial_port, TCSANOW, &tty) != 0) 
    {
        printf("Error %i from tcsetattr: %s\n", errno, strerror(errno));
    }

    // unsigned char msg[] = { 'H', 'e', 'l', 'l', 'o', '\r' };
    // write(serial_port, msg, sizeof(msg));

    // Allocate memory for read buffer, set size according to your needs
    char* read_buf = new char[256];

    // while (1)
    // {

    //     // Normally you wouldn't do this memset() call, but since we will just receive
    //     // ASCII data for this example, we'll set everything to 0 so we can
    //     // call printf() easily.
    //     memset(&read_buf, '\0', sizeof(read_buf));

    //     // Read bytes. The behaviour of read() (e.g. does it block?,
    //     // how long does it block for?) depends on the configuration
    //     // settings above, specifically VMIN and VTIME
    //     int num_bytes = read(serial_port, &read_buf, sizeof(read_buf));

    //     // n is the number of bytes read. n may be 0 if no bytes were received, and can also be negative to signal an error.
    //     if (num_bytes < 0)
    //     {
    //         printf("Error reading: %s", strerror(errno)); return 1;
    //     } 

    //     // printf("Read %i bytes.", num_bytes);
    //     // printf("Read %i bytes. Received message: %s", num_bytes, read_buf);

    //     Eigen::VectorXd arduino_output;
    //     parseArduinoData(read_buf, arduino_output);

    //     // std::cout << arduino_output << std::endl << std::endl;

    //     // printf("Read %i bytes. Duration: %d. Received message: %s", num_bytes, duration_milliseconds, read_buf);
    //     // printf("%s", read_buf);

    // }

    // allocate system matrices
    Eigen::Matrix<double, N_S, N_S>                     Ad;
    Eigen::Matrix<double, N_S, N_C>                     Bd;
    Eigen::Matrix<double, N_O, N_S>                     Cd;
    Eigen::Matrix<double, N_O, N_C>                     Dd;
    Eigen::Matrix<double, 2, N_C>                       G;
    Eigen::Matrix<double, 2, N_C>                       S; 
                                                        S = G;
    Eigen::Matrix<double, 2*mpcWindow, 4>               Sbar;
    Eigen::Matrix<double, 2, 1>                         W0; 
                                                        W0 << 1, 1;
    Eigen::Matrix<double, 2*mpcWindow, 1>               W;
    Eigen::SparseMatrix<double>                         Gbar; // 2*mpcWindow, 4*mpcWindow

    // allocate the weight matrices
    Eigen::Matrix<double, N_S, N_S>                     Q;    
    Eigen::Matrix<double, N_C, N_C>                     R;    
    Eigen::Matrix<double, 1,   1>                       RD;
    Eigen::Matrix<double, N_S*mpcWindow, N_S*mpcWindow> Qbar;    
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Rbar;
    Eigen::Matrix<double, 1*mpcWindow,   1*mpcWindow>   RbarD;

    // allocate Sx, Su, Su1, CAB matrices
    Eigen::Matrix<double, N_S*mpcWindow, N_S>           Sx;
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Su;
    Eigen::Matrix<double, N_C*mpcWindow, 1>             Su1;
    Eigen::Matrix<double, N_S*mpcWindow, N_S>           CAB;

    // allocate the initial and the reference state space
    Eigen::Matrix<double, N_S, 1>                       X;
    Eigen::Matrix<double, N_O, mpcWindow>               ref;

    // allocate QP problem matrices and vectors
    Eigen::SparseMatrix<double>                         H;
    Eigen::Matrix<double, N_C*mpcWindow, 1>             Fu;
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Fr;
    Eigen::Matrix<double, N_C*mpcWindow, N_S>           Fx;
    Eigen::Matrix<double, N_C*mpcWindow, 1>             f;
    Eigen::Matrix<double, 2*mpcWindow, 1>               lb;
    Eigen::Matrix<double, 2*mpcWindow, 1>               ub;

    // set arbitraty (to be filled) variables
    X << 0, 0, 0, 0;
    ref = Eigen::Matrix<double, N_O, mpcWindow>::Zero();
    lb = Eigen::Matrix<double, 2*mpcWindow, 1>::Ones() * Eigen::Infinity;
    Eigen::Matrix<double, 1, mpcWindow> t;
    double t0 = 0;

    // set MPC problem quantities
    setSystemVars(Ad,Bd,Cd,Dd);
    setQ_R_RD(Q, R, RD, weightQ, weightR, weightRD);
    computeQbar_Rbar_RbarD(Qbar, Rbar, RbarD, Q, R, RD);
    computeSx_Su_Su1_CAB(Sx, Su, Su1, CAB, Ad, Cd, Bd);
    computeH(H, Rbar, Su, Qbar);
    setFVars(Fu, Fr, Fx, Qbar, Rbar, Sx, Su, Su1); 
    computeGbar_Sbar_W(Gbar, Sbar, W, G, S, W0);
    setF(f, Fu, Fr, Fx, X, ref);
    
    ub = W+Sbar*X;

    std::cout << "[INFO] All matrices built successfully." << std::endl;

    // instantiate the solver
    OsqpEigen::Solver solver;

    // settings
    solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // set up the QP solver
    solver.data()->setNumberOfVariables(N_S*mpcWindow);
    solver.data()->setNumberOfConstraints(2*mpcWindow);
    if(!solver.data()->setHessianMatrix(H)) return 1;
    if(!solver.data()->setGradient(f)) return 1;
    if(!solver.data()->setLinearConstraintsMatrix(Gbar)) return 1;
    if(!solver.data()->setLowerBound(lb)) return 1;
    if(!solver.data()->setUpperBound(ub)) return 1;

    // instantiate the solver
    if(!solver.initSolver()) return 1;

    // controller input and QPSolution vector
    Eigen::Vector4d U;
    Eigen::VectorXd QPSolution;


    // enter loop
    std::cout << std::endl;
    std::cout << "[INFO] ----------------------" << std::endl << std::endl;
    std::cout << "[INFO] Entering control loop." << std::endl << std::endl;
    std::cout << "[INFO] ----------------------" << std::endl << std::endl;

    auto start = high_resolution_clock::now();
    auto stop  = high_resolution_clock::now();

    int count = 0;
    while(true)
    {
        // start = high_resolution_clock::now();

        // get state by reading from serial
        // getDataFromSerial(t0, X);
        
        // update reference trajectory
        t = linspace(t0,t0+0.005*(mpcWindow-1),mpcWindow);
        // std::cout << t << std::endl << std::endl << std::endl;
        
        updateRef(ref, t);

        ref = Eigen::Matrix<double, 4, 40>::Zero();

        // update gradient and constraints
        setF(f, Fu, Fr, Fx, X, ref);
        if(!solver.updateGradient(f)) return 1;
        if(!solver.updateUpperBound(W+Sbar*X)) return 1;

        // solve the QP problem
        if(!solver.solve()) return 1;

        
        // get the controller input
        QPSolution = solver.getSolution();
        U = QPSolution.block<N_C, 1>(0, 0);
        // std::cout << U << std::endl;
        // sendDataToSerial(U);

        
        // save data into file
        // auto x0Data = x0.data();

        // auto stop = high_resolution_clock::now();

        // std::cout << duration_cast<microseconds>(stop-start).count()/1000.0 << std::endl;
        count++;
    }

    close(serial_port);

    return 0;
}
