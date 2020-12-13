#pragma once

#include <string>

#include <fcntl.h>      // Contains file controls like O_RDWR
#include <errno.h>      // Error integer and strerror() function
#include <termios.h>    // Contains POSIX terminal control definitions
#include <unistd.h>     // write(), read(), close()

#include "ModelPredictiveControlAPI.h"

// eigen and eigen helpers
#include <Eigen/Dense>

// json api
#include <nlohmann/json.hpp>

class SerialPort
{
public:
    /** constructor
     * 
     * @param verbose     sets verbose output mode
     */
    SerialPort(bool);

    // destructor
    ~SerialPort();

    /**
     * converts serial baudrate from config
     * 
     * @param baud value from config
     */
    int get_baud(int);

    /**
     * read data from a serial connection
     * 
     * @param dt            pointer to dt
     * 
     * @param X             pointer to state vector
     * 
     * @param read_buf      buffer from serial
     */ 
    void getDataFromSerial(double                           &, 
                           Eigen::Matrix<double, N_S, 1>    &,
                           char*                            ); 

    /**
     * get data from serial port
     * 
     * @param dt            pointer to time diff
     * 
     * @param X             pointer to state vector
     * 
     * @returns             status message:
     *                          0 if read successfully
     *                         -1 if no message was read
     */                           
    bool readPort(double, Eigen::Vector4d &);

    /**
     * write data to port
     * 
     * @param U             control signal
     */
    void writePort(Eigen::MatrixXd);

    /**
     * send init byte to arduino
     */
    void sendInit();
    
    bool verbose;

    int serial_port;
    struct termios tty;

    // set serial read variables
    int num_bytes; 
    char read_buf [42];

    // config variables
    nlohmann::json cfg;

    std::string port;
};