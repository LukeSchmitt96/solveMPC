#pragma once

#include <string>

#include <fcntl.h>      // Contains file controls like O_RDWR
#include <errno.h>      // Error integer and strerror() function
#include <termios.h>    // Contains POSIX terminal control definitions
#include <unistd.h>     // write(), read(), close()

#include "ModelPredictiveControlAPI.h"

#include <Eigen/Dense>

class SerialPort
{
public:
    /** constructor
     * 
     * @param port_name     name of port to open for serial communication
     */
    SerialPort(const char*);

    // destructor
    ~SerialPort();

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
     */                           
    void readPort(double, Eigen::Vector4d &);

    /**
     * write data to port
     * 
     * @param U             control signal
     */
    void writePort(Eigen::Vector4d);

    int serial_port;
    struct termios tty;

    // set serial read variables
    int num_bytes; 
    char read_buf [32];

};