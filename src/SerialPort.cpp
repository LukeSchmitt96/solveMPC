#include "./SerialPort.h"

// C library headers
#include <stdio.h>
#include <string.h>
#include <iostream>

// Linux headers
#include <fcntl.h> // Contains file controls like O_RDWR
#include <errno.h> // Error integer and strerror() function
#include <termios.h> // Contains POSIX terminal control definitions
#include <unistd.h> // write(), read(), close()


SerialPort::SerialPort(bool verbose_)
{
    std::cout << "[SerialPort]\tSerial Port object created." << std::endl;

    verbose = verbose_;

    // get stream from config file and parse
    std::ifstream file("./config/Serial_Port.json");
    cfg = json::parse(file);

    // get port name from cfg
    port = cfg["port"].get<std::string>();

    if(verbose)
    {
        std::cout << "[SerialPort]\tAttempting connection to serial port "
                << port
                << "\"." << std::endl;
    }

    // try to open the port
    int attempt = 0;
    serial_port = open(port.c_str(), O_RDWR);
    while(serial_port < 0) 
    {
        // // stop trying after max_attempts
        // if(attempt>4)
        // {
        //     printf("[SerialPort]\t Could not open serial port after %d attempts. Check connection. Shutting down...",attempt);
        //     return;
        // }

        printf("[SerialPort]\tError %i from open: %s. Reattempting connection...\n", errno, strerror(errno));
        sleep(1);
        serial_port = open(port.c_str(), O_RDWR);
        attempt++;
    }

    printf("[SerialPort]\tSerial port opened successfully.\n");
     
    // much of the following was taken from the guide "Linux Serial Ports Using C/C++" found at
    //      https://blog.mbedded.ninja/programming/operating-systems/linux/linux-serial-ports-using-c-cpp/#overview

    // Read in existing settings, and handle any error
    // NOTE: This is important! POSIX states that the struct passed to tcsetattr()
    // must have been initialized with a call to tcgetattr() overwise behaviour
    // is undefined
    if(tcgetattr(serial_port, &tty) != 0)
    {
        printf("[SerialPort]\tError %i from tcgetattr: %s", errno, strerror(errno));
    }

    tty.c_cflag &= ~PARENB;         // Clear parity bit, disabling parity (most common)
    tty.c_cflag &= ~CSTOPB;         // Clear stop field, only one stop bit used in communication (most common)
    tty.c_cflag &= ~CSIZE;          // Clear all the size bits, then use one of the statements below
    tty.c_cflag |= CS8;             // 8 bits per byte (most common)
    tty.c_cflag &= ~CRTSCTS;        // Disable RTS/CTS hardware flow control (most common)
    tty.c_cflag |= CREAD | CLOCAL;  // Turn on READ & ignore ctrl lines (CLOCAL = 1)
    tty.c_lflag &= ~ECHO;           // Disable echo
    tty.c_lflag &= ~ECHOE;          // Disable erasure
    tty.c_lflag &= ~ECHONL;         // Disable new-line echo
    tty.c_lflag &= ~ISIG;           // Disable interpretation of INTR, QUIT and SUSP
    tty.c_iflag &= ~(IXON | IXOFF | IXANY); // Turn off s/w flow ctrl
    tty.c_iflag &= ~(IGNBRK|BRKINT|PARMRK|ISTRIP|INLCR|IGNCR|ICRNL); // Disable any special handling of received bytes
    tty.c_oflag &= ~OPOST; // Prevent special interpretation of output bytes (e.g. newline chars)
    tty.c_oflag &= ~ONLCR; // Prevent conversion of newline to carriage return/line feed
    // tty.c_oflag &= ~OXTABS; // Prevent conversion of tabs to spaces (NOT PRESENT IN LINUX)
    // tty.c_oflag &= ~ONOEOT; // Prevent removal of C-d chars (0x004) in output (NOT PRESENT IN LINUX)

    tty.c_cc[VTIME] = 30;   // Wait for up to 1s (10 deciseconds), returning as soon as any data is received.
    tty.c_cc[VMIN] = 30;

    // Set in/out baud rate. Settings in config file
    cfsetispeed(&tty, get_baud(cfg["baud"].get<int>()));
    cfsetospeed(&tty, get_baud(cfg["baud"].get<int>()));

    // Save tty settings, also checking for error
    if (tcsetattr(serial_port, TCSANOW, &tty) != 0) 
    {
        printf("[SerialPort]\tError %i from tcsetattr: %s\n", errno, strerror(errno));
    }
}


SerialPort::~SerialPort()
{
    printf("[SerialPort]\tDestructing Serial Port Communication Object.\n");
    close(serial_port);
}


void SerialPort::getDataFromSerial(double &dt, Eigen::Matrix<double, N_S, 1> &X, char* read_buf) 
{
    float ref[5] = {0.0000, 0.0000, 0.0000, 0.0000, 0.0000};

    // Function credits: https://forum.arduino.cc/index.php?topic=236162.0

    // this function takes the characters from the serial input and converts them
    // to a single floating point value using the function "atof()"
     
    // a similar approach can be used to read an integer value if "atoi()" is used

    int index = 0;
    char *ptr = NULL;
    double temp;

    ptr = strtok(read_buf, " ");  // takes a list of delimiters
    while(index < 5)
    {
        temp = atof(ptr);
        if (temp != 0) 
        {
            ref[index] = temp;
        }
        index++;
        ptr = strtok(NULL, " ");  // takes a list of delimiters
    }
    
    dt   = ref[0];

    X(0) = ref[1];
    X(1) = ref[2];
    X(2) = ref[3];
    X(3) = ref[4];
}


bool SerialPort::readPort(double            dt,
                          Eigen::Vector4d   &X)
{
    memset(&read_buf, '\0', sizeof(read_buf));
    num_bytes = read(serial_port, &read_buf, sizeof(read_buf));

    if (num_bytes > 30)         // TODO: tune this number for best performance
    {
        if(verbose){printf("\n[SerialPort]\tRead %i bytes. Received message: %s", num_bytes, read_buf);}
        getDataFromSerial(dt, X, read_buf);
        return true;
    } 
    else
    {
        // if(verbose){printf("[SerialPort]\tBad serial read, reusing last control signal.\n");}
        return false;
    }
}


void SerialPort::writePort(Eigen::MatrixXd data_send)
{
    // write(serial_port, data_send.data(), sizeof(data_send));
    write(serial_port, std::to_string(data_send(0)).c_str(), sizeof(std::to_string(data_send(0)).c_str()));
}


void SerialPort::sendInit()
{
    
}

int SerialPort::get_baud(int baud)
{
    switch (baud) {
    case 9600:
        return B9600;
    case 19200:
        return B19200;
    case 38400:
        return B38400;
    case 57600:
        return B57600;
    case 115200:
        return B115200;
    case 230400:
        return B230400;
    case 460800:
        return B460800;
    case 500000:
        return B500000;
    case 576000:
        return B576000;
    case 921600:
        return B921600;
    case 1000000:
        return B1000000;
    case 1152000:
        return B1152000;
    case 1500000:
        return B1500000;
    case 2000000:
        return B2000000;
    case 2500000:
        return B2500000;
    case 3000000:
        return B3000000;
    case 3500000:
        return B3500000;
    case 4000000:
        return B4000000;
    default: 
        return -1;
    }
}
