#pragma once
#ifndef HEADER_H
#define HEADER_H

#include <chrono>
using namespace std::chrono;

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <math.h>
#include <chrono>

#include <string.h>
#include <stdio.h>
#include <fcntl.h>      // Contains file controls like O_RDWR
#include <errno.h>      // Error integer and strerror() function
#include <termios.h>    // Contains POSIX terminal control definitions
#include <unistd.h>     // write(), read(), close()

// osqp-eigen
#include "OsqpEigen/OsqpEigen.h"

// eigen and eigen helpers
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>


/* ------------------------------------------------------------ */
/* ----------------- Parameters and Helpers ------------------- */
/* ------------------------------------------------------------ */

// // String splitting delimiter
// const std::string delimiter = ", ";

// set the MPC preview window
const int mpcWindow = 40;

// set matrix dimensions
const int N_S = 4;  // number of states
const int N_C = 4;  // number of controls
const int N_O = 4;  // number of outputs

// set weights
const double weightQ =  1.0;    // state weight
const double weightR =  1e-3;   // control weight
const double weightRD = 0.0;    // change in control weight

// timing parameters
const double Ts = 100;          // square wave period


/* ------------------------------------------------------------ */
/* ----------------- Eigen Helper Functions ------------------- */
/* ------------------------------------------------------------ */

Eigen::MatrixXd blkdiag(const Eigen::MatrixXd& a, int count)
{
    Eigen::MatrixXd bdm = Eigen::MatrixXd::Zero(a.rows() * count, a.cols() * count);
    for (int i = 0; i < count; ++i)
    {
        bdm.block(i * a.rows(), i * a.cols(), a.rows(), a.cols()) = a;
    }

    return bdm;
}


template<typename T>
Eigen::Matrix<double, 1, mpcWindow> linspace(T start_in, T end_in, int num_in)
{

    Eigen::Matrix<double, 1, mpcWindow> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced(0,i) = start + delta * i;
    }
    linspaced(0,num-1) = end;   // I want to ensure that start and end
                                // are exactly the same as the input
        
    return linspaced;
}


/* ------------------------------------------------------------ */
/* ------------------- MPC Helper Functions ------------------- */
/* ------------------------------------------------------------ */

void setSystemVars(Eigen::Matrix<double,N_S,N_S> &Ad,
                   Eigen::Matrix<double,N_S,N_C> &Bd,
                   Eigen::Matrix<double,N_O,N_S> &Cd,
                   Eigen::Matrix<double,N_O,N_C> &Dd)
{
    Ad <<   1.0001,     0.0050,     -0.0038,    -0.0002,
            0.0162,     1.0064,     -0.9551,    -0.0601,
            0.0011,     0.0004,      0.9378,     0.0009,
            0.2738,     0.1080,     -15.0372,   -0.0141;

    Bd <<  -0.0001,    -0.0050,      0.0046,     0.0003,
           -0.0022,    -0.0009,      0.1868,     0.0122,
           -0.0013,    -0.0005,      0.0724,    -0.0001,
           -0.0382,    -0.0162,      2.1065,     0.2054;

    Cd = Eigen::Matrix<double,N_O,N_S>::Identity();

    // Dd = Eigen::Matrix<double,N_O,N_C>::Zero();
    Dd <<   -0.0000,   -0.0025,      0.0016,    0.0001, 
            -0.0134,   -0.0053,      0.7694,    0.0482, 
            -0.0005,   -0.0002,      0.0266,   -0.0008, 
            -0.2266,   -0.0889,     12.4390,    0.8142; 

    std::cout << "[INFO] System variables created."   << std::endl;
    // std::cout << "Ad:" << std::endl << Ad << std::endl << std::endl;
    // std::cout << "Bd:" << std::endl << Bd << std::endl << std::endl;
    // std::cout << "Cd:" << std::endl << Cd << std::endl << std::endl;
    // std::cout << "Dd:" << std::endl << Dd << std::endl << std::endl;
}


void setQ_R_RD(Eigen::Matrix<double,N_S,N_S>    &Q,
               Eigen::Matrix<double,N_C,N_C>    &R,
               Eigen::Matrix<double,1,1>        &RD,
               double wq,
               double wr,
               double wrd)
{
    Q  = Eigen::Matrix<double, N_S, N_S>::Identity() * wq;
    R  = Eigen::Matrix<double, N_C, N_C>::Identity() * wr;
    RD = Eigen::Matrix<double,1,1>::Identity() * wrd;

    std::cout << "[INFO] Set Q, R, and RD matrices created."  << std::endl;
    // std::cout << "Q:"  << std::endl << Q  << std::endl << std::endl;
    // std::cout << "R:"  << std::endl << R  << std::endl << std::endl;
    // std::cout << "RD:" << std::endl << RD << std::endl << std::endl;
}


void computeQbar_Rbar_RbarD(Eigen::Matrix<double,N_S*mpcWindow,N_S*mpcWindow>   &Qbar,
                            Eigen::Matrix<double,N_C*mpcWindow,N_C*mpcWindow>   &Rbar,
                            Eigen::Matrix<double,1*mpcWindow,1*mpcWindow>       &RbarD,
                            Eigen::Matrix<double,N_S,N_S>                       Q,
                            Eigen::Matrix<double,N_C,N_C>                       R,
                            Eigen::Matrix<double,1,1>                           RD)
{
    Qbar = blkdiag(Q, mpcWindow);
    Rbar = blkdiag(R, mpcWindow);
    RbarD = blkdiag(RD, mpcWindow);

    std::cout << "[INFO] Lifted weight matrices created."   << std::endl;
}


void computeSx_Su_Su1_CAB(Eigen::Matrix<double, N_S*mpcWindow, N_S>             &Sx,
                          Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>   &Su,
                          Eigen::Matrix<double, N_C*mpcWindow, 1>               &Su1,
                          Eigen::Matrix<double, N_S*mpcWindow, N_S>             &CAB,
                          Eigen::Matrix<double, N_S, N_S>                       Ad,
                          Eigen::Matrix<double, N_S, N_C>                       Bd,
                          Eigen::Matrix<double, N_O, N_S>                       Cd)
{
    Sx  = Eigen::MatrixXd::Zero(Sx.rows(),  Sx.cols());
    Su  = Eigen::MatrixXd::Zero(Su.rows(),  Su.cols());
    Su1 = Eigen::MatrixXd::Zero(Su1.rows(), Su1.cols());
    CAB = Eigen::MatrixXd::Zero(CAB.rows(), CAB.cols());

    // compute Sx and CAB
    for(int i=0; i<mpcWindow; i++)
    {
        Sx.block<N_S, N_O>(i*N_O,0)  = Cd*Ad.pow(i+1);
        CAB.block<N_S,N_O>(i*N_O,0) = Cd*Ad.pow(i)*Bd;
    }

    Eigen::VectorXd CAB_Vector(Eigen::Map<Eigen::VectorXd>(CAB.data(), CAB.cols()*CAB.rows()));

    // compute Su
    for(int i=0; i<Su.rows(); i++)
    {
        for(int j=0; j<i; j++)
        {
            Su(i,j) = CAB_Vector.segment(0,i-j).sum();
        }
    }

    Su1 = Su.col(0);

    std::cout << "[INFO] Sx Su, Su1, CAB created"   << std::endl;
}


void computeH(Eigen::SparseMatrix<double>                         &H,
              Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Rbar,
              Eigen::Matrix<double, N_S*mpcWindow, N_S*mpcWindow> Qbar,
              Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Su)
{
    Eigen::MatrixXd H_temp_1, H_temp_2;
    H_temp_1 = 2*(Rbar + Su.transpose() * Qbar * Su);
    // H_temp_2 = (H_temp_1+H_temp_1.transpose())/2.0;
    H_temp_2 = H_temp_1;

    H.resize(N_S*mpcWindow, N_S*mpcWindow); 
 
    for(int i=0; i<mpcWindow*N_S; i++)
    {
        for(int j=0; j<mpcWindow*N_S; j++)
        {
            H.insert(i,j) = H_temp_2(i,j);
        }
    }



    H.makeCompressed();

    std::cout << "[INFO] Hessian H created." << std::endl;
}


void setFVars(Eigen::Matrix<double, N_C*mpcWindow, 1>               &Fu,
              Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>   &Fr, 
              Eigen::Matrix<double, N_C*mpcWindow, N_S>             &Fx, 
              Eigen::Matrix<double, N_S*mpcWindow, N_S*mpcWindow>   Qbar,
              Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>   Rbar, 
              Eigen::Matrix<double, N_S*mpcWindow, N_S>             Sx,
              Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>   Su, 
              Eigen::Matrix<double, N_C*mpcWindow, 1>               Su1)
{
    Fu = 2*(Rbar.diagonal().transpose() + Su1.transpose() * Qbar * Su).transpose();
    Fr = -2*(Qbar*Su).transpose();
    Fx = 2*(Sx.transpose() * Qbar * Su).transpose();

    std::cout << "[INFO] Components of F created."   << std::endl;
}


void computeGbar_Sbar_W(Eigen::SparseMatrix<double>                     &Gbar,
                        Eigen::Matrix<double, 2*mpcWindow, N_C>         &Sbar,
                        Eigen::Matrix<double, 2*mpcWindow, 1>           &W,
                        Eigen::Matrix<double, 2, N_C>                   G,
                        Eigen::Matrix<double, 2, N_C>                   S,
                        Eigen::Matrix<double, 2, 1>                     W0)
{
    Eigen::Matrix<double, 2*mpcWindow, N_C*mpcWindow> Gbar_temp;
    Gbar_temp = blkdiag(G, mpcWindow);
    Gbar.resize(2*mpcWindow, 4*mpcWindow);

    for(int i=0; i<2*mpcWindow; i++)
    {
        for(int j=0; j<N_C*mpcWindow; j++)
        {
            Gbar.insert(i,j) = Gbar_temp(i,j);
        }
    }

    Gbar.makeCompressed();

    for(int i; i<mpcWindow; i++)
    {
        Sbar.block<2, N_C> (i*S.rows(),0) = S;
        W.block<2, 1>(i*W0.rows(),0) = W0;
    }

    std::cout << "[INFO] Constraints created."   << std::endl;
}


void setF(Eigen::Matrix<double, N_C*mpcWindow, 1>               &f,
          Eigen::Matrix<double, N_C*mpcWindow, 1>               Fu,
          Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>   Fr, 
          Eigen::Matrix<double, N_C*mpcWindow, N_S>             Fx, 
          Eigen::Matrix<double, N_S, 1>                         X,
          Eigen::Matrix<double, N_S, mpcWindow>                 ref)
{
    Eigen::Map<Eigen::Matrix<double, N_C*mpcWindow, 1>> ref_vector(ref.data(), ref.size());
    f = Fx*X+Fr*ref_vector;
}


void updateRef(Eigen::Matrix<double, N_S, mpcWindow>            &ref, 
               Eigen::Matrix<double, 1, mpcWindow>              t)
{
    ref = Eigen::Matrix<double, N_C, mpcWindow>::Zero();

    Eigen::Matrix<double, 1, mpcWindow> ref_position;
    Eigen::Matrix<double, 1, mpcWindow> temp;
    
    temp = t.unaryExpr([](const double x) 
        { 
            return fmod(x/(Ts/2),2);
        }
    );

    ref_position = temp.unaryExpr([](const double x)
        {
            return x <= 1.0 ? 1.0 : 0.0;
        }
    );

    ref.block<1,mpcWindow>(0,0) = ref_position;
}


/* ------------------------------------------------------------ */
/* ------------------ OSQP Helper Functions ------------------- */
/* ------------------------------------------------------------ */

const int buffSize = 64;
char inputSeveral[buffSize]; // space for 31 chars and a terminator

int maxChars = 63; // a shorter limit to make it easier to see what happens
                           //   if too many chars are entered


float getDataFromSerial(char input[64]) 
{
    float ref[4] = {0.0, 0.0, 0.0, 0.0};
    
    // Function credits: https://forum.arduino.cc/index.php?topic=236162.0

    // this function takes the characters from the serial input and converts them
    // to a single floating point value using the function "atof()"
     
    // a similar approach can be used to read an integer value if "atoi()" is used

    // first read severalChars into the array inputSeveral
    // inputSeveral[0] = 0; 
    int index = 0;
    char *ptr = NULL;
    double temp;

    ptr = strtok(inputSeveral, " ");  // takes a list of delimiters
    while(index < 4)
    {
        temp = atof(ptr);
        if (temp != 0) 
        {
            ref[index] = temp;
        }
        index++;
        ptr = strtok(NULL, " ");  // takes a list of delimiters
    }
    
    std::cout << ref[0];
    std::cout << " ";
    std::cout << ref[1];
    std::cout << " ";
    std::cout << ref[2];
    std::cout << " ";
    std::cout << ref[3];
    std::cout << "\n";
    // and then convert the string into a floating point number

    // return ref;

}



/* ------------------------------------------------------------ */
/* ----------------- Comms Helper Functions ------------------- */
/* ------------------------------------------------------------ */


std::vector<double> convertStringVectortoDoubleVector(const std::vector<std::string>&);

std::vector<double> splitToDouble(std::string str, std::string token)
{
    std::cout << "in splitToDouble" << std::endl;

    double d;
    std::vector<std::string> result;
    while(str.size()){
        int index = str.find(token);
        std::cout << "index: " << index << std::endl;

        if(index!=std::string::npos)
        {
            result.push_back(str.substr(0,index));
            str = str.substr(index+token.size());
            if(str.size()==0)result.push_back(str);
        }
        else
        {
            result.push_back(str);
            str = "";
        }
    }
    return convertStringVectortoDoubleVector(result);
}


std::vector<double> convertStringVectortoDoubleVector(const std::vector<std::string>& stringVector)
{
    std::vector<double> doubleVector(stringVector.size());
    std::transform(stringVector.begin(), stringVector.end(), doubleVector.begin(), [](const std::string& val)
                 {
                     return stod(val);
                 });
    return doubleVector;
}



#endif // HEADER_H