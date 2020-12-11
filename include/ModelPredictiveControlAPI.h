#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cctype>
#include <math.h>

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
const int N_O = 1;  // number of outputs

// timing parameters
const double Ts = 100;          // square wave period

// MPC API Class
class ModelPredictiveControlAPI 
{

public:
    // constructor
    ModelPredictiveControlAPI();

    // destructor
    ~ModelPredictiveControlAPI();
    
    /**
     * Set the system values. These should be the discrete-time matrices
     *      found in MATLAB after c2d is run
     */
    void setSystemVars();

    /**
     * Set the state, control, and change in control weight matrices based
     *      on weight inputs
     */
    void setQ_R_RD();

    /**
     * Set the lifted state, control, and change in control weight matrices
     */
    void computeQbar_Rbar_RbarD();

    /**
     * Complete lifted dynamics transform matrices
     */
    void computeSx_Su_Su1_CAB();       

    /**
     * Compute hessian for the QP problem
     */    
    void computeH();

    /**
     * Compute components for the qp gradient  
     */
    void setFVars();    

    /**
     * Compute S
     */
    void computeS();      

    /**
     * Compute Sbar
     */
    void computeSbar(); 

    /**
     * Compute G
     */
    void computeG();

    /**
     * Compute constraints components for the qp problem
     */
    void computeGbar();          

    /**
     * compute gradient for the qp problem
     */    
    void setF();

    /**
     * update the reference given the time
     */
    void updateRef();

    /**
     * set LL matrix
     */
    void setLL();      
    
    /**
     * set Lu matrix
     */
    void setLu();               

    /**
     * linearizes nonlinear model into a linear model about X
     * 
     * @param ts        sampling time
     * 
     * @param U_LQR     control from LQR
     */
    void linearizeABCD(double ts, double U_LQR);

    /**
     * converts a continuous time system to discrete time
     * 
     * @param ts        sampling time
     */
    void c2d(double ts);

    // system matrices - Updated 
    Eigen::Matrix<double, N_S, N_S>                     A;
    Eigen::Matrix<double, N_S, 1>                       B;
    Eigen::Matrix<double, N_O, N_S>                     C;
    Eigen::Matrix<double, N_O, 1>                       D;
    Eigen::Matrix<double, N_S, N_S>                     Ad;
    Eigen::Matrix<double, N_S, N_C>                     Bd;
    Eigen::Matrix<double, N_O, N_S>                     Cd;
    Eigen::Matrix<double, N_O, N_C>                     Dd;

    // State Constraints 
    Eigen::Matrix<double, mpcWindow, N_S*mpcWindow>     G;
    Eigen::Matrix<double, mpcWindow, N_O>               S;


    // constraint matrices 
    Eigen::Matrix<double, 2*mpcWindow, 4>               Sbar;
    Eigen::Matrix<double, 2*mpcWindow, 1>               W0;
    Eigen::SparseMatrix<double>                         Gbar;

    // weight matrices - Updated 
    Eigen::Matrix<double, N_C, N_C>                     Q;    
    Eigen::Matrix<double, N_C, N_C>                     R;    
    Eigen::Matrix<double, N_C, N_C>                     RD;
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Qbar;    
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Rbar;
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> RbarD;        

    // Sx, Su, Su1, CAB, LL matrices  -- All updated other than Lu. Not used ? 
    Eigen::Matrix<double, N_C*mpcWindow, N_S>           Sx;
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Su;
    Eigen::Matrix<double, N_C*mpcWindow, N_C>           Su1;
    Eigen::Matrix<double, N_C*mpcWindow, N_C>           CAB;
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> LL;
    Eigen::Matrix<double, mpcWindow*N_C, N_C>           Lu;

    // state and the reference signal
    Eigen::Matrix<double, N_S, 1>                       X;
    Eigen::Matrix<double, N_O, mpcWindow>               ref;
    Eigen::Matrix<double, 1, N_O>                       K;

    // QP problem matrices and vectors
    Eigen::SparseMatrix<double>                         H;
    Eigen::Matrix<double, N_C*mpcWindow, N_C>           Fu;
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Fr;
    Eigen::Matrix<double, N_C*mpcWindow, N_S>           Fx;
    Eigen::Matrix<double, N_C*mpcWindow, 1>             f;
    Eigen::Matrix<double, 2*mpcWindow, 1>               lb;
    Eigen::Matrix<double, 2*mpcWindow, 1>               ub;    

    // timing variables
    Eigen::Matrix<double, 1, mpcWindow> t;
    double t0, dt;

    //  number of decisiona vars and constraints
    int n_variables;
    int n_constraints;  

    // output verbosity
    bool verbose;

    /**
     * Create a block diagonal matrix in the same style as MATLAB
     * 
     * @param a matrix block to be iterated
     * 
     * @param count number of times to iterate block a 
     */
    Eigen::MatrixXd blkdiag(const Eigen::MatrixXd &a, int count);

    /**
     * Create an evenly-spaced Eigen vector in the same style as MATLAB
     * 
     * @param start_in starting value of sequence
     * 
     * @param end_in end value of the sequence
     * 
     * @param num_in number of samples to generate 
     */
    Eigen::Matrix<double, 1, mpcWindow> linspace(double start_in, double end_in, int num_in);
};
