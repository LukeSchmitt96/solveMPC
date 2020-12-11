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

const int mpcWindow = 10;           // MPC preview window

// set matrix dimensions
const int N_S = 4;                  // number of states
const int N_C = 1;                  // number of controls
const int N_O = 1;                  // number of outputs

// timing parameters
const double Ts = 100;              // square wave period

// MPC API Class
class ModelPredictiveControlAPI 
{

public:
    /** constructor
     * 
     * @param verbose_              // verbosity of the API
     */
    ModelPredictiveControlAPI(bool);

    // destructor
    ~ModelPredictiveControlAPI();
    
    /**
     * Sets the verbosity of the MPC API
     * 
     * @param verbose_          // verbosity of the API
     */
    void setVerbosity(bool);

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
     * 
     * @param pos_ref           position reference in meters
     */
    void updateRef(double);

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

    /**
     * performs a controller step by solving the qp problem
     *
     * @returns solver success flag
     */
    bool controllerStep();

    // OSQPEigen Solver
    OsqpEigen::Solver solver;
    bool solverFlag;

    // system matrices
    Eigen::Matrix<double, N_S, N_S>                     A;      // system A matrix
    Eigen::Matrix<double, N_S, 1>                       B;      // system B matrix
    Eigen::Matrix<double, N_O, N_S>                     C;      // system C matrix
    Eigen::Matrix<double, N_O, 1>                       D;      // system D matrix
    Eigen::Matrix<double, N_S, N_S>                     Ad;     // system Ad matrix
    Eigen::Matrix<double, N_S, N_C>                     Bd;     // system Bd matrix
    Eigen::Matrix<double, N_O, N_S>                     Cd;     // system Cd matrix
    Eigen::Matrix<double, N_O, N_C>                     Dd;     // system Dd matrix

    // State Constraints 
    Eigen::Matrix<double, mpcWindow, N_S*mpcWindow>     G;      // QP linear constraints matrix component
    Eigen::Matrix<double, mpcWindow, N_S>               S;      // QP S
    Eigen::Matrix<double, 2*mpcWindow, 4>               Sbar;   // QP Sbar    
    Eigen::Matrix<double, 2*mpcWindow, 1>               W0;     // QP W0
    Eigen::SparseMatrix<double>                         Gbar;   // QP linear constraints matrix

    // weight matrices 
    Eigen::Matrix<double, N_O, N_O>                     Q;      // State cost
    Eigen::Matrix<double, N_O, N_O>                     R;      // Control cost
    Eigen::Matrix<double, N_O, N_O>                     RD;     // Change in control cost
    Eigen::Matrix<double, N_O*mpcWindow, N_O*mpcWindow> Qbar;   // Lifted state cost
    Eigen::Matrix<double, N_O*mpcWindow, N_O*mpcWindow> Rbar;   // Lifted control cost
    Eigen::Matrix<double, N_O*mpcWindow, N_O*mpcWindow> RbarD;  // Lifted change in control cost

    // Sx, Su, Su1, CAB, LL matrices
    Eigen::Matrix<double, N_C*mpcWindow, N_S>           Sx;
    Eigen::Matrix<double, N_O*mpcWindow, N_O*mpcWindow> Su;
    Eigen::Matrix<double, N_O*mpcWindow, N_O>           Su1;
    Eigen::Matrix<double, N_C*mpcWindow, N_C>           CAB;
    Eigen::Matrix<double, N_O*mpcWindow, N_O*mpcWindow> LL;
    Eigen::Matrix<double, mpcWindow*N_C, N_C>           Lu;

    // state and the reference signal
    Eigen::Matrix<double, N_S, 1>                       X;      // State vector
    Eigen::Matrix<double, N_O, N_O>                     U;      // Control vector
    Eigen::Matrix<double, N_O, mpcWindow>               ref;    // Regerence vector
    Eigen::Matrix<double, 1, N_S>                       K;      // LQR feedback law

    // QP problem matrices and vectors
    Eigen::SparseMatrix<double>                         H;      // QP Hessian
    Eigen::Matrix<double, N_C*mpcWindow, N_C>           Fu;     // Control component of QP gradient
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Fr;     // Reference component of QP gradient
    Eigen::Matrix<double, N_C*mpcWindow, N_S>           Fx;     // State component of QP gradient
    Eigen::Matrix<double, N_C*mpcWindow, 1>             f;      // QP gradient
    Eigen::Matrix<double, 2*mpcWindow, 1>               lb;     // QP lower bound
    Eigen::Matrix<double, 2*mpcWindow, 1>               ub;     // QP upper bound

    // timing variables
    Eigen::Matrix<double, 1, mpcWindow> t;
    double t0, dt;

    int n_variables;        // number of decision variables
    int n_constraints;      // number of constraints

    bool verbose;           // output verbosity

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
