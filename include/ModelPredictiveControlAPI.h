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
const int N_O = 4;  // number of outputs

// set weights
const double weightQ =  1.0;    // state weight
const double weightR =  1e-3;   // control weight
const double weightRD = 0.0;    // change in control weight

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
     * 
     * @param Ad pointer to system matrix A
     * 
     * @param Bd pointer to system matrix B
     * 
     * @param Cd pointer to system matrix C
     * 
     * @param Dd pointer to system matrix D
     */
    void setSystemVars(Eigen::Matrix<double,N_S,N_S>                                &,
                       Eigen::Matrix<double,N_S,N_C>                                &,
                       Eigen::Matrix<double,N_O,N_S>                                &,
                       Eigen::Matrix<double,N_O,N_C>                                &);

    /**
     * Set the state, control, and change in control weight matrices based
     *      on weight inputs
     * 
     * @param Q     state weight matrix pointer 
     * 
     * @param R     control weight matrix pointer   
     * 
     * @param RD    change in control weight matrix pointer 
     * 
     * @param wq    state weight value  
     * 
     * @param wr    control weight value    
     * 
     * @param wrd   change in control weight value
     */
    void setQ_R_RD(Eigen::Matrix<double,N_S,N_S>                                    &,
                   Eigen::Matrix<double,N_C,N_C>                                    &,
                   Eigen::Matrix<double,1,1>                                        &,
                   double                                                            ,
                   double                                                            ,
                   double                                                            );

    /**
     * Set the lifted state, control, and change in control weight matrices
     * 
     * @param Qbar      lifted state weight matrix
     * 
     * @param Rbar      lifted control weight matrix
     * 
     * @param RbarD     lifted change in control weight matrix
     * 
     * @param Q         state weight matrix pointer 
     * 
     * @param R         control weight matrix pointer   
     * 
     * @param RD        change in control weight matrix pointer 
     */
    void computeQbar_Rbar_RbarD(Eigen::Matrix<double,N_S*mpcWindow,N_S*mpcWindow>   &,
                                Eigen::Matrix<double,N_C*mpcWindow,N_C*mpcWindow>   &,
                                Eigen::Matrix<double,1*mpcWindow,1*mpcWindow>       &,
                                Eigen::Matrix<double,N_S,N_S>                       ,
                                Eigen::Matrix<double,N_C,N_C>                       ,
                                Eigen::Matrix<double,1,1>                           );

    /**
     * Complete lifted dynamics transform matrices
     * 
     * @param Sx        pointer to lifted control to state transform
     * 
     * @param Su        pointer to lifted control to output transform
     * 
     * @param Su1       pointer to first column of Su
     * 
     * @param CAB       pointer to state propogation matrix
     * 
     * @param Ad        system matrix A
     * 
     * @param Bd        system matrix B
     * 
     * @param Cd        system matrix C
     */
    void computeSx_Su_Su1_CAB(Eigen::Matrix<double, N_S*mpcWindow, N_S>             &,
                              Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>   &,
                              Eigen::Matrix<double, N_C*mpcWindow, 1>               &,
                              Eigen::Matrix<double, N_S*mpcWindow, N_S>             &,
                              Eigen::Matrix<double, N_S, N_S>                       ,
                              Eigen::Matrix<double, N_S, N_C>                       ,
                              Eigen::Matrix<double, N_O, N_S>                       );       

    /**
     * Compute hessian for the QP problem
     * 
     * @param H         pointer to QP hessian
     * 
     * @param Rbar      lifted control weight matrix
     * 
     * @param Qbar      lifted state weight matrix
     * 
     * @param Su        lifted control to output transform
     */    
    void computeH(Eigen::SparseMatrix<double>                                       &,
                  Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>               ,
                  Eigen::Matrix<double, N_S*mpcWindow, N_S*mpcWindow>               ,
                  Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>               );

    /**
     * Compute components for the qp gradient
     * 
     * @param Fu        pointer to gradient component for control
     * 
     * @param Fr        pointer to gradient component for reference signal
     * 
     * @param Fx        pointer to gradient component for state
     * 
     * @param Qbar      lifted state weight matrix
     * 
     * @param Rbar      lifted control weight matrix   
     * 
     * @param Sx        pointer to lifted control to state transform
     * 
     * @param Su        pointer to lifted control to output transform
     * 
     * @param Su1       pointer to first column of Su    
     */
    void setFVars(Eigen::Matrix<double, N_C*mpcWindow, 1>                           &,
                  Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>               &, 
                  Eigen::Matrix<double, N_C*mpcWindow, N_S>                         &, 
                  Eigen::Matrix<double, N_S*mpcWindow, N_S*mpcWindow>               ,
                  Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>               , 
                  Eigen::Matrix<double, N_S*mpcWindow, N_S>                         ,
                  Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>               , 
                  Eigen::Matrix<double, N_C*mpcWindow, 1>                           );    

    /**
     * Compute constraints components for the qp problem
     * 
     * @param Gbar      pointer to linear constraints matrix
     * 
     * @param Sbar      pointer to gradient component for constraints
     * 
     * @param W         pointer to bias component for constraints
     * 
     * @param G         lifted constraints on control
     * 
     * @param S         lifted constraints on control
     * 
     * @param W0        components of bias term for constraints
     */
    void computeGbar_Sbar_W(Eigen::SparseMatrix<double>                             &,
                            Eigen::Matrix<double, 2*mpcWindow, N_C>                 &,
                            Eigen::Matrix<double, 2*mpcWindow, 1>                   &,
                            Eigen::Matrix<double, 2, N_C>                           ,
                            Eigen::Matrix<double, 2, N_C>                           ,
                            Eigen::Matrix<double, 2, 1>                             );          

    /**
     * compute gradient for the qp problem
     * 
     * @param f         pointer to qp problem gradient
     * 
     * @param Fu        gradient component for control
     * 
     * @param Fr        gradient component for reference
     * 
     * @param Fx        gradient component for state
     * 
     * @param X         current state
     * 
     * @param ref       current reference vectors
     * 
     * @param t         time vector corresponding to reference
     */    
    void setF(Eigen::Matrix<double, N_C*mpcWindow, 1>                               &,
              Eigen::Matrix<double, N_C*mpcWindow, 1>                               ,
              Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>                   ,
              Eigen::Matrix<double, N_C*mpcWindow, N_S>                             ,
              Eigen::Matrix<double, N_S, 1>                                         ,
              Eigen::Matrix<double, N_S, mpcWindow>                                 );

    /**
     * update the reference given the time
     * 
     * @param ref       pointer to the current referenece vectors
     * 
     * @param t         time vector corresponding to reference
     */
    void updateRef(Eigen::Matrix<double, N_S, mpcWindow>                            &, 
                   Eigen::Matrix<double, 1, mpcWindow>                              );

    /**
     * linearizes nonlinear model into a linear model about X
     * 
     * @param ts        sampling time
     * 
     * @param U_LQR     control from LQR
     * 
     * @param X         state to linearize about
     */
    void linearizeABCD(double ts,
                       double U,
                       Eigen::Vector4d X);

    /**
     * converts a continuous time system to discrete time
     * 
     * @param ts        sampling time
     * 
     * @param Ad        pointer to Ad
     * 
     * @param Bd        pointer to Bd
     * 
     * @param Cd        pointer to Cd
     * 
     * @param Dd        pointer to Dd
     * 
     * @param A         system matrix A
     * 
     * @param B         system matrix B
     * 
     * @param C         system matrix C
     * 
     * @param D         system matrix D
     */
    void c2d(double                                                                 ,
             Eigen::Matrix<double, N_S, N_S>                                        &,
             Eigen::Matrix<double, N_S, 1>                                          &,
             Eigen::Matrix<double, N_O, N_S>                                        &,
             Eigen::Matrix<double, N_O, 1>                                          &,
             Eigen::Matrix<double, N_S, N_S>                                        ,
             Eigen::Matrix<double, N_S, 1>                                          ,
             Eigen::Matrix<double, N_O, N_S>                                        ,
             Eigen::Matrix<double, N_O, 1>                                          );

    // system matrices
    Eigen::Matrix<double, N_S, N_S>                     A;
    Eigen::Matrix<double, N_S, 1>                       B;
    Eigen::Matrix<double, N_O, N_S>                     C;
    Eigen::Matrix<double, N_O, 1>                       D;
    Eigen::Matrix<double, N_S, N_S>                     Ad;
    Eigen::Matrix<double, N_S, N_C>                     Bd;
    Eigen::Matrix<double, N_O, N_S>                     Cd;
    Eigen::Matrix<double, N_O, N_C>                     Dd;
    Eigen::Matrix<double, 2, N_C>                       G;
    Eigen::Matrix<double, 2, N_C>                       S; 

    // constraint matrices
    Eigen::Matrix<double, 2*mpcWindow, 4>               Sbar;
    Eigen::Matrix<double, 2, 1>                         W0; 
    Eigen::Matrix<double, 2*mpcWindow, 1>               W;
    Eigen::SparseMatrix<double>                         Gbar;    

    // weight matrices
    Eigen::Matrix<double, N_S, N_S>                     Q;    
    Eigen::Matrix<double, N_C, N_C>                     R;    
    Eigen::Matrix<double, 1,   1>                       RD;
    Eigen::Matrix<double, N_S*mpcWindow, N_S*mpcWindow> Qbar;    
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Rbar;
    Eigen::Matrix<double, 1*mpcWindow,   1*mpcWindow>   RbarD;        

    // Sx, Su, Su1, CAB matrices
    Eigen::Matrix<double, N_S*mpcWindow, N_S>           Sx;
    Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Su;
    Eigen::Matrix<double, N_C*mpcWindow, 1>             Su1;
    Eigen::Matrix<double, N_S*mpcWindow, N_S>           CAB;

    // state and the reference signal
    Eigen::Matrix<double, N_S, 1>                       X;
    Eigen::Matrix<double, N_O, mpcWindow>               ref;

    // QP problem matrices and vectors
    Eigen::SparseMatrix<double>                         H;
    Eigen::Matrix<double, N_C*mpcWindow, 1>             Fu;
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

    /**
     * Create a block diagonal matrix in the same style as MATLAB
     * 
     * @param a matrix block to be iterated
     * 
     * @param count number of times to iterate block a 
     */
    Eigen::MatrixXd blkdiag(const Eigen::MatrixXd &, int);

    /**
     * Create an evenly-spaced Eigen vector in the same style as MATLAB
     * 
     * @param start_in starting value of sequence
     * 
     * @param end_in end value of the sequence
     * 
     * @param num_in number of samples to generate 
     */
    Eigen::Matrix<double, 1, mpcWindow> linspace(double, double, int);
};
