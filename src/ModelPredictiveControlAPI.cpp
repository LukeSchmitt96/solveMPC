#include "./ModelPredictiveControlAPI.h"

ModelPredictiveControlAPI::ModelPredictiveControlAPI(bool verbose_)
{

    std::cout << "[MPC API]\t MPC API object created." << std::endl;

    verbose = verbose_;
    solverFlag = true;

    K << -50.0, -150.0, 4000.0, 250.0;
    W0 = 255.0 * Eigen::Matrix<double, 2*mpcWindow, 1>::Ones();

    // set arbitraty (to be filled) variables
    X << 0.0, 0.0, 0.0, 0.0;
    U << 0.0;
    ref = Eigen::Matrix<double, N_O, mpcWindow>::Zero();
    lb = Eigen::Matrix<double, 2*mpcWindow, 1>::Ones() * Eigen::Infinity;
    t0 = 0.0;
    dt = 0.0;


    /** some of the following functions are commented out
     *      because we are not constraining the output
    */

    setSystemVars();
    setQ_R_RD();
    computeQbar_Rbar_RbarD();
    computeSx_Su_Su1_CAB();
    setLL();
    computeH();
    // setLu();
    setFVars(); 
    // computeS();
    // computeSbar();
    // computeGbar();
    updateRef(0.0);
    setF();

    ub = W0+Sbar*X;

    std::cout << "[MPC API]\t All QP matrices built successfully." << std::endl;

    n_variables = N_O*mpcWindow;
    n_constraints = 2*mpcWindow;


    // solver settings
    solver.settings()->setVerbosity(verbose);
    solver.settings()->setWarmStart(true);

    solver.data()->setNumberOfVariables(n_variables);
    solver.data()->setNumberOfConstraints(n_constraints);

    lb = Eigen::Matrix<double, 2*mpcWindow, 1>::Ones() * Eigen::Infinity;
    ub = Eigen::Matrix<double, 2*mpcWindow, 1>::Ones() * Eigen::Infinity;



    // this sets the qp problem as unconstrained!!
    // TEMPORARY REMOVE THIS IN FINAL CODE --------------------------------------------------------
    // TEMPORARY REMOVE THIS IN FINAL CODE --------------------------------------------------------
    
    lb = Eigen::Matrix<double, 2*mpcWindow, 1>::Ones() * Eigen::Infinity;
    ub = Eigen::Matrix<double, 2*mpcWindow, 1>::Ones() * -Eigen::Infinity;
    Eigen::Matrix<double, 2*mpcWindow, N_S> Gbar_temp = Eigen::Matrix<double, 2*mpcWindow, N_S>::Zero();

    Gbar.resize(2*mpcWindow, N_C*mpcWindow);
    for(int i = 0; i<Gbar_temp.rows(); i++)
    {
        for(int j=0; j<Gbar_temp.cols(); j++)
        {
            Gbar.insert(i,j) = Gbar_temp(i,j);
        }
    }
    Gbar.makeCompressed();

    std::cout << n_constraints << std::endl;
    std::cout << Gbar_temp.size() << std::endl;

    if(!solver.data()->setLinearConstraintsMatrix(Gbar))    {solverFlag = false; return;};
    if(!solver.data()->setLowerBound(lb))                   {solverFlag = false; return;};
    if(!solver.data()->setUpperBound(ub))                   {solverFlag = false; return;};

    // TEMPORARY REMOVE THIS IN FINAL CODE --------------------------------------------------------
    // TEMPORARY REMOVE THIS IN FINAL CODE --------------------------------------------------------




    if(!solver.data()->setHessianMatrix(H))                 {solverFlag = false; return;};
    if(!solver.data()->setGradient(f))                      {solverFlag = false; return;};


    // instantiate the solver
    if(!solver.initSolver()){solverFlag = false; return;};

}


ModelPredictiveControlAPI::~ModelPredictiveControlAPI()
{
    printf("[MPC API]\t Destructing MPC API object...\n");
}


void ModelPredictiveControlAPI::setVerbosity(bool verbose_)
{
    verbose = verbose_;
    std::cout << "[MPC API]\t Verbosity set to " << verbose << std::endl;
}


bool ModelPredictiveControlAPI::controllerStep()
{
    // update time vector
    t0 += dt;

    // set reference (for now just hold at 10.0)
    updateRef(10.0);

    // calculate new F based on ref, U and X
    setF();

    // update F
    if(!solver.updateGradient(f)) return false;

    // // update upper bound based on X
    // if(!solver.updateUpperBound(W0+Sbar*X)) return false;

    // solve the QP problem
    if(!solver.solve()) return false;

    // get the next controller input
    U += solver.getSolution().block<N_C, 1>(0, 0);

    return true;
}


void ModelPredictiveControlAPI::setSystemVars()
{
    Ad <<   1.00001778004125,       0.00505320517549851,    -0.0014579033084379,    -9.14586695703811e-05,
            0.00631068462774803,    1.01889078789663,       -0.525176404883372,     -0.0330113264471781,
            0.000300335964611163,   0.000898728268207066,    0.978168996256657,      0.00345977076140498,
            0.10666260164028,       0.319290150442394,      -7.75759238566293,       0.444855988055254;


    Bd <<   -1.77800412462243e-05,
            -0.00631068462774803,
            -0.000300335964611163,
            -0.10666260164028;

    Cd << 1, 0, 0, 0;

    Dd << 0; 

    if(verbose)
    {
        std::cout << "[MPC API]\t System variables created."   << std::endl;

        std::cout << "Ad rows: " << Ad.rows() << "\tAd cols: " << Ad.cols()  << std::endl;
        std::cout << "Ad:" << std::endl << Ad << std::endl << std::endl;

        std::cout << "Bd rows: " << Bd.rows() << "\tBd cols: " << Bd.cols()  << std::endl;
        std::cout << "Bd:" << std::endl << Bd << std::endl << std::endl;

        std::cout << "Cd rows: " << Cd.rows() << "\tCd cols: " << Cd.cols()  << std::endl;
        std::cout << "Cd:" << std::endl << Cd << std::endl << std::endl;

        std::cout << "Dd rows: " << Dd.rows() << "\tDd cols: " << Dd.cols()  << std::endl;
        std::cout << "Dd:" << std::endl << Dd << std::endl << std::endl;
    }
}

void ModelPredictiveControlAPI::setQ_R_RD()
{
    Q << 50.0; 
    R << 1.0 / 30;
    RD << 5.0;

    if(verbose)
    {
        std::cout << "[MPC API]\t Set Q, R, and RD matrices created."  << std::endl;

        std::cout << "Q rows: " << Q.rows() << "\tQ cols: " << Q.cols()  << std::endl;
        std::cout << "Q:"  << std::endl << Q  << std::endl << std::endl;

        std::cout << "R rows: " << R.rows() << "\tR cols: " << R.cols()  << std::endl;
        std::cout << "R:"  << std::endl << R  << std::endl << std::endl;

        std::cout << "RD rows: " << RD.rows() << "\tRD cols: " << RD.cols()  << std::endl;
        std::cout << "RD:" << std::endl << RD << std::endl << std::endl;
    }
}


void ModelPredictiveControlAPI::computeQbar_Rbar_RbarD()
{
    Qbar = blkdiag(Q, mpcWindow);
    Rbar = blkdiag(R, mpcWindow);
    RbarD = blkdiag(RD, mpcWindow);

    if(verbose)
    {
        std::cout << "[MPC API]\t Lifted weight matrices created."   << std::endl;

        std::cout << "Qbar rows: " << Qbar.rows() << "\tQbar cols: " << Qbar.cols()  << std::endl;
        std::cout << "Qbar:"  << std::endl << Qbar  << std::endl << std::endl;

        std::cout << "Rbar rows: " << Rbar.rows() << "\tRbar cols: " << Rbar.cols()  << std::endl;
        std::cout << "Rbar:"  << std::endl << Rbar  << std::endl << std::endl;

        std::cout << "RbarD rows: " << RbarD.rows() << "\tRbarD cols: " << RbarD.cols()  << std::endl;
        std::cout << "RbarD:" << std::endl << RbarD << std::endl << std::endl;
    }
}


void ModelPredictiveControlAPI::computeSx_Su_Su1_CAB()
{
    Sx  = Eigen::MatrixXd::Zero(Sx.rows(),  Sx.cols());
    Su  = Eigen::MatrixXd::Zero(Su.rows(),  Su.cols());
    Su1 = Eigen::MatrixXd::Zero(Su1.rows(), Su1.cols());
    CAB = Eigen::MatrixXd::Zero(CAB.rows(), CAB.cols());

    for(int i=0; i<mpcWindow; i++)
    {
        Sx.block<N_O, N_S>(i*N_O,0) = Cd*Ad.pow(i+1);
        CAB.block<N_C,N_O>(i*N_O,0) = Cd*Ad.pow(i)*Bd;
    }

    Eigen::VectorXd CAB_Vector(Eigen::Map<Eigen::VectorXd>(CAB.data(), CAB.cols()*CAB.rows()));

    // compute Su
    for(int i=0; i<Su.rows(); i++)
    {
        // TODO: IS THIS RIGHT??
        for(int j=0; j<=i; j++)
        {
            Su(i,j) = CAB_Vector.segment(0,i-j+1).sum();
        }
    }

    Su1 = Su.leftCols<N_O>();

    if(verbose)
    {
        std::cout << "[MPC API]\t Sx Su, Su1, CAB created"       << std::endl;

        std::cout << "Sx rows: " << Sx.rows() << "\tSx cols: " << Sx.cols()  << std::endl;
        std::cout << "Sx:"    << std::endl << Sx    << std::endl << std::endl;

        std::cout << "CAB rows: " << CAB.rows() << "\tCAB cols: " << CAB.cols()  << std::endl;
        std::cout << "CAB:"   << std::endl << CAB   << std::endl << std::endl;

        std::cout << "Su rows: " << Su.rows() << "\tSu cols: " << Su.cols()  << std::endl;
        std::cout << "Su:"    << std::endl << Su    << std::endl << std::endl;

        std::cout << "Su1 rows: " << Su1.rows() << "\tSu1 cols: " << Su1.cols()  << std::endl;
        std::cout << "Su1:"   << std::endl << Su1   << std::endl << std::endl;
    }
}


void ModelPredictiveControlAPI::computeH()
{
    Eigen::MatrixXd H_temp_1, H_temp_2;
    H_temp_1 = 2*(LL.transpose() * Rbar * LL + RbarD + Su.transpose() * Qbar * Su);
    H_temp_2 = (H_temp_1 + H_temp_1.transpose())/2.0;


    H.resize(N_O*mpcWindow, N_O*mpcWindow); 
    for(int i=0; i<H.rows(); i++)
    {
        for(int j=0; j<H.cols(); j++)
        {
            H.insert(i,j) = H_temp_2(i,j);
        }
    }

    H.makeCompressed();

    if(verbose)
    {
        std::cout << "[MPC API]\t Hessian H created." << std::endl;
        std::cout << "H rows: " << H_temp_1.rows() << "\tH cols: " << H_temp_1.cols()  << std::endl;
        std::cout << "H:" << std::endl << H_temp_1 << std::endl  << std::endl;
    }
    
}

void ModelPredictiveControlAPI::setLu()
{
    for(int i=0; i<mpcWindow; i++)
    {
        // TODO: Ask Z about this
        Lu.block<N_C,N_C>(i,0) = (mpcWindow-i+2)*Eigen::Matrix<double,N_C,N_C>::Identity();
    }

    if(verbose)
    {
        std::cout << "[MPC API]\t Lu created." << std::endl;
        std::cout << "Lu rows: " << Lu.rows() << "\tLu cols: " << Lu.cols()  << std::endl;
        std::cout << "Lu:" << std::endl << Lu << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI:: setLL()
{
    // Eigen::Matrix<double, N_S, N_S> tempEye = Eigen::Matrix<double, N_S, N_S>::Identity();

    // std::cout << tempEye << std::endl;

    // Eigen::Matrix<double, N_S*mpcWindow, N_S*mpcWindow > tempRep = tempEye.replicate(mpcWindow,mpcWindow);


    // std::cout << tempRep << std::endl;
    
    LL = Eigen::MatrixXd(Eigen::Matrix<double, N_O*mpcWindow, N_O*mpcWindow>::Ones().triangularView<Eigen::Lower>());

    if(verbose)
    {
        std::cout << "[MPC API]\t LL created."   << std::endl;
        std::cout << "LL rows: " << LL.rows() << "\tLL cols: " << LL.cols()  << std::endl;
        std::cout << "LL:" << std::endl << LL << std::endl  << std::endl;
    }
    
} 


void ModelPredictiveControlAPI::setFVars()
{
    Fu = 2*((LL.transpose() * Rbar.transpose()).diagonal().transpose() + Su1.transpose() * Qbar * Su).transpose();
    // Fu = 2*((LL.transpose() * Rbar.transpose()).transpose() + Su1.transpose() * Qbar * Su).transpose();
    Fr = -2*(Qbar*Su).transpose();
    Fx = 2*(Sx.transpose() * Qbar * Su).transpose();

    if(verbose)
    {
        std::cout << "[MPC API]\t Components of F created."   << std::endl;

        std::cout << "Fu rows: " << Fu.rows() << "\tFu cols: " << Fu.cols()  << std::endl;
        std::cout << "Fu:" << std::endl << Fu << std::endl  << std::endl;

        std::cout << "Fr rows: " << Fr.rows() << "\tFr cols: " << Fr.cols()  << std::endl;
        std::cout << "Fr:" << std::endl << Fr << std::endl  << std::endl;

        std::cout << "Fx rows: " << Fx.rows() << "\tFx cols: " << Fx.cols()  << std::endl;
        std::cout << "Fx:" << std::endl << Fx << std::endl  << std::endl;
    }
    
}


// void ModelPredictiveControlAPI::computeS()
// {
//     for(int i=0; i<mpcWindow; i++)
//     {
//         S.block<1,N_S>(i,0) = -K*Ad.pow(i+1);
//     }

//     if(verbose)
//     {
//         std::cout << "[MPC API]\t S created."   << std::endl;
//         std::cout << "S:" << std::endl << S << std::endl  << std::endl;
//     }
// }


// void ModelPredictiveControlAPI::computeSbar()
// {
//     Sbar << S, -S;

//     if(verbose)
//     {
//         std::cout << "[MPC API]\t Sbar created."   << std::endl;
//         std::cout << "Sbar:" << std::endl << Sbar << std::endl  << std::endl;
//     }
// }


// void ModelPredictiveControlAPI::computeG()
// {
//     Eigen::Matrix<double, N_S, N_C> AiB;
//     AiB = Eigen::Matrix<double, N_S, N_C>::Zero();

//     Eigen::Matrix<double, N_S*mpcWindow, N_C> AB;
    
//     for(int i; i<mpcWindow; i++)
//     {
//         AiB += Ad.pow(i) * Bd;
//         AB.block<N_S, N_S>(i*N_S, 0) = AiB;
//     }


//     for(int i=0; i<mpcWindow; i++)
//     {
//         for(int j=0; j<=i; j++)
//         {
//             G.block<1,N_S>(i*N_S, 4*j-3) = -K*(Eigen::Matrix<double, N_S, N_S>::Identity() - AB.block<N_S, N_S>(N_S*(i-j), 0)); 
//         }
//     }

//     if(verbose)
//     {
//         std::cout << "[MPC API]\t G created."   << std::endl;
//         std::cout << "G:" << std::endl << G << std::endl  << std::endl;
//     }
// }


// void ModelPredictiveControlAPI::computeGbar()
// {
//     Eigen::Matrix<double, 2*mpcWindow, N_C*mpcWindow> Gbar_temp;
//     Gbar_temp << G, -G;
//     Gbar.resize(2*mpcWindow, N_C*mpcWindow);

//     for(int i=0; i<2*mpcWindow; i++)
//     {
//         for(int j=0; j<N_C*mpcWindow; j++)
//         {
//             Gbar.insert(i,j) = Gbar_temp(i,j);
//         }
//     }

//     Gbar.makeCompressed();

//     if(verbose)
//     {
//         std::cout << "[MPC API]\t Gbar created."   << std::endl;
//         std::cout << "Gbar:" << std::endl << Gbar_temp << std::endl  << std::endl;
//     }
// }


void ModelPredictiveControlAPI::setF()
{

    f = Fx*X + Fu*U + Fr*ref.transpose();

    if(verbose)
    {
        std::cout << "[MPC API]\tf created."   << std::endl;

        std::cout << "f rows: " << f.rows() << "\tf cols: " << f.cols()  << std::endl;
        std::cout << "f:" << std::endl << f << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI::updateRef(double pos_ref)
{
    // t = linspace(t0, t0 + dt/1000.0*(mpcWindow-1), mpcWindow);
    // ref = Eigen::Matrix<double, N_C, mpcWindow>::Zero();

    // Eigen::Matrix<double, 1, mpcWindow> ref_position;
    // Eigen::Matrix<double, 1, mpcWindow> temp;
    
    // temp = t.unaryExpr([](const double x) 
    //     { 
    //         return fmod(x/(Ts/2),2);
    //     }
    // );

    // ref_position = temp.unaryExpr([](const double x)
    //     {
    //         return x <= 1.0 ? 1.0 : 0.0;
    //     }
    // );

    ref.block<1,mpcWindow>(0,0) = pos_ref * Eigen::Matrix<double, N_C, mpcWindow>::Ones();

    if(verbose)
    {
        std::cout << "[MPC API]\tref: " << ref << std::endl;
    }
}


void ModelPredictiveControlAPI::linearizeABCD(double ts,
                                              double U_LQR)
{
    double sinx3 = sin(X(2));
    double cosx3 = cos(X(2));
    
    Eigen::Vector4d va4;
    va4(0) = 0.;
    va4(1) = cosx3/(0.80128205128*cosx3*cosx3 - 0.2086653922);
    va4(2) = (0.0989561664*cosx3 - 0.0001557504*X(3)*X(3)*cosx3*cosx3 + sinx3*(0.0001557504*sinx3*X(3)*X(3) 
            - 0.00026676593*U_LQR + 0.0001248*X(1)))/(0.0001557504*cosx3*cosx3 - 0.00059808672) 
            - (5746175536355785.*cosx3*sinx3*(0.0989561664*sinx3 - cosx3*(0.0001557504*sinx3*X(3)*X(3) 
            - 0.00026676593*U_LQR + 0.0001248*X(1))))/(18446744073709551616.*pow(0.0001557504*cosx3*cosx3 - 0.00059808672,2));
    va4(3) = (X(3)*cosx3*sinx3)/(2*cosx3*cosx3 - 0.52082881893);
    
    Eigen::Vector4d va2;
    va2(0) = 0.;
    va2(1) = 136358332192861.0/(18446744073709551616.0*(0.0001557504*cosx3*cosx3 - 0.00059808672));
    va2(2) = - (3125.*((1911.*cosx3)/15625. + (3408958304821525.*(0.0989561664*cosx3 
                - 0.0001557504*X(3)*X(3)*cosx3*cosx3 + sinx3*(0.0001557504*sinx3*X(3)*X(3) - 0.00026676593*U_LQR 
                + (39.0*X(1))/312500.)))/(7.1827194e+14*cosx3*cosx3 - 2.7581882e+15) 
                + (19588472815622334031693326272125.*cosx3*sinx3*(0.0989561664*sinx3 
                - cosx3*((1521.0*sinx3*pow(X(3),2))/9765625. - 0.00026676593*U_LQR 
                + 0.0001248*X(1))))/(85070591730234615865843651857942052864.*pow(0.0001557504*cosx3*cosx3 
                - 0.00059808672,2))))/(39.*cosx3) - (3125.*sinx3*((1911.*sinx3)/15625. 
                + (3408958304821525.*(0.0989561664*sinx3 - cosx3*((1521.0*sinx3*pow(X(3),2))/9765625.0 
                - 0.00026676593*U_LQR + 0.0001248*X(1))))/(7.1827194e+14*cosx3*cosx3 - 2.7581882e+15)))/(39.*cosx3*cosx3);
    va2(3) = (X(3)*sinx3)/(0.11846153945*cosx3*cosx3 - 0.03084909163);

    double va4_u = -cosx3/(1.71277846634*cosx3*cosx3 - 0.44603219682);
    double va2_u = -1/(0.101449188*cosx3*cosx3 - 0.02641883125);

}

// void ModelPredictiveControlAPI::c2d(double                                                 ts,
//                                     Eigen::Matrix<double, N_S, N_S>                        &Ad,
//                                     Eigen::Matrix<double, N_S, 1>                          &Bd,
//                                     Eigen::Matrix<double, N_O, N_S>                        &Cd,
//                                     Eigen::Matrix<double, N_O, 1>                          &Dd,
//                                     Eigen::Matrix<double, N_S, N_S>                        A,
//                                     Eigen::Matrix<double, N_S, 1>                          B,
//                                     Eigen::Matrix<double, N_O, N_S>                        C,
//                                     Eigen::Matrix<double, N_O, 1>                          D)
// {
//     // build an exponential matrix
//     Eigen::Matrix<double, N_S, N_S + 1> em_upper;
//     em_upper << A, B;

//     // Need to stack zeros under the a and b matrices
//     Eigen::Matrix<double, N_S, N_S + 1> em_lower;
//     em_lower << Eigen::Matrix4d::Zero((B.rows(), A.rows())), Eigen::Matrix4d::Zero((B.rows(), B.rows()));

//     Eigen::MatrixXd em;
//     em << em_upper, 
//           em_lower;
    
//     // do zoh
//     Eigen::MatrixXd ms_temp;
//     ms_temp << linalg.expm(dt * em);

//     // Dispose of the lower rows
//     Eigen::MatrixXd ms;
//     ms = ms_temp[:A.shape[0], :];

//     Ad = ms[:, 0:A.shape[1]];
//     Bd = ms[:, A.shape[1]:];

//     Cd = C;
//     Dd = D;
// }

Eigen::MatrixXd ModelPredictiveControlAPI::blkdiag(const Eigen::MatrixXd& a, int count)
{
    Eigen::MatrixXd bdm = Eigen::MatrixXd::Zero(a.rows() * count, a.cols() * count);
    for (int i = 0; i < count; ++i)
    {
        bdm.block(i * a.rows(), i * a.cols(), a.rows(), a.cols()) = a;
    }

    return bdm;
}

Eigen::Matrix<double, 1, mpcWindow> ModelPredictiveControlAPI::linspace(double start_in, double end_in, int num_in)
{

    Eigen::Matrix<double, 1, mpcWindow> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num; ++i)
    {
        linspaced(0,i) = start + delta * i;
    }
    //linspaced(0,num-1) = end;   // I want to ensure that start and end
                                // are exactly the same as the input
        
    return linspaced;
}
