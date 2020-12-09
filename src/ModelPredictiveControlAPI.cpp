#include "./ModelPredictiveControlAPI.h"

ModelPredictiveControlAPI::ModelPredictiveControlAPI()
{
    S = G;
    W0 << 1, 1;

    // set arbitraty (to be filled) variables
    X << 0, 0, 0, 0;
    ref = Eigen::Matrix<double, N_O, mpcWindow>::Zero();
    lb = Eigen::Matrix<double, 2*mpcWindow, 1>::Ones() * Eigen::Infinity;
    t0 = 0;
    dt = 0;

    setSystemVars(Ad,Bd,Cd,Dd);
    // ModelPredictiveControlAPI::setSystemVars(Ad,Bd,Cd,Dd);

    setQ_R_RD(Q, R, RD, weightQ, weightR, weightRD);
    computeQbar_Rbar_RbarD(Qbar, Rbar, RbarD, Q, R, RD);
    computeSx_Su_Su1_CAB(Sx, Su, Su1, CAB, Ad, Cd, Bd);
    computeH(H, Rbar, Su, Qbar);
    setFVars(Fu, Fr, Fx, Qbar, Rbar, Sx, Su, Su1); 
    computeGbar_Sbar_W(Gbar, Sbar, W, G, S, W0);
    setF(f, Fu, Fr, Fx, X, ref);

    ub = W+Sbar*X;

    std::cout << "[MPC API]\t All QP matrices built successfully." << std::endl;

    n_variables = N_S*mpcWindow;
    n_constraints = 2*mpcWindow;
}

ModelPredictiveControlAPI::~ModelPredictiveControlAPI()
{
    printf("[MPC API]\t Destructing MPC API object...");
}

void ModelPredictiveControlAPI::setSystemVars(Eigen::Matrix<double,N_S,N_S> &Ad,
                                              Eigen::Matrix<double,N_S,N_C> &Bd,
                                              Eigen::Matrix<double,N_O,N_S> &Cd,
                                              Eigen::Matrix<double,N_O,N_C> &Dd)
{
    Ad <<   1.0001,     0.0050,     -0.0014,    -0.0001,
            0.0063,     1.0188,     -0.5056,    -0.0328,
            0.0003,     0.0008,      0.9783,     0.0035,
            0.1062,     0.3178,     -7.7149,     0.4474;

    Bd <<  -0.0001,    -0.0001,      0.0014,     0.0001,
           -0.0063,    -0.0188,      0.5027,     0.0314,
           -0.0003,    -0.0009,      0.0239,    -0.0015,
           -0.1062,    -0.3185,      8.4938,     0.5309;

    Cd = Eigen::Matrix<double,N_O,N_S>::Identity();

    Dd = Eigen::Matrix<double,N_O,N_C>::Zero();
    /*Dd <<   -0.0000,   -0.0025,      0.0016,    0.0001, 
            -0.0134,   -0.0053,      0.7694,    0.0482, 
            -0.0005,   -0.0002,      0.0266,   -0.0008, 
            -0.2266,   -0.0889,     12.4390,    0.8142; */

    std::cout << "[MPC API]\t System variables created."   << std::endl;
    // std::cout << "Ad:" << std::endl << Ad << std::endl << std::endl;
    // std::cout << "Bd:" << std::endl << Bd << std::endl << std::endl;
    // std::cout << "Cd:" << std::endl << Cd << std::endl << std::endl;
    // std::cout << "Dd:" << std::endl << Dd << std::endl << std::endl;
}


void ModelPredictiveControlAPI::setQ_R_RD(Eigen::Matrix<double,N_S,N_S>    &Q,
                                          Eigen::Matrix<double,N_C,N_C>    &R,
                                          Eigen::Matrix<double,1,1>        &RD,
                                          double wq,
                                          double wr,
                                          double wrd)
{
    Q  = Eigen::Matrix<double, N_S, N_S>::Identity() * wq;
    R  = Eigen::Matrix<double, N_C, N_C>::Identity() * wr;
    RD = Eigen::Matrix<double,1,1>::Identity() * wrd;

    std::cout << "[MPC API]\t Set Q, R, and RD matrices created."  << std::endl;
    // std::cout << "Q:"  << std::endl << Q  << std::endl << std::endl;
    // std::cout << "R:"  << std::endl << R  << std::endl << std::endl;
    // std::cout << "RD:" << std::endl << RD << std::endl << std::endl;
}


void ModelPredictiveControlAPI::computeQbar_Rbar_RbarD(Eigen::Matrix<double,N_S*mpcWindow,N_S*mpcWindow>   &Qbar,
                                                       Eigen::Matrix<double,N_C*mpcWindow,N_C*mpcWindow>   &Rbar,
                                                       Eigen::Matrix<double,1*mpcWindow,1*mpcWindow>       &RbarD,
                                                       Eigen::Matrix<double,N_S,N_S>                       Q,
                                                       Eigen::Matrix<double,N_C,N_C>                       R,
                                                       Eigen::Matrix<double,1,1>                           RD)
{
    Qbar = blkdiag(Q, mpcWindow);
    Rbar = blkdiag(R, mpcWindow);
    RbarD = blkdiag(RD, mpcWindow);

    std::cout << "[MPC API]\t Lifted weight matrices created."   << std::endl;
}


void ModelPredictiveControlAPI::computeSx_Su_Su1_CAB(Eigen::Matrix<double, N_S*mpcWindow, N_S>             &Sx,
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

    std::cout << "[MPC API]\t Sx Su, Su1, CAB created"   << std::endl;
}


void ModelPredictiveControlAPI::computeH(Eigen::SparseMatrix<double>                         &H,
                                         Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Rbar,
                                         Eigen::Matrix<double, N_S*mpcWindow, N_S*mpcWindow> Qbar,
                                         Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow> Su)
{
    Eigen::MatrixXd H_temp;
    H_temp = 2*(Rbar + Su.transpose() * Qbar * Su);

    H.resize(N_S*mpcWindow, N_S*mpcWindow); 

    for(int i=0; i<mpcWindow*N_S; i++)
    {
        for(int j=0; j<mpcWindow*N_S; j++)
        {
            H.insert(i,j) = H_temp(i,j);
        }
    }

    H.makeCompressed();

    std::cout << "[MPC API]\t Hessian H created." << std::endl;
}


void ModelPredictiveControlAPI::setFVars(Eigen::Matrix<double, N_C*mpcWindow, 1>               &Fu,
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

    std::cout << "[MPC API]\t Components of F created."   << std::endl;
}


void ModelPredictiveControlAPI::computeGbar_Sbar_W(Eigen::SparseMatrix<double>                     &Gbar,
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

    std::cout << "[MPC API]\t Constraints created."   << std::endl;
}


void ModelPredictiveControlAPI::setF(Eigen::Matrix<double, N_C*mpcWindow, 1>               &f,
                                     Eigen::Matrix<double, N_C*mpcWindow, 1>               Fu,
                                     Eigen::Matrix<double, N_C*mpcWindow, N_C*mpcWindow>   Fr, 
                                     Eigen::Matrix<double, N_C*mpcWindow, N_S>             Fx, 
                                     Eigen::Matrix<double, N_S, 1>                         X,
                                     Eigen::Matrix<double, N_S, mpcWindow>                 ref)
{
    Eigen::Map<Eigen::Matrix<double, N_C*mpcWindow, 1>> ref_vector(ref.data(), ref.size());
    f = Fx*X+Fr*ref_vector;
}


void ModelPredictiveControlAPI::updateRef(Eigen::Matrix<double, N_S, mpcWindow>            &ref, 
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

void ModelPredictiveControlAPI::linearizeABCD(double ts,
                                              double U_LQR,
                                              Eigen::Vector4d X)
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
