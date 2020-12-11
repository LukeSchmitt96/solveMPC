#include "./ModelPredictiveControlAPI.h"

ModelPredictiveControlAPI::ModelPredictiveControlAPI()
{
    verbose = false;

    K << -50, -150, 4000, 250;
    W0 = 255.0 * Eigen::Matrix<double, 2*mpcWindow, 1>::Ones();

    // set arbitraty (to be filled) variables
    X << 0, 0, 0, 0;
    ref = Eigen::Matrix<double, N_O, mpcWindow>::Zero();
    lb = Eigen::Matrix<double, 2*mpcWindow, 1>::Ones() * Eigen::Infinity;
    t0 = 0;
    dt = 0;

    setSystemVars();
    setQ_R_RD();
    computeQbar_Rbar_RbarD();
    computeSx_Su_Su1_CAB();
    computeH();
    setLu();
    setFVars(); 
    computeS();
    computeSbar();
    computeGbar();
    setF();
    setLL();

    ub = W0+Sbar*X;

    std::cout << "[MPC API]\t All QP matrices built successfully." << std::endl;

    n_variables = N_S*mpcWindow;
    n_constraints = 2*mpcWindow;
}

ModelPredictiveControlAPI::~ModelPredictiveControlAPI()
{
    printf("[MPC API]\t Destructing MPC API object...");
}

void ModelPredictiveControlAPI::setSystemVars()
{
    Ad <<   1.0001,     0.0102,     -0.0048,    -0.0003,
            0.0091,     1.0272,     -0.7927,    -0.0501,
            0.0010,     0.0029,      0.9295,     0.0049,
            0.1538,     0.4609,     -11.2071,    0.1605;

    Bd <<  -0.0001,   
           -0.0091,  
           -0.0010,   
           -0.1538;

    Cd << 1, 0, 0, 0;

    Dd = << 1; 

    if(verbose)
    {
        std::cout << "[MPC API]\t System variables created."   << std::endl;
        std::cout << "Ad:" << std::endl << Ad << std::endl << std::endl;
        std::cout << "Bd:" << std::endl << Bd << std::endl << std::endl;
        std::cout << "Cd:" << std::endl << Cd << std::endl << std::endl;
        std::cout << "Dd:" << std::endl << Dd << std::endl << std::endl;
    }
}

void ModelPredictiveControlAPI::setQ_R_RD()
{
    Q << 50; 
    R << 1 / 30;
    RD << 5

    if(verbose)
    {
        std::cout << "[MPC API]\t Set Q, R, and RD matrices created."  << std::endl;
        std::cout << "Q:"  << std::endl << Q  << std::endl << std::endl;
        std::cout << "R:"  << std::endl << R  << std::endl << std::endl;
        std::cout << "RD:" << std::endl << RD << std::endl << std::endl;
    }
}


void ModelPredictiveControlAPI::computeQbar_Rbar_RbarD()
{
    Qbar = blkdiag(Q, mpcWindow);
    Rbar = blkdiag(R, mpcWindow);
    RbarD = blkdiag(RD, mpcWindow);

    std::cout << "[MPC API]\t Lifted weight matrices created."   << std::endl;
    if(verbose)
    {
        std::cout << "Qbar:"  << std::endl << Qbar  << std::endl << std::endl;
        std::cout << "Rbar:"  << std::endl << Rbar  << std::endl << std::endl;
        std::cout << "RbarD:" << std::endl << RbarD << std::endl << std::endl;
    }
}


void ModelPredictiveControlAPI::computeSx_Su_Su1_CAB()
{
    Sx  = Eigen::MatrixXd::Zero(Sx.rows(),  Sx.cols());
    Su  = Eigen::MatrixXd::Zero(Su.rows(),  Su.cols());
    Su1 = Eigen::MatrixXd::Zero(Su1.rows(), Su1.cols());
    CAB = Eigen::MatrixXd::Zero(CAB.rows(), CAB.cols());

    // TODO: compute Sx and CAB --- Sx = [Sx;Cd*Ad^ii];
    // CAB = [CAB;CAiB]; 
    /*CAiB = CAiB + Cd * Ad ^ (ii - 1) * Bd;
    AiB = AiB + Ad ^ (ii - 1) * Bd;*/
    for(int i=0; i<mpcWindow; i++)
    {
        Sx.block<N_C, N_S>(i*N_O,0) = Cd*Ad.pow(i+1);
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

    Su1 = Su.leftCols<1>();

    if(verbose)
    {
        std::cout << "[MPC API]\t Sx Su, Su1, CAB created"       << std::endl;
        std::cout << "Sx:"    << std::endl << Sx    << std::endl << std::endl;
        std::cout << "CAB:"   << std::endl << CAB   << std::endl << std::endl;
        std::cout << "Su:"    << std::endl << Su    << std::endl << std::endl;
        std::cout << "Su1:"   << std::endl << Su1  << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI::computeH()
{
    Eigen::MatrixXd H_temp;
    H_temp = 2*(LL.transpose() * Rbar * LL + RbarD + Su.transpose() * Qbar * Su);

    H.resize(N_S*mpcWindow, N_S*mpcWindow); 

    for(int i=0; i<mpcWindow*N_S; i++)
    {
        for(int j=0; j<mpcWindow*N_S; j++)
        {
            H.insert(i,j) = H_temp(i,j);
        }
    }

    H.makeCompressed();

    if(verbose)
    {
        std::cout << "[MPC API]\t Hessian H created." << std::endl;
        std::cout << "H:" << std::endl << H_temp << std::endl  << std::endl;
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
        std::cout << "Lu:" << std::endl << Lu << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI:: setLL()
{
    Eigen::Matrix<double, N_S, N_S> tempEye = Eigen::Matrix<double, N_S, N_S>::Identity();
    Eigen::Matrix<double, N_S*mpcWindow, N_S*mpcWindow > tempRep = tempEye.replicate(mpcWindow,mpcWindow);
    
    std::cout << "[MPC API]\t LL created."   << std::endl;
    LL = Eigen::MatrixXd(tempRep.triangularView<Eigen::Lower>());
} 


void ModelPredictiveControlAPI::setFVars()
{
    Fu = 2*((Rbar*Lu).transpose() + Su1.transpose() * Qbar * Su).transpose();
    Fr = -2*(Qbar*Su).transpose();
    Fx = 2*(Sx.transpose() * Qbar * Su).transpose();

    if(verbose)
    {
        std::cout << "[MPC API]\t Components of F created."   << std::endl;
        std::cout << "Fu:" << std::endl << Fu << std::endl  << std::endl;
        std::cout << "Fr:" << std::endl << Fr << std::endl  << std::endl;
        std::cout << "Fx:" << std::endl << Fx << std::endl  << std::endl;
    }
    
}


void ModelPredictiveControlAPI::computeS()
{
    for(int i=0; i<mpcWindow; i++)
    {
        S.block<1,N_C>(i,0) = -K*Ad.pow(i+1);
    }

    if(verbose)
    {
        std::cout << "[MPC API]\t S created."   << std::endl;
        std::cout << "S:" << std::endl << S << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI::computeSbar()
{
    Sbar << S, -S;

    if(verbose)
    {
        std::cout << "[MPC API]\t Sbar created."   << std::endl;
        std::cout << "Sbar:" << std::endl << Sbar << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI::computeG()
{
    Eigen::Matrix<double, N_S, N_S> AiB;
    AiB = Eigen::Matrix<double, N_S, N_S>::Zero();

    Eigen::Matrix<double, N_S*mpcWindow, N_S> AB;
    
    for(int i; i<mpcWindow; i++)
    {
        AiB += Ad.pow(i) * Bd;
        AB.block<N_S, N_S>(i*N_S, 0) = AiB;
    }


    for(int i=0; i<mpcWindow; i++)
    {
        for(int j=0; j<=i; j++)
        {
            G.block<1,N_S>(i*N_S, 4*j-3) = -K*(Eigen::Matrix<double, N_S, N_S>::Identity() - AB.block<N_S, N_S>(N_S*(i-j), 0)); 
        }
    }

    if(verbose)
    {
        std::cout << "[MPC API]\t G created."   << std::endl;
        std::cout << "G:" << std::endl << G << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI::computeGbar()
{
    Eigen::Matrix<double, 2*mpcWindow, N_C*mpcWindow> Gbar_temp;
    Gbar_temp << G, -G;
    Gbar.resize(2*mpcWindow, N_C*mpcWindow);

    for(int i=0; i<2*mpcWindow; i++)
    {
        for(int j=0; j<N_C*mpcWindow; j++)
        {
            Gbar.insert(i,j) = Gbar_temp(i,j);
        }
    }

    Gbar.makeCompressed();

    if(verbose)
    {
        std::cout << "[MPC API]\t Gbar created."   << std::endl;
        std::cout << "Gbar:" << std::endl << Gbar_temp << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI::setF()
{
    Eigen::Map<Eigen::Matrix<double, N_C*mpcWindow, 1>> ref_vector(ref.data(), ref.size());
    f = Fx*X+Fr*ref_vector;

    if(verbose)
    {
        std::cout << "[MPC API]\t f created."   << std::endl;
        std::cout << "f:" << std::endl << f << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI::updateRef()
{
    t = linspace(t0, t0 + dt/1000.0*(mpcWindow-1), mpcWindow);
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
