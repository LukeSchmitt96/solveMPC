#include "./ModelPredictiveControlAPI.h"

ModelPredictiveControlAPI::ModelPredictiveControlAPI(bool verbose_)
{

    std::cout << "[MPC API]\tMPC API object created." << std::endl;

    verbose = verbose_;
    solverFlag = true;

    // get stream from config file and parse
    std::ifstream file("./config/MPC_API.json");
    cfg = json::parse(file);

    // set K from cfg
    K = from_json(cfg["K"], K.rows(), K.cols());

    // set xref from cfg
    xref = cfg["xref"].get<double>();

    // set arbitraty/to be filled variables
    X << 0.0, 0.0, 0.0, 0.0;
    U << 0.0;
    t0 = 0.0;
    dt = 0.0;

    // set the matrices neccessary to solve the qp problem
    setSystemVars();
    setCosts();
    setLiftedCosts();
    setTransformations();
    setLL();
    setH();
    setLu();
    setFVars();
    setLinearConstraints();
    setUpperBound();
    updateRef(xref);        // assumes constant reference
    setF();

    // lb will be constant as our problem constraints are formulated as Ax <= b
    lb = Eigen::Matrix<double, 2*mpcWindow, 1>::Ones() * -std::numeric_limits<double>::max();
    ub = W0 + Sbar*X + Ku*U;

    std::cout << "[MPC API]\tAll QP matrices built successfully." << std::endl;

    n_variables = N_O*mpcWindow;
    n_constraints = 2*mpcWindow;

    // solver settings
    solver.settings()->setVerbosity(verbose);
    solver.settings()->setWarmStart(true);

    solver.data()->setNumberOfVariables(n_variables);
    solver.data()->setNumberOfConstraints(n_constraints);

    if(!solver.data()->setHessianMatrix(H))                 {solverFlag = false; return;};
    if(!solver.data()->setGradient(f))                      {solverFlag = false; return;};
    if(!solver.data()->setLinearConstraintsMatrix(Gbar))    {solverFlag = false; return;};
    if(!solver.data()->setLowerBound(lb))                   {solverFlag = false; return;};
    if(!solver.data()->setUpperBound(ub))                   {solverFlag = false; return;};

    // instantiate the solver
    if(!solver.initSolver()){solverFlag = false; return;};
}


ModelPredictiveControlAPI::~ModelPredictiveControlAPI()
{
    printf("[MPC API]\tDestructing MPC API object...\n");
}


void ModelPredictiveControlAPI::setVerbosity(bool verbose_)
{
    verbose = verbose_;
    std::cout << "[MPC API]\tVerbosity set to " << verbose << std::endl;
}


bool ModelPredictiveControlAPI::controllerStep()
{
    // update time vector
    t0 += dt;

    // recalculate reference (for now just hold at xref)
    updateRef(xref);

    // calculate new F based on ref, U and X
    setF();

    // recalculate upper bound constraints
    setUpperBound();

    // update F
    if(!solver.updateGradient(f)) return false;

    // update upper bound based on X            W0+S*x+Ku*U,
    if(!solver.updateUpperBound(W0 + Sbar*X + Ku*U)) return false;

    // solve the QP problem
    if(!solver.solve()) return false;

    // get the next controller input
    U += solver.getSolution().block<N_C, 1>(0, 0);

    return true;
}


void ModelPredictiveControlAPI::setSystemVars()
{
    Ad = from_json(cfg["Ad"], Ad.rows(), Ad.cols());
    Bd = from_json(cfg["Bd"], Bd.rows(), Bd.cols());
    Cd = from_json(cfg["Cd"], Cd.rows(), Cd.cols());
    Dd = from_json(cfg["Dd"], Dd.rows(), Dd.cols());

    if(verbose)
    {
        std::cout << "[MPC API]\tSystem variables created."   << std::endl;

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

void ModelPredictiveControlAPI::setCosts()
{
    Q  = from_json(cfg["Q"],  Q.rows(),  Q.cols());
    R  = from_json(cfg["R"],  R.rows(),  R.cols());
    RD = from_json(cfg["RD"], RD.rows(), RD.cols());

    if(verbose)
    {
        std::cout << "[MPC API]\tSet Q, R, and RD matrices created."  << std::endl;

        std::cout << "Q rows: " << Q.rows() << "\tQ cols: " << Q.cols()  << std::endl;
        std::cout << "Q:"  << std::endl << Q  << std::endl << std::endl;

        std::cout << "R rows: " << R.rows() << "\tR cols: " << R.cols()  << std::endl;
        std::cout << "R:"  << std::endl << R  << std::endl << std::endl;

        std::cout << "RD rows: " << RD.rows() << "\tRD cols: " << RD.cols()  << std::endl;
        std::cout << "RD:" << std::endl << RD << std::endl << std::endl;
    }
}


void ModelPredictiveControlAPI::setLiftedCosts()
{
    Qbar = blkdiag(Q, mpcWindow);
    Rbar = blkdiag(R, mpcWindow);
    RbarD = blkdiag(RD, mpcWindow);

    if(verbose)
    {
        std::cout << "[MPC API]\tLifted weight matrices created."   << std::endl;

        std::cout << "Qbar rows: " << Qbar.rows() << "\tQbar cols: " << Qbar.cols()  << std::endl;
        std::cout << "Qbar:"  << std::endl << Qbar  << std::endl << std::endl;

        std::cout << "Rbar rows: " << Rbar.rows() << "\tRbar cols: " << Rbar.cols()  << std::endl;
        std::cout << "Rbar:"  << std::endl << Rbar  << std::endl << std::endl;

        std::cout << "RbarD rows: " << RbarD.rows() << "\tRbarD cols: " << RbarD.cols()  << std::endl;
        std::cout << "RbarD:" << std::endl << RbarD << std::endl << std::endl;
    }
}


void ModelPredictiveControlAPI::setTransformations()
{
    CAiB = Eigen::MatrixXd::Zero(N_S, N_O);
    Su_full = Eigen::MatrixXd::Zero(Su_full.rows(), Su_full.cols());

    for(int i=0; i<10; i++) {S.block<1,N_S>(i,0) = K;}

    for(int i=0; i<mpcWindow; i++)
    {
        Sx.block<N_O, N_S>(i*N_O,0) = Cd*Ad.pow(i+1);
        CAB.block<N_C,N_O>(i*N_O,0) = Cd*Ad.pow(i)*Bd;

        CAiB += Eigen::Matrix<double, N_S, N_S>::Identity() * Ad.pow(i) * Bd;
        CAB_full.block<N_S, N_O>(i*N_S,0) = CAiB;
    }
    Eigen::VectorXd CAB_Vector(Eigen::Map<Eigen::VectorXd>(CAB.data(), CAB.cols() * CAB.rows()));

    for(int i=0; i<mpcWindow; i++)
    {
        for(int j=0; j<=i; j++)
        {
            Su(i,j) = CAB_Vector.segment(0,i-j+1).sum();
            Su_full.block<N_S, 1>(N_S*i, j) = CAB_full.block(4*(i-j), 0, N_S, N_O);
        }
    }

    Su1 = Su.leftCols<N_O>();
    Su_full1 = Su_full.leftCols<N_O>();
    Sbar << S, -S;

    if(verbose)
    {
        std::cout << "[MPC API]\tTransformation matrices created" << std::endl;

        std::cout << "Sx rows: " << Sx.rows() << "\tSx cols: " << Sx.cols()  << std::endl;
        std::cout << "Sx:"    << std::endl << Sx    << std::endl << std::endl;

        std::cout << "CAB rows: " << CAB.rows() << "\tCAB cols: " << CAB.cols()  << std::endl;
        std::cout << "CAB:"   << std::endl << CAB   << std::endl << std::endl;

        std::cout << "CAiB rows: " << CAiB.rows() << "\tCAiB cols: " << CAiB.cols()  << std::endl;
        std::cout << "CAiB:"   << std::endl << CAiB   << std::endl << std::endl;

        std::cout << "CAB_full rows: " << CAB_full.rows() << "\tCAB_full cols: " << CAB_full.cols()  << std::endl;
        std::cout << "CAB_full:"   << std::endl << CAB_full   << std::endl << std::endl;

        std::cout << "S rows: " << S.rows() << "\tS cols: " << S.cols()  << std::endl;
        std::cout << "S:" << std::endl << S << std::endl  << std::endl;

        std::cout << "Su rows: " << Su.rows() << "\tSu cols: " << Su.cols()  << std::endl;
        std::cout << "Su:"    << std::endl << Su    << std::endl << std::endl;

        std::cout << "Su1 rows: " << Su1.rows() << "\tSu1 cols: " << Su1.cols()  << std::endl;
        std::cout << "Su1:"   << std::endl << Su1   << std::endl << std::endl;

        std::cout << "Su_full rows: " << Su_full.rows() << "\tSu1 cols: " << Su_full.cols()  << std::endl;
        std::cout << "Su_full:"   << std::endl << Su_full   << std::endl << std::endl;

        std::cout << "Su_full1 rows: " << Su_full1.rows() << "\tSu1 cols: " << Su_full1.cols()  << std::endl;
        std::cout << "Su_full1:"   << std::endl << Su_full1   << std::endl << std::endl;

        std::cout << "Sbar rows: " << Sbar.rows() << "\tSbar cols: " << Sbar.cols()  << std::endl;
        std::cout << "Sbar:" << std::endl << Sbar << std::endl  << std::endl;        
    }
}


void ModelPredictiveControlAPI::setH()
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
        std::cout << "[MPC API]\tHessian H created." << std::endl;
        std::cout << "H rows: " << H_temp_1.rows() << "\tH cols: " << H_temp_1.cols()  << std::endl;
        std::cout << "H:" << std::endl << H_temp_1 << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI::setLu()
{
    for(int i=0; i<mpcWindow; i++)
    {
        Lu.block<N_C,N_C>(i,0) = (mpcWindow-i+2)*Eigen::Matrix<double,N_C,N_C>::Identity();
    }

    if(verbose)
    {
        std::cout << "[MPC API]\tLu created." << std::endl;
        std::cout << "Lu rows: " << Lu.rows() << "\tLu cols: " << Lu.cols()  << std::endl;
        std::cout << "Lu:" << std::endl << Lu << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI:: setLL()
{    
    LL = Eigen::MatrixXd(Eigen::Matrix<double, N_O*mpcWindow, N_O*mpcWindow>::Ones().triangularView<Eigen::Lower>());

    if(verbose)
    {
        std::cout << "[MPC API]\tLL created."   << std::endl;
        std::cout << "LL rows: " << LL.rows() << "\tLL cols: " << LL.cols()  << std::endl;
        std::cout << "LL:" << std::endl << LL << std::endl  << std::endl;
    }
} 


void ModelPredictiveControlAPI::setFVars()
{
    Fu = 2*((LL.transpose() * Rbar.transpose()).diagonal().transpose() + Su1.transpose() * Qbar * Su).transpose();
    Fr = -2*(Qbar*Su).transpose();
    Fx = 2*(Sx.transpose() * Qbar * Su).transpose();

    if(verbose)
    {
        std::cout << "[MPC API]\tComponents of F created."   << std::endl;

        std::cout << "Fu rows: " << Fu.rows() << "\tFu cols: " << Fu.cols()  << std::endl;
        std::cout << "Fu:" << std::endl << Fu << std::endl  << std::endl;

        std::cout << "Fr rows: " << Fr.rows() << "\tFr cols: " << Fr.cols()  << std::endl;
        std::cout << "Fr:" << std::endl << Fr << std::endl  << std::endl;

        std::cout << "Fx rows: " << Fx.rows() << "\tFx cols: " << Fx.cols()  << std::endl;
        std::cout << "Fx:" << std::endl << Fx << std::endl  << std::endl;
    }
    
}


void ModelPredictiveControlAPI::setLinearConstraints()
{
    Eigen::Matrix<double, 2 * mpcWindow, N_C * mpcWindow> Gbar_temp;
    Eigen::Matrix<double, mpcWindow, N_C * mpcWindow> Gbar_temp_upper;
    Eigen::Matrix<double, mpcWindow, N_C * mpcWindow> Gbar_temp_lower;

    Gbar_temp_upper = Eigen::Matrix<double, mpcWindow, N_C * mpcWindow>::Ones().triangularView<Eigen::Lower>();
    Gbar_temp_lower = Eigen::Matrix<double, mpcWindow, N_C * mpcWindow>::Ones().triangularView<Eigen::Lower>();

    Gbar_temp << Gbar_temp_upper*K(0), Gbar_temp_lower*-K(0);

    Gbar.resize(Gbar_temp.rows(), Gbar_temp.cols());

    for(int i=0; i<Gbar_temp.rows(); i++)
    {
        for(int j=0; j<Gbar_temp.cols(); j++)
        {
            Gbar.insert(i,j) = Gbar_temp(i,j);
        }
    }

    Gbar.makeCompressed();


    if(verbose)
    {
        std::cout << "[MPC API]\tLinear constraints matrix created."   << std::endl;

        std::cout << "Gbar rows: " << Gbar_temp.rows() << "\tGbar cols: " << Gbar_temp.cols()  << std::endl;
        std::cout << "Gbar:" << std::endl << Gbar_temp << std::endl  << std::endl;
    }
}


void ModelPredictiveControlAPI::setUpperBound()
{
    Eigen::Matrix<double, mpcWindow, N_O> Ku_temp;

    Ku_temp = Eigen::Matrix<double, mpcWindow, N_O>::Ones()*K(0);

    Ku << -Ku_temp, Ku_temp;

    W0 = 255.0 * Eigen::Matrix<double, 2*mpcWindow, 1>::Ones();
}


void ModelPredictiveControlAPI::setF()
{
    f = Fx*X + Fu*U + Fr*ref.transpose();
}


void ModelPredictiveControlAPI::updateRef(double pos_ref)
{
    ref.block<1,mpcWindow>(0,0) = pos_ref * Eigen::Matrix<double, N_C, mpcWindow>::Ones();

    if(verbose)
    {
        std::cout << "[MPC API]\tref: " << ref << std::endl;
    }
}


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
    return linspaced;
}

Eigen::MatrixXd ModelPredictiveControlAPI::from_json(const nlohmann::json& jsonObject, int rows, int cols)
{
    Eigen::MatrixXd matrix(rows, cols);

    nlohmann::json jsonArray;
    if ( jsonObject.is_array( ) )
    {
        if ( jsonObject.empty( ) )
        {
            return matrix;
        }
        jsonArray = jsonObject;
    }
    else if ( jsonObject.is_number( ) )
    {
        jsonArray.push_back( jsonObject );
    }
    else
    {
        throw nlohmann::detail::type_error::create( 0, "" );
    }

    nlohmann::json jsonArrayOfArrays;
    if ( jsonArray.front( ).is_array( ) )  // provided matrix
    {
        jsonArrayOfArrays = jsonArray;
    }
    else  // provided vector
    {
        if ( rows == 1 )  // expected row vector
        {
            jsonArrayOfArrays.push_back( jsonArray );
        }
        else if ( cols == 1 )  // expected column vector
        {
            for ( unsigned int i = 0; i < jsonArray.size( ); ++i )
            {
                jsonArrayOfArrays.push_back( { jsonArray.at( i ) } );
            }
        }
        else  // expected matrix
        {
            std::cerr << "Expected a matrix, received a vector." << std::endl;
            throw nlohmann::detail::type_error::create( 0, "" );
        }
    }

    const unsigned int providedRows = jsonArrayOfArrays.size( );
    const unsigned int providedCols = jsonArrayOfArrays.front( ).size( );
    if ( ( rows >= 0 && int( providedRows ) != rows ) || ( cols >= 0 && int( providedCols ) != cols ) )
    {
        std::cerr << "Expected matrix of size " << rows << "x" << cols
                << ", received matrix of size " << providedRows << "x" << providedCols << "." << std::endl;
        throw nlohmann::detail::type_error::create( 0, "" );
    }

    matrix.resize( providedRows, providedCols );
    for ( unsigned int r = 0; r < providedRows; ++r )
    {
        if ( jsonArrayOfArrays.at( r ).size( ) != providedCols )
        {
            std::cerr << "Unconsistent matrix size: some rows have different number of columns." << std::endl;
            throw nlohmann::detail::type_error::create( 0, "" );
        }
        for ( unsigned int c = 0; c < providedCols; ++c )
        {
            matrix( r, c ) = jsonArrayOfArrays.at( r ).at( c );
        }
    }

    return matrix;
}
