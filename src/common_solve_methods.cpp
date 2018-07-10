
#include <bezier-com-traj/common_solve_methods.hh>
#include <solver/eiquadprog-fast.hpp>
#include <solver/glpk-wrapper.hpp>
#include <bezier-com-traj/waypoints/waypoints_definition.hh>

namespace bezier_com_traj
{
std::vector<waypoint6_t> ComputeDiscretizedWaypoints(const std::vector<waypoint6_t>& wps, const std::vector<spline::Bern<double> >& berns,  int numSteps)
{
    double dt = 1./double(numSteps);
    std::vector<waypoint6_t> res;
    for (int i =0; i < numSteps + 1; ++i)
    {
        waypoint6_t w = initwp<waypoint6_t>();
        for (int j = 0; j <5; ++j)
        {
            double b = berns[j](i*dt);
            w.first +=b*(wps[j].first );
            w.second+=b*(wps[j].second);
        }
        res.push_back(w);
    }
    return res;
}

MatrixXX initMatrixA(const int dimH, const std::vector<waypoint6_t>& wps, const long int extraConstraintSize, const bool useAngMomentum)
{
    int dimPb = useAngMomentum ? 6 : 3;
    return MatrixXX::Zero(dimH * wps.size() + extraConstraintSize, dimPb);
}

MatrixXX initMatrixD(const int colG, const std::vector<waypoint6_t>& wps, const long int extraConstraintSize, const bool useAngMomentum)
{
    int dimPb = useAngMomentum ? 6 : 3;
    return MatrixXX::Zero(6 + extraConstraintSize, wps.size() * colG +  dimPb);
}

void addKinematic(Ref_matrixXX A, Ref_vectorX b, Cref_matrixX3 Kin, Cref_vectorX kin, const long int otherConstraintIndex)
{
    int dimKin = (int)(kin.rows());
    if (dimKin == 0) return;
    A.block(A.rows() - dimKin - otherConstraintIndex ,0, dimKin, 3) = Kin;
    b.segment(A.rows() - dimKin - otherConstraintIndex, dimKin) = kin;
}

void addAngularMomentum(Ref_matrixXX A, Ref_vectorX b, Cref_matrixX3 Ang, Cref_vectorX ang)
{
    int dimAng = (int)(ang.rows());
    if(dimAng > 0)
    {
        A.block(A.rows() - dimAng  , 3, dimAng, 3) = Ang;
        b.tail(dimAng) = ang;
    }
}

int removeZeroRows(Ref_matrixXX& A, Ref_vectorX& b)
{
    Eigen::Matrix<bool, Eigen::Dynamic, 1> empty = (A.array() == 0).rowwise().all();
    size_t last = A.rows() - 1;
    for (size_t i = 0; i < last + 1;)
    {
        if (empty(i))
        {
            assert(A.row(i).norm() == 0);
            if(b(i) <0)
            {
                //std::cout << "empty row for A while correponding b is negative. Problem is infeasible" << b(i) << std::endl;
                return -1;
            }
            A.row(i).swap(A.row(last));
            b.row(i).swap(b.row(last));
            empty.segment<1>(i).swap(empty.segment<1>(last));
            --last;
        }
        else
            ++i;
    }
    return (int)last+1;
}

int Normalize(Ref_matrixXX A, Ref_vectorX b)
{
    int zeroindex = removeZeroRows(A,b);
    if(zeroindex > 0)
    {
        Eigen::VectorXd norm = A.block(0,0,zeroindex,A.cols()).rowwise().norm();
        A.block(0,0,zeroindex,A.cols()).rowwise().normalize();
        b.head(zeroindex) = b.head(zeroindex).cwiseQuotient(norm);
    }
    return zeroindex;
}



std::pair<MatrixXX, VectorX> compute6dControlPointInequalities(const ContactData& cData, const std::vector<waypoint6_t>& wps, const std::vector<waypoint6_t>& wpL, const bool useAngMomentum, bool& fail)
{
    MatrixXX A;
    VectorX  b;
    // gravity vector
    const point_t& g = cData.contactPhase_->m_gravity;
    // compute GIWC
    assert (cData.contactPhase_->getAlgorithm() ==  centroidal_dynamics::EQUILIBRIUM_ALGORITHM_PP);
    centroidal_dynamics::MatrixXX Hrow; VectorX h;
    cData.contactPhase_->getPolytopeInequalities(Hrow,h);
    MatrixXX H = -Hrow;
    H.rowwise().normalize();
    int dimH = (int)(H.rows());
    MatrixXX mH = cData.contactPhase_->m_mass * H;
    // init and fill Ab matrix
    long int dimKin =  cData.kin_.size();
    long int dimAng = useAngMomentum  ? (long int)(cData.ang_.size()) : 0;
    A = initMatrixA(dimH, wps, dimKin + dimAng, useAngMomentum);
    b = VectorX::Zero(A.rows());
    point6_t bc = point6_t::Zero(); bc.head(3) = g; // constant part of Aub, Aubi = mH * (bc - wsi)
    int i = 0;
    std::vector<waypoint6_t>::const_iterator wpLcit = wpL.begin();
    for (std::vector<waypoint6_t>::const_iterator wpcit = wps.begin(); wpcit != wps.end(); ++wpcit)
    {
        A.block(i*dimH,0, dimH, 3) = mH * wpcit->first;
        b.segment(i*dimH, dimH)    = mH * (bc - wpcit->second);
        if(useAngMomentum)
        {
            A.block  (i*dimH,3, dimH, 3) = H * wpLcit->first;
            b.segment(i*dimH, dimH)     += H * (-wpLcit->second);
            ++wpLcit;
        }
        ++i;
    }
    addKinematic(A,b, cData.Kin_,cData.kin_, dimAng);
    if (useAngMomentum)
        addAngularMomentum(A,b,cData.Kin_, cData.kin_);
    fail = false;
    // normalization removes 0 value rows, but resizing
    // must actually be done with matrices and not the references
    /*int zeroindex = Normalize(A,b);
    if(zeroindex < 0)
        fail = true;
    else
    {
        A.conservativeResize(zeroindex, A.cols());
        b.conservativeResize(zeroindex, 1);
        fail = false;
    }*/
    return std::make_pair(A,b);
}

ResultData solve(Cref_matrixXX A, Cref_vectorX ci0, Cref_matrixXX D, Cref_vectorX d, Cref_matrixXX H,
                 Cref_vectorX g, Cref_vectorX initGuess, const solvers::SolverType solver)
{
    return solvers::solve(A,ci0,D,d,H,g,initGuess,solver);
}

ResultData solve(Cref_matrixXX A, Cref_vectorX b, Cref_matrixXX H, Cref_vectorX g, Cref_vectorX initGuess, const solvers::SolverType solver)
{
    MatrixXX D = MatrixXX::Zero(0,A.cols());
    VectorX d  = VectorX::Zero(0);
    return solvers::solve(A,b,D,d,H,g,initGuess,solver);
}


ResultData solve(const std::pair<MatrixXX, VectorX>& Ab,const std::pair<MatrixXX, VectorX>& Hg,  const VectorX& init, const solvers::SolverType solver)
{
    return solve(Ab.first,Ab.second,Hg.first,Hg.second, init,solver);
}


ResultData solve(const std::pair<MatrixXX, VectorX>& Ab,const std::pair<MatrixXX, VectorX>& Dd,const std::pair<MatrixXX, VectorX>& Hg,  const VectorX& init, const solvers::SolverType solver)
{
    return solve(Ab.first,Ab.second,Dd.first, Dd.second, Hg.first,Hg.second, init, solver);
}

} // namespace bezier_com_traj
