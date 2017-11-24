
#include <bezier-com-traj/common_solve_methods.hh>


using namespace centroidal_dynamics;


Matrix3 skew(point_t_tC x)
{
    Matrix3 res = Matrix3::Zero();
    res(0,1) = - x(2); res(0,2) =   x(1);
    res(1,0) =   x(2); res(1,2) = - x(0);
    res(2,0) = - x(1); res(2,1) =   x(0);
    return res;
}

waypoint_t centroidal_dynamics::initwp()
{
    waypoint_t w;
    w.first  = matrix6_t::Zero();
    w.second = point6_t::Zero();
    return w;
}


std::vector<waypoint_t> ComputeDiscretizedWaypoints(const std::vector<waypoint_t>& wps, const std::vector<spline::Bern<double> >& berns,  int numSteps)
{
    double dt = 1./numSteps;
    std::vector<waypoint_t> res;
    for (int i =0; i < numSteps + 1; ++i)
    {
        waypoint_t w = initwp();
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

MatrixXX initMatrixA(const int dimH, const std::vector<waypoint_t>& wps, const int extraConstraintSize, const bool useAngMomentum)
{
    int dimPb = useAngMomentum ? 6 : 3;
    return MatrixXX::Zero(dimH * wps.size() + extraConstraintSize, dimPb);
}

void addKinematic(Ref_matrixXX& A, Ref_vectorX& b, Cref_matrixX3 Kin, Cref_vectorX kin, const int otherConstraintIndex)
{
    int dimKin = (int)(kin.rows());
    if (dimKin == 0) return;
    A.block(A.rows() - dimKin - otherConstraintIndex ,0, dimKin, 3) = Kin;
    b.tail(dimKin) = kin;
}

void addAngularMomentum(Ref_matrixXX& A, Ref_vectorX& b, Cref_matrixX3 Ang, Cref_vectorX ang)
{
    int dimAng = (int)(ang.rows());
    A.block(A.rows() - dimAng  , 3, dimAng, 3) = Ang;
    b.tail(dimAng) = ang;
}
/*
 *
 *
def __compute_uixs(self, l0, T, num_step = -1):
        alpha = 1. / (T)
        wps = [ui(l0, alpha)  for ui in uis]
        if num_step > 0:
            dt = (1./float(num_step))
            wps_bern = [ [ (b(i*dt)*wps[idx][0], b(i*dt)*wps[idx][1]) for idx,b in enumerate(b4)] for i in range(num_step + 1) ]
            wps = [reduce(lambda a, b : (a[0] + b[0], a[1] + b[1]), wps_bern_i) for wps_bern_i in wps_bern]
        return wps
 *
 * ups = self.__compute_uixs(l0, T ,num_steps)
        AL, bL = self._init_matrices_AL_bL(ups, A, b)
        dimH  = self._H.shape[0]
        #final matrix has num rows equal to initial matrix rows + angular momentum constraints
        # the angular momentum constraints are added AFTER the eventual kinematic ones
        for i, (uxi, usi) in enumerate(ups):
            AL[i*dimH : (i+1)*dimH, 3:]  = self._H.dot(uxi) #constant part of A, Ac = Ac * wxi
            bL[i*dimH : (i+1)*dimH    ] += self._H.dot(-usi)

        if self._angular_momentum_constraints != None:
            dimL = self._angular_momentum_constraints[0].shape[0]
            AL[-dimL:,3:] = self._angular_momentum_constraints[0][:]
            bL[-dimL:   ] = self._angular_momentum_constraints[1][:]

        AL, bL = normalize(AL,bL)
 */


size_t removeZeroRows(Ref_matrixXX& A, Ref_vectorX& b)
{
    Eigen::Matrix<bool, Eigen::Dynamic, 1> empty = (A.array() == 0).rowwise().all();
    size_t last = A.rows() - 1;
    for (size_t i = 0; i < last + 1;)
    {
        if (empty(i))
        {
            assert(b(i) * b(i) < 10e-6);
            A.row(i).swap(A.row(last));
            b.row(i).swap(b.row(last));
            empty.segment<1>(i).swap(empty.segment<1>(last));
            --last;
        }
        else
            ++i;
    }
    return last+1;
}

size_t Normalize(Ref_matrixXX& A, Ref_vectorX& b)
{
    size_t zeroindex = removeZeroRows(A,b);
    Eigen::VectorXd norm = A.block(0,0,zeroindex,A.cols()).rowwise().norm();
    A.block(0,0,zeroindex,A.cols()).rowwise().normalize();
    b.head(zeroindex) = b.head(zeroindex).cwiseQuotient(norm);
    return zeroindex;
}

std::pair<MatrixXX, VectorX> compute6dControlPointInequalities(const ContactData& cData, const std::vector<waypoint_t>& wps, const std::vector<waypoint_t>& wpL, const bool useAngMomentum, double T, double timeStep)
{
    std::pair<MatrixXX, VectorX> res;
    Ref_matrixXX A = res.first;
    Ref_vectorX  b = res.second;
    // gravity vector
    const point_t& g = cData.contactPhase_->m_gravity;
    // compute GIWC
    MatrixXX H; VectorX h;
    cData.contactPhase_->getPolytopeInequalities(H,h);
    H = -H;
    int dimH = (int)(H.rows());
    MatrixXX mH = cData.contactPhase_->m_mass * H;
    // init and fill Ab matrix
    int dimKin =  cData.kin_ == point_t::Zero() ? 0 : (int)(cData.kin_.rows());
    int dimAng = (cData.ang_ != point_t::Zero() && useAngMomentum)  ? (int)(cData.ang_.rows()) : 0;
    A = initMatrixA(dimH, wps, dimKin + dimAng, useAngMomentum);
    b = VectorX::Zero(A.rows());
    point6_t bc = point6_t::Zero(); bc.head(3) = g; // constant part of Aub, Aubi = mH * (bc - wsi)
    int i = 0;
    std::vector<waypoint_t>::const_iterator wpLcit = wpL.begin();
    for (std::vector<waypoint_t>::const_iterator wpcit = wps.begin(); wpcit != wps.end(); ++wpcit)
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
    // normalization removes 0 value rows, but resizing
    // must actually be done with matrices and not the references
    size_t zeroindex = Normalize(A,b);
    res.first.conservativeResize(zeroindex, A.cols());
    res.second.conservativeResize(zeroindex, 1);
    return res;
}
