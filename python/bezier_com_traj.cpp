#include "bezier-com-traj/solve.hh"
#include <eigenpy/memory.hpp>
#include <eigenpy/eigenpy.hpp>

#include <boost/python.hpp>

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::ContactData)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::ProblemData)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::ResultDataCOMTraj)

namespace bezier_com_traj
{
using namespace boost::python;

ResultDataCOMTraj* zeroStepCapturability(Equilibrium& eq, const Vector3& com ,const Vector3& dCom ,const Vector3& l0, const bool useAngMomentum
                              , const double timeDuration, const double timestep)
{
    bezier_com_traj::ContactData data;
    data.contactPhase_ = eq;
    bezier_com_traj::ProblemData pData;
    pData.c0_ = com;
    pData.dc0_ = dCom;
    pData.contacts_.push_back(data);
    std::vector<double> Ts;
    Ts.pushback(timeDuration);
    ResultDataCOMTraj res = solve0step( pData, Ts, timeStep);
    ResultDataCOMTraj* ret = new ResultDataCOMTraj();
}

boost::python::tuple wrapComputeEquilibriumRobustness(Equilibrium& self, const Vector3& com, const Vector3& acc)
{
    double robustness;
    LP_status status = self.computeEquilibriumRobustness(com, acc, robustness);
    return boost::python::make_tuple(status, robustness);
}

boost::python::tuple wrapGetPolytopeInequalities(Equilibrium& self)
{
    MatrixXX H;
    VectorX h;
    self.getPolytopeInequalities(H,h);
    MatrixXXColMajor _H = H;
    return boost::python::make_tuple(_H, h);
}



BOOST_PYTHON_MODULE(bezier_com_traj)
{
    /** BEGIN eigenpy init**/
    eigenpy::enableEigenPy();

    eigenpy::enableEigenPySpecific<Vector3,Vector3>();
    eigenpy::enableEigenPySpecific<VectorX,VectorX>();
    /*eigenpy::exposeAngleAxis();
    eigenpy::exposeQuaternion();*/

    /** END eigenpy init**/

    /** BEGIN enum types **/
  #ifdef CLP_FOUND
    enum_<SolverLP>("SolverLP")
            .value("SOLVER_LP_QPOASES", SOLVER_LP_QPOASES)
            .value("SOLVER_LP_CLP", SOLVER_LP_CLP)
            .export_values();
  #else
    enum_<SolverLP>("SolverLP")
            .value("SOLVER_LP_QPOASES", SOLVER_LP_QPOASES)
            .export_values();
  #endif


    enum_<EquilibriumAlgorithm>("EquilibriumAlgorithm")
            .value("EQUILIBRIUM_ALGORITHM_LP", EQUILIBRIUM_ALGORITHM_LP)
            .value("EQUILIBRIUM_ALGORITHM_LP2", EQUILIBRIUM_ALGORITHM_LP2)
            .value("EQUILIBRIUM_ALGORITHM_DLP", EQUILIBRIUM_ALGORITHM_DLP)
            .value("EQUILIBRIUM_ALGORITHM_PP", EQUILIBRIUM_ALGORITHM_PP)
            .value("EQUILIBRIUM_ALGORITHM_IP", EQUILIBRIUM_ALGORITHM_IP)
            .value("EQUILIBRIUM_ALGORITHM_DIP", EQUILIBRIUM_ALGORITHM_DIP)
            .export_values();

    enum_<LP_status>("LP_status")
            .value("LP_STATUS_UNKNOWN", LP_STATUS_UNKNOWN)
            .value("LP_STATUS_OPTIMAL", LP_STATUS_OPTIMAL)
            .value("LP_STATUS_INFEASIBLE", LP_STATUS_INFEASIBLE)
            .value("LP_STATUS_UNBOUNDED", LP_STATUS_UNBOUNDED)
            .value("LP_STATUS_MAX_ITER_REACHED", LP_STATUS_MAX_ITER_REACHED)
            .value("LP_STATUS_ERROR", LP_STATUS_ERROR)
            .export_values();

    /** END enum types **/

    bool (Equilibrium::*setNewContacts)
            (const MatrixX3ColMajor&, const MatrixX3ColMajor&, const double, const EquilibriumAlgorithm) = &Equilibrium::setNewContacts;

    class_<Equilibrium>("Equilibrium", init<std::string, double, unsigned int, optional <SolverLP, bool, const unsigned int, const bool> >())
            .def("useWarmStart", &Equilibrium::useWarmStart)
            .def("setUseWarmStart", &Equilibrium::setUseWarmStart)
            .def("getName", &Equilibrium::getName)
            .def("getAlgorithm", &Equilibrium::getAlgorithm)
            .def("setNewContacts", setNewContacts)
            .def("computeEquilibriumRobustness", wrapComputeQuasiEquilibriumRobustness)
            .def("computeEquilibriumRobustness", wrapComputeEquilibriumRobustness)
            .def("getPolytopeInequalities", wrapGetPolytopeInequalities)
    ;
}

} // namespace centroidal_dynamics
