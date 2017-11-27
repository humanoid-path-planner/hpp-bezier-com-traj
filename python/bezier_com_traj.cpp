#include "bezier-com-traj/solve.hh"
#include <eigenpy/memory.hpp>
#include <eigenpy/eigenpy.hpp>

#include <boost/python.hpp>
#include <boost/python/exception_translator.hpp>

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::ContactData)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::ProblemData)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::ResultData)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::ResultDataCOMTraj)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::bezier_t)

namespace bezier_com_traj
{
using namespace boost::python;

ResultDataCOMTraj* zeroStepCapturability(centroidal_dynamics::Equilibrium* eq, const Vector3& com ,const Vector3& dCom ,const Vector3& l0, const bool useAngMomentum
                              , const double timeDuration, const double timeStep)
{
    bezier_com_traj::ContactData data;
    data.contactPhase_ = eq;
    bezier_com_traj::ProblemData pData;
    pData.c0_  = com;
    pData.dc0_ = dCom;
    pData.l0_  = l0;
    pData.contacts_.push_back(data);
    pData.useAngularMomentum_ = useAngMomentum;
    std::vector<double> Ts;
    Ts.push_back(timeDuration);
    ResultDataCOMTraj  res = solve0step(pData, Ts, timeStep);
    return new ResultDataCOMTraj(res);
}




struct res_data_exception : std::exception
{
  char const* what() const throw() { return "attributes not accessible for false resData"; }
};

void translate(const res_data_exception & e)
{
    // Use the Python 'C' API to set up an exception object
    PyErr_SetString(PyExc_RuntimeError, e.what());
}


VectorX get_xD(const ResultData& res)
{
    if (res.x.size() > 0)
        return res.x;
    std::cout << "x is not defined" << std::endl;
    throw res_data_exception();
}

double get_costD(const ResultData& res)
{
    return res.cost_;
}

bool get_succD(const ResultData& res)
{
    return res.success_;
}

bezier_t* getC_of_t(const ResultDataCOMTraj& res)
{
    const bezier_t * curve = res.constC_of_t();
    if(curve)
        return new bezier_t(*curve);
    std::cout << "x is not defined" << std::endl;
    throw res_data_exception();
}

bezier_t* getDL_of_t(const ResultDataCOMTraj& res)
{
    const bezier_t * curve = res.constDL_of_t();
    if(curve)
        return new bezier_t(*curve);
    std::cout << "x is not defined" << std::endl;
    throw res_data_exception();
}

VectorX get_x(const ResultDataCOMTraj& res)
{
    if (res.x.size() > 0)
        return res.x;
    std::cout << "x is not defined" << std::endl;
    throw res_data_exception();
}

double get_cost(const ResultDataCOMTraj& res)
{
    return res.cost_;
}

bool get_succ(const ResultDataCOMTraj& res)
{
    return res.success_;
}

BOOST_PYTHON_MODULE(bezier_com_traj)
{
    using namespace boost::python;
    register_exception_translator<res_data_exception>(&translate);
    /** BEGIN eigenpy init**/
    eigenpy::enableEigenPy();

    eigenpy::enableEigenPySpecific<Vector3,Vector3>();
    eigenpy::enableEigenPySpecific<VectorX,VectorX>();
    /*eigenpy::exposeAngleAxis();
    eigenpy::exposeQuaternion();*/

    class_<ResultDataCOMTraj>("ResultDataCOMTraj", init<> ())
                .add_property("c_of_t",  make_function(&getC_of_t,
                                                       return_value_policy<manage_new_object>()))
                .add_property("dL_of_t",make_function(&getDL_of_t,
                                                       return_value_policy<manage_new_object>()))
                .add_property("success", &get_succ)
                .add_property("cost",    &get_cost)
                .add_property("x",       &get_x)
            ;


    class_<ResultData>("ResultData", init<> ())
                .add_property("success", &get_succD)
                .add_property("cost",    &get_costD)
                .add_property("x",       &get_xD)
            ;

    def("zeroStepCapturability", &zeroStepCapturability, return_value_policy<manage_new_object>());

    /** END eigenpy init**/

}

} // namespace bezier_com_traj
