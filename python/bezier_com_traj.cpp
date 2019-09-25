#include <hpp/bezier-com-traj/solve.hh>
#include <hpp/bezier-com-traj/solver/solver-abstract.hpp>
#include <hpp/bezier-com-traj/solve_end_effector.hh>
#include <eigenpy/memory.hpp>
#include <eigenpy/eigenpy.hpp>

#include <boost/python.hpp>
#include <boost/python/exception_translator.hpp>

EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::ContactData)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::ProblemData)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::ResultData)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::ResultDataCOMTraj)
EIGENPY_DEFINE_STRUCT_ALLOCATOR_SPECIALIZATION(bezier_com_traj::bezier_t)

namespace bezier_com_traj {
using namespace boost::python;
typedef double real;

ResultDataCOMTraj* zeroStepCapturability(centroidal_dynamics::Equilibrium* eq, const Vector3& com, const Vector3& dCom,
                                         const Vector3& l0, const bool useAngMomentum, const double timeDuration,
                                         const double timeStep) {
  bezier_com_traj::ContactData data;
  data.contactPhase_ = eq;
  bezier_com_traj::ProblemData pData;
  pData.c0_ = com;
  pData.dc0_ = dCom;
  pData.l0_ = l0;
  pData.contacts_.push_back(data);
  pData.useAngularMomentum_ = useAngMomentum;
  std::vector<double> Ts;
  Ts.push_back(timeDuration);
  ResultDataCOMTraj res = solve0step(pData, Ts, timeStep);
  return new ResultDataCOMTraj(res);
}

ResultDataCOMTraj* zeroStepCapturabilityWithKinConstraints(centroidal_dynamics::Equilibrium* eq, const Vector3& com,
                                                           const Vector3& dCom, const Vector3& l0,
                                                           const bool useAngMomentum, const double timeDuration,
                                                           const double timeStep, const MatrixXX& Kin,
                                                           const MatrixXX& kin) {
  bezier_com_traj::ContactData data;
  data.Kin_ = Kin;
  data.kin_ = kin;
  data.contactPhase_ = eq;
  bezier_com_traj::ProblemData pData;
  pData.c0_ = com;
  pData.dc0_ = dCom;
  pData.l0_ = l0;
  pData.contacts_.push_back(data);
  pData.useAngularMomentum_ = useAngMomentum;
  std::vector<double> Ts;
  Ts.push_back(timeDuration);
  ResultDataCOMTraj res = solve0step(pData, Ts, timeStep);
  return new ResultDataCOMTraj(res);
}

struct res_data_exception : std::exception {
  char const* what() const throw() { return "attributes not accessible for false resData"; }
};

void translate(const res_data_exception& e) {
  // Use the Python 'C' API to set up an exception object
  PyErr_SetString(PyExc_RuntimeError, e.what());
}

struct contact_data_exception : std::exception {
  char const* what() const throw() { return "attribute not defined yet for ContactData"; }
};

void translateContactData(const contact_data_exception& e) {
  // Use the Python 'C' API to set up an exception object
  PyErr_SetString(PyExc_RuntimeError, e.what());
}

VectorX get_xD(const ResultData& res) {
  if (res.x.size() > 0) return res.x;
  std::cout << "x is not defined" << std::endl;
  throw res_data_exception();
}

double get_costD(const ResultData& res) { return res.cost_; }

bool get_succD(const ResultData& res) { return res.success_; }

bezier_t* getC_of_t(const ResultDataCOMTraj& res) { return new bezier_t(res.c_of_t_); }

bezier_t* getDL_of_t(const ResultDataCOMTraj& res) { return new bezier_t(res.dL_of_t_); }

VectorX get_x(const ResultDataCOMTraj& res) {
  if (res.x.size() > 0) return res.x;
  std::cout << "x is not defined" << std::endl;
  throw res_data_exception();
}

double get_cost(const ResultDataCOMTraj& res) { return res.cost_; }

bool get_succ(const ResultDataCOMTraj& res) { return res.success_; }

/** BEGIN CONTACT DATA **/
centroidal_dynamics::Equilibrium* getContactPhase_(const ContactData& data) {
  if (data.contactPhase_) return data.contactPhase_;
  std::cout << "contactPhase_ is not assigned" << std::endl;
  throw contact_data_exception();
}

void setContactPhase_(ContactData& data, centroidal_dynamics::Equilibrium* eq) { data.contactPhase_ = eq; }

boost::python::tuple get_Ang(const ContactData& res) {
  if (res.ang_.size() == 0) {
    std::cout << " no angular momentum constraints assigned  " << std::endl;
    throw contact_data_exception();
  }
  return boost::python::make_tuple(res.Ang_, res.ang_);
}

boost::python::tuple get_Kin(const ContactData& res) {
  if (res.kin_.size() == 0) {
    std::cout << " no kinematic constraints assigned  " << std::endl;
    throw contact_data_exception();
  }
  return boost::python::make_tuple(res.Kin_, res.kin_);
}

void set_Kin(ContactData& res, const MatrixX3& val, const VectorX& val2) {
  if (val2.size() != val.rows()) {
    std::cout << " Kinematic inequality matrix sizes do not match  " << std::endl;
    throw contact_data_exception();
  }
  res.Kin_ = val;
  res.kin_ = val2;
}

void set_Ang(ContactData& res, const MatrixX3& val, const VectorX& val2) {
  if (val2.size() != val.rows()) {
    std::cout << " Angular inequality matrix sizes do not match  " << std::endl;
    throw contact_data_exception();
  }
  res.Ang_ = val;
  res.ang_ = val2;
}

/** END CONTACT DATA **/

/** BEGIN Constraints**/
int get_Flag(const Constraints& res) { return (int)res.flag_; }
bool get_ConstrainAcc(const Constraints& res) { return res.constraintAcceleration_; }
double get_MaxAcc(const Constraints& res) { return res.maxAcceleration_; }
double get_ReduceH(const Constraints& res) { return res.reduce_h_; }

void set_Flag(Constraints& res, const int val) { res.flag_ = (ConstraintFlag)val; }
void set_ConstrainAcc(Constraints& res, const bool val) { res.constraintAcceleration_ = val; }

void set_MaxAcc(Constraints& res, const double val) { res.maxAcceleration_ = val; }

void set_ReduceH(Constraints& res, const double val) { res.reduce_h_ = val; }

/** END Constraints **/

/** BEGIN ProblemData**/
point_t get_c0_(const ProblemData& res) { return res.c0_; }
point_t get_dc0_(const ProblemData& res) { return res.dc0_; }

point_t get_ddc0_(const ProblemData& res) { return res.ddc0_; }

point_t get_c1_(const ProblemData& res) { return res.c1_; }

point_t get_dc1_(const ProblemData& res) { return res.dc1_; }

point_t get_ddc1_(const ProblemData& res) { return res.ddc1_; }

void set_c0_(ProblemData& res, const point_t& val) { res.c0_ = val; }

void set_dc0_(ProblemData& res, const point_t& val) { res.dc0_ = val; }

void set_ddc0_(ProblemData& res, const point_t& val) { res.ddc0_ = val; }

void set_c1_(ProblemData& res, const point_t& val) { res.c1_ = val; }

void set_dc1_(ProblemData& res, const point_t& val) { res.dc1_ = val; }

void set_ddc1_(ProblemData& res, const point_t& val) { res.ddc1_ = val; }

bool get_useAngularMomentum_(const ProblemData& res) { return res.useAngularMomentum_; }
void set_useAngularMomentum_(ProblemData& res, const bool val) { res.useAngularMomentum_ = val; }

CostFunction get_costFunction_(const ProblemData& res) { return res.costFunction_; }

void set_costFunction_(ProblemData& res, const CostFunction val) { res.costFunction_ = val; }

GIWCRepresentation get_GIWC_representation_(const ProblemData& res) { return res.representation_; }

void set_GIWC_representation_(ProblemData& res, const GIWCRepresentation val) { res.representation_ = val; }

Constraints* get_constraints_(ProblemData& res) { return &res.constraints_; }

void set_constraints_(ProblemData& res, const Constraints& val) { res.constraints_ = val; }

std::vector<ContactData> get_contacts_(const ProblemData& res) { return res.contacts_; }

void addContact(ProblemData& res, const ContactData& val) { res.contacts_.push_back(ContactData(val)); }

void clearContacts(ProblemData& res) { res.contacts_.clear(); }

/** END ProblemData **/

/** BEGIN computeCOMTraj **/

ResultDataCOMTraj* computeCOMTrajPointer(const ProblemData& pData, const VectorX& Ts, const double timeStep) {
  ResultDataCOMTraj res = computeCOMTraj(pData, Ts, timeStep);
  return new ResultDataCOMTraj(res);
}

ResultDataCOMTraj* computeCOMTrajPointerChooseSolver(const ProblemData& pData, const VectorX& Ts,
                                                     const double timeStep, const solvers::SolverType solver) {
  ResultDataCOMTraj res = computeCOMTraj(pData, Ts, timeStep, solver);
  return new ResultDataCOMTraj(res);
}

/** END computeCOMTraj **/
/** BEGIN end effector **/

struct DummyPath {
  point3_t operator()(double /*u*/) const { return point3_t::Zero(); }
};

typedef std::pair<Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic>,
                  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> >
    linear_points_t;
typedef Eigen::Matrix<real, 3, Eigen::Dynamic> point_list_t;

struct MatrixVector {
  linear_points_t res;
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> A() { return res.first; }
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> b() { return res.second; }
};

ResultDataCOMTraj* computeEndEffector(const ProblemData& pData, const double time) {
  ResultDataCOMTraj res = solveEndEffector<DummyPath>(pData, DummyPath(), time, 0);
  return new ResultDataCOMTraj(res);
}

MatrixVector* computeEndEffectorConstraintsPython(const ProblemData& pData, const double time) {
  std::vector<bezier_t::point_t> pi = computeConstantWaypoints(pData, time);
  MatrixVector* res = new MatrixVector();
  res->res = computeEndEffectorConstraints(pData, time, pi);
  return res;
}

MatrixVector* computeEndEffectorVelocityCostPython(const ProblemData& pData, const double time) {
  std::vector<bezier_t::point_t> pi = computeConstantWaypoints(pData, time);
  MatrixVector* res = new MatrixVector();
  res->res = computeVelocityCost(pData, time, pi);
  return res;
}

MatrixVector* computeEndEffectorDistanceCostPython(const ProblemData& pData, const double time, const int numPoints,
                                                   point_list_t pts_l) {
  std::vector<bezier_t::point_t> pi = computeConstantWaypoints(pData, time);
  // transform the matrice 3xN in a std::vector<point3_t> of size N :
  std::vector<point3_t> pts_path;
  for (size_t c = 0; c < pts_l.cols(); ++c) {
    pts_path.push_back(pts_l.block<3, 1>(0, c));
  }
  MatrixVector* res = new MatrixVector();
  res->res = computeDistanceCostFunction(numPoints, pData, time, pts_path);
  return res;
}

Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> computeEndEffectorConstantWaypoints(const ProblemData& pData,
                                                                                        const double time) {
  std::vector<bezier_t::point_t> pi = computeConstantWaypoints(pData, time);
  Eigen::Matrix<real, Eigen::Dynamic, Eigen::Dynamic> res(3, pi.size());
  int col = 0;
  for (std::vector<bezier_t::point_t>::const_iterator cit = pi.begin(); cit != pi.end(); ++cit, ++col)
    res.block<3, 1>(0, col) = *cit;
  return res;
}

/** END end effector **/

BOOST_PYTHON_MODULE(hpp_bezier_com_traj) {
  using namespace boost::python;
  register_exception_translator<res_data_exception>(&translate);
  register_exception_translator<contact_data_exception>(&translateContactData);
  /** BEGIN eigenpy init**/
  eigenpy::enableEigenPy();

  eigenpy::enableEigenPySpecific<point_t, point_t>();
  eigenpy::enableEigenPySpecific<Vector3, Vector3>();
  eigenpy::enableEigenPySpecific<VectorX, VectorX>();
  eigenpy::enableEigenPySpecific<MatrixX3, MatrixX3>();
  eigenpy::enableEigenPySpecific<MatrixX3, MatrixX3>();

  /** END eigenpy init**/
  /*eigenpy::exposeAngleAxis();
  eigenpy::exposeQuaternion();*/

  class_<ResultDataCOMTraj>("ResultDataCOMTraj", init<>())
      .add_property("c_of_t", make_function(&getC_of_t, return_value_policy<manage_new_object>()))
      .add_property("dL_of_t", make_function(&getDL_of_t, return_value_policy<manage_new_object>()))
      .add_property("success", &get_succ)
      .add_property("cost", &get_cost)
      .add_property("x", &get_x);

  class_<ResultData>("ResultData", init<>())
      .add_property("success", &get_succD)
      .add_property("cost", &get_costD)
      .add_property("x", &get_xD);

  class_<ContactData>("ContactData", init<centroidal_dynamics::Equilibrium*>()[with_custodian_and_ward<1, 2>()])
      .add_property("contactPhase_",
                    make_function(&getContactPhase_, return_value_policy<reference_existing_object>()),
                    &setContactPhase_)
      .add_property("Kin_", &get_Kin)
      .add_property("Ang_", &get_Ang)
      .def("setKinematicConstraints", &set_Kin)
      .def("setAngularConstraints", &set_Ang);

  class_<ProblemData>("ProblemData", init<>())
      .add_property("c0_", &get_c0_, &set_c0_)
      .add_property("dc0_", &get_dc0_, &set_dc0_)
      .add_property("ddc0_", &get_ddc0_, &set_ddc0_)
      .add_property("c1_", &get_c1_, &set_c1_)
      .add_property("dc1_", &get_dc1_, &set_dc1_)
      .add_property("ddc1_", &get_ddc1_, &set_ddc1_)
      .add_property("useAngularMomentum_", &get_useAngularMomentum_, &set_useAngularMomentum_)
      .add_property("costFunction_", &get_costFunction_, &set_costFunction_)
      .add_property("GIWCrepresentation_", &get_GIWC_representation_, &set_GIWC_representation_)
      .add_property("constraints_", make_function(&get_constraints_, return_value_policy<reference_existing_object>()),
                    &set_constraints_)
      .def("clearContacts", clearContacts)
      .def("addContact", addContact);

  class_<Constraints>("Constraints", init<>())
      .add_property("flag_", &get_Flag, &set_Flag)
      .add_property("constrainAcceleration_", &get_ConstrainAcc, &set_ConstrainAcc)
      .add_property("maxAcceleration_", &get_MaxAcc, &set_MaxAcc)
      .add_property("reduce_h_", &get_ReduceH, &set_ReduceH);

  enum_<CostFunction>("CostFunction")
      .value("ACCELERATION", ACCELERATION)
      .value("DISTANCE_TRAVELED", DISTANCE_TRAVELED)
      .value("TARGET_END_VELOCITY", TARGET_END_VELOCITY)
      .value("UNKNOWN_COST", UNKNOWN_COST)
      .export_values();

  enum_<GIWCRepresentation>("GIWCRepresentation")
      .value("DOUBLE_DESCRIPTION", DOUBLE_DESCRIPTION)
      .value("FORCE", FORCE)
      .value("UNKNOWN_REPRESENTATION", UNKNOWN_REPRESENTATION)
      .export_values();

  enum_<solvers::SolverType>("SolverType")
      .value("SOLVER_QUADPROG", solvers::SOLVER_QUADPROG)
  //.value("SOLVER_QUADPROG_SPARSE", solvers::SOLVER_QUADPROG_SPARSE)
#ifdef USE_GLPK_SOLVER
      .value("SOLVER_GLPK", solvers::SOLVER_GLPK)
#endif
      .export_values();

  enum_<ConstraintFlag>("ConstraintFlag")
      .value("INIT_POS", INIT_POS)
      .value("INIT_VEL", INIT_VEL)
      .value("INIT_ACC", INIT_ACC)
      .value("INIT_JERK", INIT_JERK)
      .value("END_POS", END_POS)
      .value("END_VEL", END_VEL)
      .value("END_ACC", END_ACC)
      .value("END_JERK", END_JERK)
      .value("ONE_FREE_VAR", ONE_FREE_VAR)
      .value("TWO_FREE_VAR", TWO_FREE_VAR)
      .value("THREE_FREE_VAR", THREE_FREE_VAR)
      .value("FOUR_FREE_VAR", FOUR_FREE_VAR)
      .value("FIVE_FREE_VAR", FIVE_FREE_VAR)
      .value("UNKNOWN", UNKNOWN)
      .export_values();

  class_<MatrixVector>("MatrixVector", no_init)
      .def_readonly("A", &MatrixVector::A)
      .def_readonly("b", &MatrixVector::b);

  def("zeroStepCapturability", &zeroStepCapturability, return_value_policy<manage_new_object>());
  def("zeroStepCapturability", &zeroStepCapturabilityWithKinConstraints, return_value_policy<manage_new_object>());
  def("computeCOMTraj", &computeCOMTrajPointer, return_value_policy<manage_new_object>());
  def("computeCOMTraj", &computeCOMTrajPointerChooseSolver, return_value_policy<manage_new_object>());
  def("computeEndEffector", &computeEndEffector, return_value_policy<manage_new_object>());
  def("computeEndEffectorConstraints", &computeEndEffectorConstraintsPython, return_value_policy<manage_new_object>());
  def("computeEndEffectorVelocityCost", &computeEndEffectorVelocityCostPython,
      return_value_policy<manage_new_object>());
  def("computeEndEffectorDistanceCost", &computeEndEffectorDistanceCostPython,
      return_value_policy<manage_new_object>());
  def("computeEndEffectorConstantWaypoints", &computeEndEffectorConstantWaypoints);
}

}  // namespace bezier_com_traj
