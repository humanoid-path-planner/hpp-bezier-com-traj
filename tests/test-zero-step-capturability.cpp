/*
 * Copyright 2018, LAAS-CNRS
 * Author: Steve Tonneau
 */

#include <vector>
#include <iostream>
#include <centroidal-dynamics-lib/centroidal_dynamics.hh>
#include <bezier-com-traj/solve.hh>
#include <math.h>

using namespace centroidal_dynamics;
using namespace Eigen;
using namespace std;

#define PERF_PP "Polytope Projection"
#define PERF_LP_PREPARATION "Computation of GIWC generators"
#define PERF_LP_COIN "Compute Equilibrium Robustness with LP coin"
#define PERF_LP_OASES "Compute Equilibrium Robustness with LP oases"
#define PERF_LP2_COIN "Compute Equilibrium Robustness with LP2 coin"
#define PERF_LP2_OASES "Compute Equilibrium Robustness with LP2 oases"
#define PERF_DLP_COIN "Compute Equilibrium Robustness with DLP coin"
#define PERF_DLP_OASES "Compute Equilibrium Robustness with DLP oases"

#define EPS 1e-3  // required precision




void generateContacts(unsigned int N_CONTACTS, double MIN_CONTACT_DISTANCE, double LX, double LY,
                      RVector3 &CONTACT_POINT_LOWER_BOUNDS,
                      RVector3 &CONTACT_POINT_UPPER_BOUNDS,
                      RVector3 &RPY_LOWER_BOUNDS,
                      RVector3 &RPY_UPPER_BOUNDS,
                      MatrixX3& p, MatrixX3& N)
{
  MatrixXX contact_pos = MatrixXX::Zero(N_CONTACTS, 3);
  MatrixXX contact_rpy = MatrixXX::Zero(N_CONTACTS, 3);
  p.setZero(4*N_CONTACTS,3); // contact points
  N.setZero(4*N_CONTACTS,3); // contact normals

  // Generate contact positions and orientations
  bool collision;
  for(unsigned int i=0; i<N_CONTACTS; i++)
  {
    while(true) // generate contact position
    {
      uniform(CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS, contact_pos.row(i));
      if(i==0)
        break;
      collision = false;
      for(unsigned int j=0; j<i-1; j++)
        if((contact_pos.row(i)-contact_pos.row(j)).norm() < MIN_CONTACT_DISTANCE)
          collision = true;
      if(collision==false)
        break;
    }

//     generate contact orientation
    uniform(RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS, contact_rpy.row(i));
    generate_rectangle_contacts(LX, LY, contact_pos.row(i), contact_rpy.row(i),
                                p.middleRows<4>(i*4), N.middleRows<4>(i*4));
//    printf("Contact surface %d position (%.3f,%.3f,%.3f) ", i, contact_pos(i,0), contact_pos(i,1), contact_pos(i,2));
//    printf("Orientation (%.3f,%.3f,%.3f)\n", contact_rpy(i,0), contact_rpy(i,1), contact_rpy(i,2));
  }

//  for(int i=0; i<p.rows(); i++)
//  {
//    printf("Contact point %d position (%.3f,%.3f,%.3f) ", i, p(i,0), p(i,1), p(i,2));
//    printf("Normal (%.3f,%.3f,%.3f)\n", N(i,0), N(i,1), N(i,2));
//  }
}

double fRandom(double fMin, double fMax)
{
    double f = (double)std::rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

bool findStaticEqCOMPos(Equilibrium * eq, const RVector3& com_LB, const RVector3& com_UB, Vector3& c0)
{
    int trials = 0;
    while (trials < 1000)
    {
        ++trials;
        double x = fRandom(com_LB(0), com_UB(0));
        double y = fRandom(com_LB(1), com_UB(1));
        double z = fRandom(com_LB(2), com_UB(2));
        c0 << x , y, z;
        bool isStable(false);
        eq->checkRobustEquilibrium(c0,isStable);
        if(isStable)
        {
            return true;
        }
    }
    return false;
}

Vector6 computew(const Equilibrium* eq, const bezier_com_traj::Vector3 c, const Vector3 ddc,
                 const Vector3 dL)
{
    Vector6 w; Vector3 w1;
    w1 = eq->m_mass * (ddc -eq->m_gravity);
    w.head(3) = w1;
    w.tail(3) = c.cross(w1) + dL;
    return w;
}

bool checkTrajectory(const Vector3& c0, const std::string& solver, const Equilibrium* eq, const bezier_com_traj::ResultDataCOMTraj& resData, const double T, const int num_steps = 100)
{
    // retrieve H
    centroidal_dynamics::MatrixXX Hrow; VectorX h;
    eq->getPolytopeInequalities(Hrow,h);
    MatrixXX H = -Hrow;
    H.rowwise().normalize();
    const bezier_com_traj::bezier_t c_of_t  = resData.c_of_t_;
    const bezier_com_traj::bezier_t dL_of_t = resData.dL_of_t_;
    bezier_com_traj::bezier_t ddc_of_t = c_of_t.compute_derivate(2);
    for (int i = 0; i < num_steps; ++i)
    {
        double dt = double(i) / float(num_steps) * T;
        c_of_t(dt);
        dL_of_t(dt);
        Vector6 w = computew(eq, c_of_t(dt),ddc_of_t(dt),dL_of_t(dt));
        VectorX res = H*w;
        for(long j=0; j<res.size(); j++)
          if(res(j)>0.0001)
          {
            std::cout << "check trajectory failed for solver " << solver  << " " <<  res(j)  << "; iteration " << dt<< "; numsteps " << num_steps << " c_of_t(dt) " << (c_of_t)(dt).transpose() << " dc_of_t(dt) " << c_of_t.derivate(dt,1).transpose() << " dL_of_t(dt) " << (dL_of_t)(dt).transpose()  << std::endl;
            std::cout << "c0 " << c0.transpose() << std::endl;
            std::cout << "H row " << H.row(j) << std::endl;
            std::cout << "w " << w << std::endl;
            std::cout << "j " << j << std::endl;
            //throw;
            std::cout << "mult ddc " << ddc_of_t.mult_T_ << std::endl;
            std::cout << "mult c " << c_of_t.mult_T_ << std::endl;
            std::cout << "T c" << c_of_t.T_ << std::endl;
            std::cout << "T " << T << std::endl;
            std::cout << "T ddc " << ddc_of_t.T_ << std::endl;
            std::cout << "1 / T * T " <<  1/ (T * T) << std::endl;
            return false;
          }
    }
    return true;
}

int main()
{
  srand((unsigned int)time(NULL));
  RVector3 CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS;
  RVector3 RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS;

  /************************************** USER PARAMETERS *******************************/
  unsigned int N_TESTS = 500;
  double mass = 55.0;
  double mu = 0.3;  // friction coefficient
  unsigned int generatorsPerContact = 4;
  unsigned int N_CONTACTS = 2;
  double MIN_CONTACT_DISTANCE = 0.3;
  double LX = 0.5*0.2172;        // half contact surface size in x direction
  double LY = 0.5*0.138;         // half contact surface size in y direction
  CONTACT_POINT_LOWER_BOUNDS << 0.0,  0.0,  0.0;
  CONTACT_POINT_UPPER_BOUNDS << 0.5,  0.5,  0.5;
  double gamma = atan(mu);   // half friction cone angle
  RPY_LOWER_BOUNDS << -2*gamma, -2*gamma, -M_PI;
  RPY_UPPER_BOUNDS << +2*gamma, +2*gamma, +M_PI;
  double X_MARG = 0.07;
  double Y_MARG = 0.07;
  /************************************ END USER PARAMETERS *****************************/

  MatrixX3 p, N;
  RVector3 com_LB, com_UB;
  Equilibrium solver_PP ("PP", mass, generatorsPerContact, SOLVER_LP_QPOASES,false,10,false);
  int succContinuous = 0, succDiscretize = 0, succdL = 0, succDiscretizedL = 0 , succKin = 0, succdLKin = 0, succdLAng = 0;
  int numSol = 0;
  for(unsigned n_test=0; n_test<N_TESTS; n_test++)
  {
    generateContacts(N_CONTACTS, MIN_CONTACT_DISTANCE, LX, LY,
                     CONTACT_POINT_LOWER_BOUNDS, CONTACT_POINT_UPPER_BOUNDS,
                     RPY_LOWER_BOUNDS, RPY_UPPER_BOUNDS, p, N);

    // compute upper and lower bounds of com positions to test
    com_LB(0) = p.col(0).minCoeff()-X_MARG;
    com_UB(0) = p.col(0).maxCoeff()+X_MARG;
    com_LB(1) = p.col(1).minCoeff()-Y_MARG;
    com_UB(1) = p.col(1).maxCoeff()+Y_MARG;
    com_LB(2) = p.col(2).minCoeff()+0.3;
    com_UB(2) = p.col(2).maxCoeff()+1.5;

    Vector3 c0;
    solver_PP.setNewContacts(p,N,mu,EQUILIBRIUM_ALGORITHM_PP);
    const double DISCRETIZATION_STEP = 0.05;
    if(findStaticEqCOMPos(&solver_PP, com_LB, com_UB,c0))
    {
        numSol +=1;
        for(int j = 0; j < 100; ++j)
        {
            bezier_com_traj::ContactData data;
            data.contactPhase_ = &solver_PP;
            bezier_com_traj::ProblemData pData;
            pData.c0_ = c0;
            pData.dc0_ << fRandom(-1.,1.) , fRandom(-1.,1.) , fRandom(-1.,1.);
            //pData.dc0_ *= 1;
            pData.l0_ = Vector3(0.,0.,0.);
            pData.contacts_.push_back(data);            
            double T;
            for (int k = -1; k < 2; ++k)
            //for (int k = 0; k < 1; ++k)
            {
                pData.useAngularMomentum_ = false;
                pData.contacts_.front().kin_ = VectorX::Zero(0);
                bool succCont = false, succDisc = false, succdLbool= false, succDiscdL = false ,
                        succKinBool = false, succdLKinBool = false, succdLAngBool = false;
                T = 1. + 0.3 * k;
                std::vector<double> Ts;
                Ts.push_back(T);
                bezier_com_traj::ResultDataCOMTraj rData = bezier_com_traj::solve0step(pData,Ts);
                if(rData.success_)
                {                    
                    assert ((rData.c_of_t_(0.) - pData.c0_).norm() < 0.0001);
                    succCont = true;
                    succContinuous += 1;
                    std::string solverName("continuous");
                    checkTrajectory(c0, solverName, &solver_PP,rData, T);
                }
                else
                {
                    succCont = false;
                }
               // try discretize
                bezier_com_traj::ResultDataCOMTraj rData0 = bezier_com_traj::solve0step(pData,Ts,DISCRETIZATION_STEP);
                if(rData0.success_)
                {
                    assert ((rData0.c_of_t_(0.) - c0).norm() < 0.0001);
                    succDisc = true;
                    succDiscretize += 1;
                    std::string solverName("discretize");
                    checkTrajectory(c0, solverName, &solver_PP,rData0,T, int(T / DISCRETIZATION_STEP));
                }
                else
                {
                    if(succCont)
                        std::cout << "error: Solver discretize failed while a solution was found for the continuous case" << std::endl;
                    succDisc = false;
                }

                pData.useAngularMomentum_ = true;
                bezier_com_traj::ResultDataCOMTraj rData1 = bezier_com_traj::solve0step(pData,Ts);
                if(rData1.success_)
                {
                    assert ((rData1.c_of_t_(0.) - c0).norm() < 0.0001);
                    succdLbool = true;
                    succdL += 1;
                    std::string solverName("continuous AngularMomentum");
                    checkTrajectory(c0, solverName, &solver_PP,rData1, T);
                }
                else
                {
                    //if(succCont)
                    //    std::cout << "error: Solver Ang momentum failed while a solution was found without angular momentum" << std::endl;
                    succDisc = false;
                }

                bezier_com_traj::ResultDataCOMTraj rData2 = bezier_com_traj::solve0step(pData,Ts,DISCRETIZATION_STEP);
                if(rData2.success_)
                {
                    assert ((rData2.c_of_t_(0.) - c0).norm() < 0.0001);
                    succDiscdL = true;
                    succDiscretizedL += 1;
                    std::string solverName("discretize AngularMomentum");
                    checkTrajectory(c0, solverName, &solver_PP,rData2,T, int(T / DISCRETIZATION_STEP ));
                }
                else
                {
                    //if(succCont || succDisc ||succdLbool)
                    //    std::cout << "error: Solver discretize with angular momentum failed while a solution was found for another case" << std::endl;
                    succDiscdL = false;
                }
                pData.contacts_.front().Kin_ = Eigen::Matrix3d::Identity();
                pData.contacts_.front().kin_ = Vector3::Constant(3,0.5);
                bezier_com_traj::ResultDataCOMTraj rData3 = bezier_com_traj::solve0step(pData,Ts);
                pData.contacts_.front().Ang_ = Eigen::Matrix3d::Identity();
                pData.contacts_.front().ang_ = Vector3::Constant(3,0.1);
                if(rData3.success_)
                {
                    assert ((rData3.c_of_t_(0.) - c0).norm() < 0.0001);
                    succdLKinBool = true;
                    succdLKin += 1;
                    std::string solverName("kinematic and angular constraint");
                    checkTrajectory(c0, solverName, &solver_PP,rData3,T);
                }

                pData.contacts_.front().kin_ = VectorX::Zero(0);
                bezier_com_traj::ResultDataCOMTraj rData31 = bezier_com_traj::solve0step(pData,Ts);
                if(rData31.success_)
                {
                    assert ((rData31.c_of_t_(0.) - c0).norm() < 0.0001);
                    succdLAngBool = true;
                    succdLAng += 1;
                    std::string solverName("angular constraint");
                    checkTrajectory(c0, solverName, &solver_PP,rData31,T);
                }

                pData.useAngularMomentum_ = false;
                pData.contacts_.front().ang_ = VectorX::Zero(0);
                pData.contacts_.front().kin_ = Vector3::Constant(3,0.5);
                bezier_com_traj::ResultDataCOMTraj rData4 = bezier_com_traj::solve0step(pData,Ts);
                if(rData4.success_)
                {
                    assert ((rData4.c_of_t_(0.) - c0).norm() < 0.0001);
                    succKinBool = true;
                    succKin += 1;
                    std::string solverName("kinematic constraint");
                    checkTrajectory(c0, solverName, &solver_PP,rData4,T);
                }
            }
        }
    }
  }

  std::cout << "sucesses continunous" << succContinuous << std::endl;
  std::cout << "sucesses continuous with Kinematic constraints " << succKin << std::endl;
  std::cout << "sucesses discretize " << succDiscretize << std::endl;
  std::cout << "sucesses continunous with angular momentum" << succdL << std::endl;
  std::cout << "sucesses discretize with angular momentum" << succDiscretizedL << std::endl;
  std::cout << "sucesses continunous with angular momentum and angular constraints" << succdLAng << std::endl;
  std::cout << "sucesses continunous with angular momentum and kinematic and angular constraints" << succdLKin << std::endl;
  std::cout << "sucesses discretize with angular momentum" << succDiscretizedL << std::endl;
  std::cout << "eq static point found " << numSol << std::endl;


  return 0;
}
