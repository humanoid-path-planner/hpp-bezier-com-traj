#ifndef TEST_HELPER_HH
#define TEST_HELPER_HH

#include <bezier-com-traj/solve.hh>
#include <bezier-com-traj/common_solve_methods.hh>
#include <centroidal-dynamics-lib/centroidal_dynamics.hh>

using bezier_com_traj::MatrixXX;
using bezier_com_traj::MatrixX3;
using bezier_com_traj::Matrix3;
using bezier_com_traj::VectorX;
using bezier_com_traj::Vector3;

#define MASS 50.
#define MU 0.5
#define LX 0.2172        // contact surface size in x direction
#define LY 0.138         // contact surface size in y direction
#define KIN_X_MIN -0.3
#define KIN_X_MAX  1
#define KIN_Y_MIN -0.5
#define KIN_Y_MAX  0.5

typedef std::pair<MatrixX3,VectorX> ConstraintsPair;

std::pair<MatrixX3, MatrixX3> generateKinematicsConstraints(){
    // generate a simple polytgone  : faces aligned along x,y,z axis
    // size : x [-0.5,0.5] ; y = [-0.3,1] ; z [0,0.8]
    MatrixX3 N(6,3);
    MatrixX3 V(6,3);

    N.block<1,3>(0,0) = Vector3(-1,0,0);
    V.block<1,3>(0,0) = Vector3(KIN_X_MIN,KIN_Y_MAX,0);
    N.block<1,3>(1,0) = Vector3(0,1,0);
    V.block<1,3>(1,0) = Vector3(KIN_X_MIN,KIN_Y_MAX,0);
    N.block<1,3>(2,0) = Vector3(1,0,0);
    V.block<1,3>(2,0) = Vector3(KIN_X_MAX,KIN_Y_MIN,0);
    N.block<1,3>(3,0) = Vector3(0,-1,0);
    V.block<1,3>(3,0) = Vector3(KIN_X_MAX,KIN_Y_MIN,0);
    N.block<1,3>(4,0) = Vector3(0,0,-1);
    V.block<1,3>(4,0) = Vector3(0,0,0);
    N.block<1,3>(5,0) = Vector3(0,0,1);
    V.block<1,3>(5,0) = Vector3(0,0,0.8);

    return std::make_pair(N,V);
}


std::pair<MatrixXX, VectorX> generateKinematicsConstraints(Matrix3 endEffRotation, Vector3 endEffTranslation){

    std::pair<MatrixX3, MatrixX3> NV = generateKinematicsConstraints();
    MatrixX3 N = NV.first;
    MatrixX3 V = NV.second;
    size_t numFaces = N.rows();
    MatrixX3 A(numFaces,3);
    VectorX b(numFaces);
    VectorX n,v;

    for(size_t i = 0 ; i < numFaces ; ++i){
        n = endEffRotation * (N.block<1,3>(i,0).transpose());
        v = endEffRotation * (V.block<1,3>(i,0).transpose()) + endEffTranslation;
        A.block<1,3>(i,0) = n;
        b[i] = v.dot(n);
    }

    return std::make_pair(A,b);
}

std::pair<MatrixX3, MatrixX3> computeRectangularContacts(MatrixX3 normals, MatrixX3 positions, double size_X,double size_Y){
    // TODO : consider normal != z (see code in rbprm :: stability.cc (or add it as dependency ?)

    BOOST_CHECK(normals.rows() == positions.rows());
    MatrixX3 rec_normals(normals.rows()*4,3);
    MatrixX3 rec_positions(normals.rows()*4,3);

    double lx = size_X/2.;
    double ly = size_Y/2.;
    MatrixX3 p(4,3);
    p << lx,  ly, 0,
         lx, -ly, 0,
        -lx, -ly, 0,
        -lx,  ly, 0;


    for (size_t ic = 0 ; ic < normals.rows() ; ++ic){
        for (size_t i = 0 ; i < 4 ; ++i){
            rec_normals.block<1,3>(ic*4+i,0) = normals.block<1,3>(ic,0);
            rec_positions.block<1,3>(ic*4+i,0) = positions.block<1,3>(ic,0) + p.block<1,3>(i,0);
        }
    }
    return std::make_pair(rec_normals,rec_positions);
}

centroidal_dynamics::Equilibrium ComputeContactCone(MatrixX3 normals, MatrixX3 positions){
    centroidal_dynamics::Equilibrium contactCone("test-quasiStatic", MASS,4,centroidal_dynamics::SOLVER_LP_QPOASES,true,10,false);
    centroidal_dynamics::EquilibriumAlgorithm alg = centroidal_dynamics::EQUILIBRIUM_ALGORITHM_PP;
    contactCone.setNewContacts(positions,normals,MU,alg);
    return contactCone;
}

std::pair<MatrixXX, VectorX> generateStabilityConstraints(centroidal_dynamics::Equilibrium contactPhase,Vector3 acc = Vector3::Zero()){
    const Vector3& g = contactPhase.m_gravity;
    const Matrix3 gSkew = bezier_com_traj::skew(g);
    const Matrix3 accSkew = bezier_com_traj::skew(acc);
    // compute GIWC
    centroidal_dynamics::MatrixXX Hrow;
    VectorX h;
    contactPhase.getPolytopeInequalities(Hrow,h);
    MatrixXX H = -Hrow;
    H.rowwise().normalize();
    int dimH = (int)(H.rows());
    MatrixXX mH = contactPhase.m_mass * H;
    // constraints : mH[:,3:6] g^  x <= h + mH[:,0:3]g
    // A = mH g^
    // b = h + mHg
    MatrixX3 A = mH.block(0,3,dimH,3) * (gSkew - accSkew);
    VectorX b = h+mH.block(0,0,dimH,3)*(g - acc);
    return std::make_pair(A,b);
}

std::pair<MatrixXX, VectorX> generateStabilityConstraints(MatrixX3 normals, MatrixX3 positions,Vector3 acc = Vector3::Zero()){
    std::pair<MatrixX3, MatrixX3> contacts = computeRectangularContacts(normals,positions,LX,LY);
    centroidal_dynamics::Equilibrium contactPhase = ComputeContactCone(contacts.first,contacts.second);
    return generateStabilityConstraints(contactPhase,acc);
}

std::pair<Matrix3, Vector3> computeCost(){
    Matrix3 H = Matrix3::Identity();
    Vector3 g = Vector3::Zero();
    return std::make_pair(H,g);
}

std::pair<MatrixX3,VectorX> generateConstraints(MatrixX3 normals, MatrixX3 positions,Matrix3 endEffRotation, Vector3 endEffTranslation){
    std::pair<MatrixX3,VectorX> Ab = generateKinematicsConstraints(endEffRotation,endEffTranslation);
    std::pair<MatrixX3,VectorX> Cd = generateStabilityConstraints(normals,positions);
    size_t numIneq = Ab.first.rows() + Cd.first.rows();
    MatrixXX M(numIneq,3);
    VectorX  n(numIneq);
    M.block(0,0,Ab.first.rows(),3) = Ab.first;
    M.block(Ab.first.rows(),0,Cd.first.rows(),3) = Cd.first;
    n.segment(0,Ab.first.rows()) = Ab.second;
    n.segment(Ab.first.rows(),Cd.first.rows()) = Cd.second;
    return std::make_pair(M,n);
}


double fRandom(double fMin, double fMax)
{
    double f = (double)std::rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}



ConstraintsPair stackConstraints(const ConstraintsPair& Ab,const ConstraintsPair& Cd){
    size_t numIneq = Ab.first.rows() + Cd.first.rows();
    MatrixX3 M(numIneq,3);
    VectorX  n(numIneq);
    M.block(0,0,Ab.first.rows(),3) = Ab.first;
    M.block(Ab.first.rows(),0,Cd.first.rows(),3) = Cd.first;
    n.segment(0,Ab.first.rows()) = Ab.second;
    n.segment(Ab.first.rows(),Cd.first.rows()) = Cd.second;
    return std::make_pair(M,n);
}

#endif // TEST_HELPER_HH
