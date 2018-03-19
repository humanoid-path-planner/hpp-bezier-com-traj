/*
 * Copyright 2017, LAAS-CNRS
 * Author: Steve Tonneau
 */

#include <bezier-com-traj/utils.hh>
namespace bezier_com_traj
{

Matrix3 skew(point_t_tC x)
{
    Matrix3 res = Matrix3::Zero();
    res(0,1) = - x(2); res(0,2) =   x(1);
    res(1,0) =   x(2); res(1,2) = - x(0);
    res(2,0) = - x(1); res(2,1) =   x(0);
    return res;
}

std::vector<spline::Bern<double> > ComputeBersteinPolynoms(const unsigned int degree)
{
    std::vector<spline::Bern<double> > res;
    for (unsigned int i =0; i <= (unsigned int)degree; ++i)
        res.push_back(spline::Bern<double>(degree,i));
    return res;
}

void printQHullFile(const std::pair<MatrixXX, VectorX>& Ab,VectorX intPoint,const std::string& fileName,bool clipZ){
     std::ofstream file;
     using std::endl;
     std::string path("/local/fernbac/bench_iros18/constraints_obj/");
     path.append(fileName);
     file.open(path.c_str(),std::ios::out | std::ios::trunc);
     file<<"3 1"<<endl;
     file<<"\t "<<intPoint[0]<<"\t"<<intPoint[1]<<"\t"<<intPoint[2]<<endl;
     file<<"4"<<endl;
     clipZ ? file<<Ab.first.rows()+2<<endl : file<<Ab.first.rows()<<endl;
     for(int i = 0 ; i < Ab.first.rows() ; ++i){
         file<<"\t"<<Ab.first(i,0)<<"\t"<<Ab.first(i,1)<<"\t"<<Ab.first(i,2)<<"\t"<<-Ab.second[i]-0.001<<endl;
     }
     if(clipZ){
         file<<"\t"<<0<<"\t"<<0<<"\t"<<1.<<"\t"<<-3.<<endl;
         file<<"\t"<<0<<"\t"<<0<<"\t"<<-1.<<"\t"<<-1.<<endl;
     }
     file.close();
}


} // namespace bezier_com_traj
