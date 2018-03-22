
namespace bezier_com_traj
{

template <typename Point>
std::vector< std::pair<double,Point> > computeDiscretizedWaypoints
    (const ProblemData& pData, double T,const T_time& timeArray)
{
    typedef std::pair<double,Point> coefs_t;
    //int numStep = int(T / timeStep);
    std::vector<coefs_t> wps;
    std::vector<Point> pi = computeConstantWaypoints(pData,T);
    // evaluate curve work with normalized time !
    double t;
    for (std::size_t i = 0 ; i<timeArray.size() ; ++i )
    {
        t = std::min(timeArray[i] / T, 1.);
        wps.push_back(evaluateCurveAtTime(pData,pi,t));
    }
    return wps;
}


template <typename Point>
std::vector< std::pair<double,Point> > computeDiscretizedAccelerationWaypoints
    (const ProblemData& pData, double T,const T_time& timeArray)
{
    typedef std::pair<double,Point> coefs_t;
    std::vector<coefs_t> wps;
    std::vector<Point> pi = computeConstantWaypoints(pData,T);
    double t;
    for (std::size_t i = 0 ; i<timeArray.size() ; ++i ){
        t = timeArray[i] / T;
        if(t>1)
            t=1.;
        wps.push_back(evaluateAccelerationCurveAtTime(pData,pi,T,t));
    }
    return wps;
}

} // end namespace bezier_com_traj


