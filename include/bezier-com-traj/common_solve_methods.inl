
namespace bezier_com_traj
{

template <typename Point>
std::vector< std::pair<double,Point> > computeDiscretizedWaypoints
    (const ProblemData& pData, double T,const T_time& timeArray)
{
    typedef std::pair<double,Point> coefs_t;
    std::vector<coefs_t> wps;
    std::vector<Point> pi = computeConstantWaypoints(pData,T);
    // evaluate curve work with normalized time !
    for (CIT_time cit = timeArray.begin(); cit != timeArray.end(); ++cit)
    {
        double t = std::min(cit->first / T, 1.);
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
    // evaluate curve work with normalized time !
    for (CIT_time cit = timeArray.begin(); cit != timeArray.end(); ++cit)
    {
        double t = std::min(cit->first / T, 1.);
        wps.push_back(evaluateAccelerationCurveAtTime(pData,pi,T,t));
    }
    return wps;
}

} // end namespace bezier_com_traj


