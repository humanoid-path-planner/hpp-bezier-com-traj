/*
 * Copyright 2018, LAAS-CNRS
 * Author: Steve Tonneau
 */

#ifndef BEZIER_COM_TRAJ_FLAGS_H
#define BEZIER_COM_TRAJ_FLAGS_H

#include <bezier-com-traj/config.hh>

namespace bezier_com_traj
{
    enum BEZIER_COM_TRAJ_DLLAPI CostFunction{
        ACCELERATION        = 0x00001,
        DISTANCE_TRAVELED   = 0x00002,
        TARGET_END_VELOCITY = 0x00004,
        UNKNOWN_COST        = 0x00008
      };

    enum BEZIER_COM_TRAJ_DLLAPI ConstraintFlag{
        INIT_POS = 0x00001,
        INIT_VEL = 0x00002,
        INIT_ACC = 0x00004,
        END_POS  = 0x00008,
        END_VEL  = 0x00010,
        END_ACC  = 0x00020,
        UNKNOWN  = 0x00040
      };

    inline ConstraintFlag operator~(ConstraintFlag a)
    {return static_cast<ConstraintFlag>(~static_cast<const int>(a));}

    inline ConstraintFlag operator|(ConstraintFlag a, ConstraintFlag b)
    {return static_cast<ConstraintFlag>(static_cast<const int>(a) | static_cast<const int>(b));}

    inline ConstraintFlag operator&(ConstraintFlag a, ConstraintFlag b)
    {return static_cast<ConstraintFlag>(static_cast<const int>(a) & static_cast<const int>(b));}

    inline ConstraintFlag operator^(ConstraintFlag a, ConstraintFlag b)
    {return static_cast<ConstraintFlag>(static_cast<const int>(a) ^ static_cast<const int>(b));}

    inline ConstraintFlag& operator|=(ConstraintFlag& a, ConstraintFlag b)
    {return (ConstraintFlag&)((int&)(a) |= static_cast<const int>(b));}

    inline ConstraintFlag& operator&=(ConstraintFlag& a, ConstraintFlag b)
    {return (ConstraintFlag&)((int&)(a) &= static_cast<const int>(b));}

    inline ConstraintFlag& operator^=(ConstraintFlag& a, ConstraintFlag b)
    {return (ConstraintFlag&)((int&)(a) ^= static_cast<const int>(b));}

} // end namespace bezier_com_traj

#endif
