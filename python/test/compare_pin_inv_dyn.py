"""
Created on Thu Sep  1 16:54:39 2016

@author: stonneau
"""

from math import atan, pi

import numpy as np
from centroidal_dynamics import Equilibrium, EquilibriumAlgorithm
from curves import bezier
from numpy import array, asarray, asmatrix, matrix, zeros
from numpy import cross as X
from pinocchio_inv_dyn.multi_contact.bezier.bezier_0_step_capturability import (
    BezierZeroStepCapturability,
    compute_CWC,
)
from pinocchio_inv_dyn.multi_contact.utils import (
    check_static_eq,
    find_static_equilibrium_com,
    generate_contacts,
)

__EPS = 1e-5
np.set_printoptions(precision=2, suppress=True, linewidth=100)


#####################
# EQUILIBRIUM CHECK #
#####################
def compute_w(
    c, ddc, dL=array([0.0, 0.0, 0.0]), m=54.0, g_vec=array([0.0, 0.0, -9.81])
):
    w1 = m * (ddc - g_vec)
    return array(w1.tolist() + (X(c, w1) + dL).tolist())


def is_stable(
    H,
    c=array([0.0, 0.0, 0.0]),
    ddc=array([0.0, 0.0, 0.0]),
    dL=array([0.0, 0.0, 0.0]),
    m=54.0,
    g_vec=array([0.0, 0.0, -9.81]),
    robustness=0.0,
):
    w = compute_w(c, ddc, dL, m, g_vec)
    res = np.max(H.dot(w))
    if res > -robustness:
        print("offset ", res, dL)
        w = compute_w(c, ddc, array([0.0, 0.0, 0.0]), m, g_vec)
        print("offset2 ", np.max(H.dot(w)))
    return res <= -robustness


def allZeros(t):
    return zeros(3)


def __check_trajectory(p0, p1, p2, p3, T, H, mass, g, time_step=0.1, dL=allZeros):
    if time_step < 0:
        time_step = 0.01
    resolution = int(T / time_step)
    wps = [p0, p1, p2, p3]
    wps = matrix([pi.tolist() for pi in wps]).transpose()
    c_t = bezier(wps)
    ddc_t = c_t.compute_derivate(2)

    def c_tT(t):
        return asarray(c_t(t / T)).flatten()

    def ddc_tT(t):
        return 1.0 / (T * T) * asarray(ddc_t(t / T)).flatten()

    def dL_tT(t):
        return 1.0 / (T) * asarray(dL(t / T)).flatten()

    for i in range(resolution):
        t = T * float(i) / float(resolution)
        if not (
            is_stable(
                H,
                c=c_tT(t),
                ddc=ddc_tT(t),
                dL=dL_tT(t),
                m=mass,
                g_vec=g,
                robustness=-0.00001,
            )
        ):
            if t > 0.1:
                raise ValueError("trajectory is not stale ! at ", t)
            else:
                print(
                    is_stable(
                        H,
                        c=c_tT(t),
                        ddc=ddc_tT(t),
                        dL=asarray(dL(t / T)).flatten(),
                        m=mass,
                        g_vec=g,
                        robustness=-0.00001,
                    )
                )
                print("failed at 0")


###################
# LP BEZIER TESTS #
###################


def test_continuous_cpp_vs_continuous_py(N_CONTACTS=2, solver="qpoases", verb=0):
    mu = 0.5
    # friction coefficient
    lx = 0.1
    # half foot size in x direction
    ly = 0.07
    # half foot size in y direction
    # First, generate a contact configuration
    CONTACT_POINT_UPPER_BOUNDS = [0.5, 0.5, 0.5]
    CONTACT_POINT_LOWER_BOUNDS = [-0.5, -0.5, 0.0]
    gamma = atan(mu)
    # half friction cone angle
    RPY_LOWER_BOUNDS = [-2 * gamma, -2 * gamma, -pi]
    RPY_UPPER_BOUNDS = [+2 * gamma, +2 * gamma, +pi]
    MIN_CONTACT_DISTANCE = 0.3
    global mass
    global g_vector
    X_MARG = 0.07
    Y_MARG = 0.07

    succeeded = False
    while not succeeded:
        p, N = generate_contacts(
            N_CONTACTS,
            lx,
            ly,
            mu,
            CONTACT_POINT_LOWER_BOUNDS,
            CONTACT_POINT_UPPER_BOUNDS,
            RPY_LOWER_BOUNDS,
            RPY_UPPER_BOUNDS,
            MIN_CONTACT_DISTANCE,
            False,
        )
        X_LB = np.min(p[:, 0] - X_MARG)
        X_UB = np.max(p[:, 0] + X_MARG)
        Y_LB = np.min(p[:, 1] - Y_MARG)
        Y_UB = np.max(p[:, 1] + Y_MARG)
        Z_LB = np.max(p[:, 2] + 0.3)
        Z_UB = np.max(p[:, 2] + 1.5)
        # (H,h) = compute_GIWC(p, N, mu, False);
        H = -compute_CWC(p, N, mass, mu)
        h = zeros(H.shape[0])
        (succeeded, c0) = find_static_equilibrium_com(
            mass, [X_LB, Y_LB, Z_LB], [X_UB, Y_UB, Z_UB], H, h
        )

    rng = np.random.default_rng()
    dc0 = rng.uniform(-1, 1, size=3)

    Z_MIN = np.max(p[:, 2]) - 0.1
    Ineq_kin = zeros([3, 3])
    Ineq_kin[2, 2] = -1
    ineq_kin = zeros(3)
    ineq_kin[2] = -Z_MIN

    bezierSolver = BezierZeroStepCapturability(
        "ss",
        c0,
        dc0,
        p,
        N,
        mu,
        g_vector,
        mass,
        verb=verb,
        solver=solver,
        kinematic_constraints=[Ineq_kin, ineq_kin],
    )
    eqCpp = Equilibrium("dyn_eq2", mass, 4)
    eqCpp.setNewContacts(
        asmatrix(p), asmatrix(N), mu, EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_PP
    )
    # bezierSolver = BezierZeroStepCapturability("ss", c0, dc0, p, N, mu, g_vector,
    #                                            mass, verb=verb, solver=solver,
    # kinematic_constraints = None);
    # stabilitySolver = StabilityCriterion("ss", c0, dc0, p, N, mu, g_vector, mass,
    #                                       verb=verb, solver=solver);
    window_times = (
        [1] + [0.1 * i for i in range(1, 10)] + [0.1 * i for i in range(11, 21)]
    )  # try nominal time first
    # try nominal time first
    # window_times =  [0.2*i for i in range(1,5)] + [0.2*i for i in range(6,11)]
    # window_times = [1]+ [0.4*i for i in range(1,4)] #try nominal time first
    # window_times = [1]+ [0.4*i for i in range(3,6)] #try nominal time first
    # window_times = [0.7]
    found = False
    time_step_check = -0.2
    for i, el in enumerate(window_times):
        if found:
            break
        res = bezierSolver.can_I_stop(T=el, time_step=time_step_check)
        if res.is_stable:
            found = True
            # print("continuous found at ", el)
            __check_trajectory(
                bezierSolver._p0,
                bezierSolver._p1,
                res.c,
                res.c,
                el,
                bezierSolver._H,
                bezierSolver._mass,
                bezierSolver._g,
                time_step=time_step_check,
                dL=bezier(matrix([p_i.tolist() for p_i in res.wpsdL]).transpose()),
            )
            if i != 0:
                print("continuous Failed to stop at 1, but managed to stop at ", el)

    found = False
    time_step_check = 0.05
    for i, el in enumerate(window_times):
        if found:
            break
        res2 = bezierSolver.can_I_stop(T=el, time_step=time_step_check, l0=zeros(3))
        if res2.is_stable:
            found = True
            # print("ang_momentum found at ", el)
            __check_trajectory(
                bezierSolver._p0,
                bezierSolver._p1,
                res2.c,
                res2.c,
                el,
                bezierSolver._H,
                # bezierSolver._mass, bezierSolver._g, time_step =
                #               time_step_check, dL = res2.dL_of_t)
                bezierSolver._mass,
                bezierSolver._g,
                time_step=time_step_check,
                dL=bezier(matrix([p_i.tolist() for p_i in res2.wpsdL]).transpose()),
            )
            if i != 0:
                print("ang_momentum Failed to stop at 1, but managed to stop at ", el)
    # res2 = None
    # try:
    # res2 = stabilitySolver.can_I_stop();
    # except Exception as e:
    # pass

    if res2.is_stable != res.is_stable:
        if res.is_stable:
            print("continuous won")
        else:
            print("ang_momentum won")

    return res2.is_stable, res.is_stable, res2, res, c0, dc0, H, h, p, N


if __name__ == "__main__":
    g_vector = np.array([0.0, 0.0, -9.81])
    mass = 75
    # mass of the robot
    from matplotlib import rcParams

    rcParams.update({"font.size": 11})
    mine_won = 0
    mine_lose = 0
    total_stop = 0
    total_not_stop = 0
    total_disagree = 0
    margin_i_win_he_lose = []  # remaining speed
    margin_he_wins_i_lost = []  # remaining acceleration
    curves_when_i_win = []
    # times_disagree = []
    # times_agree_stop = []

    import matplotlib.pyplot as plt

    def __plot_3d_points(ax, points, c="b"):
        xs = [point[0] for point in points]
        ys = [point[1] for point in points]
        zs = [point[2] for point in points]
        ax.scatter(xs[:1], ys[:1], zs[:1], c="r")
        ax.scatter(xs[1:-1], ys[1:-1], zs[1:-1], c=c)
        ax.scatter(xs[-1:], ys[-1:], zs[-1:], c="g")
        ax.set_xlabel("X Label", fontsize=11)
        ax.set_ylabel("Y Label", fontsize=11)
        ax.set_zlabel("Z Label", fontsize=11)

    def plot_support_polygon(H, h, p, N, ax, c0):
        from pinocchio_inv_dyn.multi_contact.utils import compute_support_polygon

        # (H,h) = compute_GIWC(p, N, mu);
        global mass
        global g_vector
        (B_sp, b_sp) = compute_support_polygon(
            H, h, mass, g_vector, eliminate_redundancies=False
        )
        X_MIN = np.min(p[:, 0])
        X_MAX = np.max(p[:, 0])
        X_MIN -= 0.5 * (X_MAX - X_MIN)
        X_MAX += 0.5 * (X_MAX - X_MIN)
        Y_MIN = np.min(p[:, 1])
        Y_MAX = np.max(p[:, 1])
        Y_MIN -= 0.5 * (Y_MAX - Y_MIN)
        Y_MAX += 0.5 * (Y_MAX - Y_MIN)
        num_steps = 50
        dx = (X_MAX - X_MIN) / float(num_steps)
        dy = (Y_MAX - Y_MIN) / float(num_steps)
        # points = [
        # (X_MIN + dx * i, Y_MAX + dy * j, 0.0)
        # for i in range(num_steps + 1)
        # for j in range(num_steps + 1)
        # if check_static_eq(
        # H, h, mass, array([X_MIN + dx * i, Y_MAX + dy * j, 0.0]), g_vector
        # )
        # ]
        # points = [c0] + [
        # [X_MIN + dx * i, Y_MIN + dy * j, -0.5]
        # for i in range(num_steps + 1)
        # for j in range(num_steps + 1)
        # if check_static_eq(
        # H, h, mass, [X_MIN + dx * i, Y_MAX + dy * j, 0.0], g_vector
        # )
        # ]
        points = [c0] + [
            [X_MIN + dx * i, Y_MIN + dy * j, 0]
            for i in range(num_steps + 1)
            for j in range(num_steps + 1)
        ]
        pts2 = []
        for pt in points:
            if check_static_eq(H, h, mass, pt, g_vector):
                pts2 += [pt]
        __plot_3d_points(ax, pts2, c="r")
        # __plot_3d_points(ax, points2, c="r")
        __plot_3d_points(ax, p, c="r")
        # for i in range(num_steps):
        # for j in range(num_steps):
        # plot_inequalities(B_sp, b_sp, [X_MIN,X_MAX], [Y_MIN,Y_MAX],
        #                   ax=ax, color='b', lw=4, is_3d=False);
        # plot_inequalities(B_sp, b_sp, [X_MIN,X_MAX], [Y_MIN,Y_MAX],
        #                   ax=ax, color='b', lw=4, is_3d=False);
        # plt.show();

    def plot_win_curve(n=-1, num_pts=20):
        global curves_when_i_win
        if n > len(curves_when_i_win) - 1 or n < 0:
            print("n bigger than num curves or equal to -1, plotting last curve")
            n = len(curves_when_i_win) - 1
        (
            c0,
            dc0,
            c_end,
            dc_end,
            t_max,
            c_of_t,
            dc_of_t,
            ddc_of_t,
            H,
            h,
            p,
            N,
            dl_of_t,
            L_of_t,
        ) = curves_when_i_win[n]
        print("c0 ", c0)
        print("Is c0 stable ? ", check_static_eq(H, h, mass, c0, g_vector))
        print("Is end stable ? ", check_static_eq(H, h, mass, c_of_t(t_max), g_vector))

        w = np.zeros(6)
        w[2] = -mass * 9.81
        w[3:] = mass * np.cross(c_of_t(t_max), g_vector)
        print("max ", np.max(np.dot(H, w) - h))

        X_MIN = np.min(p[:, 0])
        X_MAX = np.max(p[:, 0])
        X_MIN -= 0.1 * (X_MAX - X_MIN)
        X_MAX += 0.1 * (X_MAX - X_MIN)
        Y_MIN = np.min(p[:, 1])
        Y_MAX = np.max(p[:, 1])
        print("Is XMIN ? ", X_MIN)
        print("Is XMAX ? ", X_MAX)
        print("Is YMIN ? ", Y_MIN)
        print("Is YMAX ? ", Y_MAX)
        delta = t_max / float(num_pts)
        num_pts += 1
        fig = plt.figure()
        ax = fig.add_subplot(221, projection="3d")
        # ax = fig.add_subplot(221)
        __plot_3d_points(ax, [c_of_t(i * delta) for i in range(num_pts)])
        __plot_3d_points(
            ax, [c0 + (c_end - c0) * i * delta for i in range(num_pts)], c="y"
        )
        plot_support_polygon(H, h, p, N, ax, c0)
        ax = fig.add_subplot(222, projection="3d")
        __plot_3d_points(ax, [dc_of_t(i * delta) for i in range(num_pts)])
        __plot_3d_points(
            ax, [dc0 + (dc_end - dc0) * i * delta for i in range(num_pts)], c="y"
        )
        ax = fig.add_subplot(223, projection="3d")
        # __plot_3d_points(ax, [ddc_of_t(i * delta) for i in range(num_pts)])
        # ax = fig.add_subplot(224, projection='3d')
        __plot_3d_points(ax, [L_of_t(i * delta) for i in range(num_pts)], c="y")
        # __plot_3d_points(ax, [-dc0* i * delta for i in range(num_pts)], c = "y")
        ax = fig.add_subplot(224, projection="3d")
        __plot_3d_points(ax, [dl_of_t(i * delta) for i in range(num_pts)])
        # ax = fig.add_subplot(121, projection='3d')
        # __plot_3d_points(ax, [ddc_of_t(i * delta) for i in range(num_pts)])
        # ax = fig.add_subplot(122, projection='3d')
        # __plot_3d_points(ax, [-dc0* i * delta for i in range(num_pts)])
        # print("cross product ", X(-dc0,ddc_of_t(0.5) - ddc_of_t(0) )
        #       / norm(X(-dc0,ddc_of_t(0.5) - ddc_of_t(0) )))
        # print("init acceleration ", ddc_of_t(0))
        print("init velocity ", dc_of_t(0))
        print("end velocity ", dc_of_t(t_max))
        # print("cross product ", X(-dc0,ddc_of_t(t_max) - ddc_of_t(0) )
        #       / norm(X(-dc0,ddc_of_t(t_max) - ddc_of_t(0))))

        # plt.show()

    def plot_n_win_curves(n=-1, num_pts=50):
        global curves_when_i_win
        if n > len(curves_when_i_win) - 1 or n < 0:
            print("n bigger than num curves or equal to -1, plotting last curve")
            n = len(curves_when_i_win) - 1
        for i in range(n):
            plot_win_curve(i, num_pts)
        plt.show()

    num_tested = 0.0
    for i in range(1000):
        num_tested = i - 1
        (
            mine,
            theirs,
            r_mine,
            r_theirs,
            c0,
            dc0,
            H,
            h,
            p,
            N,
        ) = test_continuous_cpp_vs_continuous_py()
        # mine, theirs, r_mine, r_theirs, c0, dc0, H,h, p, N = \
        # test_momentum_cpp_vs_momentum_py()
        # print("H test", H.shape)
        if mine != theirs:
            total_disagree += 1
            if mine:
                # times_disagree +=[r_mine.t]
                # raise ValueError (" BITE ME " )
                margin_i_win_he_lose += [r_theirs.dc]
                curves_when_i_win += [
                    (
                        c0[:],
                        dc0[:],
                        r_theirs.c[:],
                        r_theirs.dc[:],
                        r_mine.t,
                        r_mine.c_of_t,
                        r_mine.dc_of_t,
                        r_mine.ddc_of_t,
                        H[:],
                        h[:],
                        p[:],
                        N,
                        r_mine.dL_of_t,
                        r_mine.L_of_t,
                    )
                ]
                # print("margin when he lost: ", norm(r_theirs.dc))
            # else:
            # times_disagree +=[r_theirs.t]
            if mine:
                mine_won += 1
            else:
                mine_lose += 1
        elif mine or theirs:
            total_stop += 1
            # times_agree_stop+=[r_mine.t]
            margin_he_wins_i_lost += [r_theirs.ddc_min]

            # margin_i_win_he_lose+=[r_theirs.dc]
            # curves_when_i_win+=[(c0[:], dc0[:], r_theirs.c[:],
            # r_theirs.dc[:], r_mine.t, r_mine.c_of_t,
            # r_mine.dc_of_t, r_mine.ddc_of_t, H[:], h[:], p[:], N)]
            # print("margin when he wins: ", r_theirs.ddc_min)
        else:
            total_not_stop += 1

    print("% of stops", 100.0 * float(total_stop) / num_tested, total_stop)
    print(
        "% of total_disagree",
        100.0 * float(total_disagree) / num_tested,
        total_disagree,
    )
    if total_disagree > 0:
        print("% of wins", 100.0 * float(mine_won) / total_disagree)
        print("% of lose", 100.0 * float(mine_lose) / total_disagree)
