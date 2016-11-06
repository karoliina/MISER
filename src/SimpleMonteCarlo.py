from RandomNumberGenerators import Ran
from test_functions import hyperbatman, batman_area, batman, batman_upper, batman_lower
import numpy as np
import matplotlib.pyplot as plt

def estimate_pi(R, N):
    """Estimates the value of pi by generating random points within a square with
    the lower left corner at (-R, R) and the upper right corner at (R, R), and counting
    the number of them that fall inside a circle centered at (0, 0) with radius R.

    :param R: float - the radius of the circle
    :param N: the number of points
    :return: float - the estimated value of pi
    """
    pts = np.zeros((N, 2))
    ran = Ran(5331)
    Nc = 0 # points inside a circle of radius R
    for i in range(N):
        x = ran.rand(-R, R)
        y = ran.rand(-R, R)
        pts[i,0] = x
        pts[i,1] = y
        if x**2 + y**2 <= R**2:
            Nc += 1
    return 4*Nc/N, pts

def mc_integrate(func, low, high, N):
    """Integrates the function func over the given rectangular region
    using a simple Monte Carlo method with N sample points.

    :param func: function - the function to integrate, R^n -> R
    :param low: array-like - the lower-left coordinates of the region
    :param high: array-like - the upper-right coordinates of the region
    :param N: integer - the number of n-dimensional sample points to use
    :return: integral - float, the value of the integral
             error - float, the one standard deviation error estimate
    """
    ran = Ran(5331)
    n = len(low)
    avg = 0
    avg2 = 0
    pts = []
    # compute the volume of the region
    V = high[0] - low[0]
    for j in range(1, n):
        V *= high[j] - low[j]
    for i in range(N):
        pt = np.zeros(n)
        for j in range(n):
            pt[j] = ran.rand(low[j], high[j])
        fval = func(pt)
        avg += fval
        avg2 += fval**2
        pts.append(pt)
    avg /= N
    avg2 /= N
    error = V*np.sqrt((avg2 - avg**2)/N)
    return avg*V, error, pts

def get_points(pts):
    N = len(pts)
    if N == 0:
        return None
    n = len(pts[0])
    out_pts = []
    for j in range(n):
        tmp = np.zeros(N)
        for i in range(N):
            tmp[i] = pts[i][j]
        out_pts.append(tmp.copy())
    return out_pts

"""
# batman
low = [-7, -3]
high = [7, 3]
ave, var, pts = mc_integrate(batman, low, high, int(5e4))
pts = get_points(pts)
fig, ax = plt.subplots(figsize=(9,5))
plt.xlim([low[0], high[0]])
plt.ylim([low[1], high[1]])
patches = []

# for t in miser.terminals:
#     x = t[0]
#     y = t[1]
#     width = t[2] - t[0]
#     height = t[3] - t[1]
#     patches.append(Rectangle((x,y), width, height, facecolor=(0,0,0,0), edgecolor='#eeeeee'))

# batman curve
x = np.linspace(low[0], high[0], 1000)
batman_up = [batman_upper(pt) for pt in x]
batman_low = [batman_lower(pt) for pt in x]

ax.fill_between(x, batman_low, batman_up, facecolor='#fde311')
ax.fill_between(x, batman_up, 7*np.ones(1000), facecolor='#000000')
ax.fill_between(x, -7*np.ones(1000), batman_low, facecolor='#000000')

# points
plt.xlabel('x')
plt.ylabel('y')
plt.plot(x, batman_up, 'k-', linewidth=2)
plt.plot(x, batman_low, 'k-', linewidth=2)
plt.plot(pts[0], pts[1], 'r,')

plt.show()

# TESTING - multidimensional simple Monte Carlo
# the test function
# func = lambda x: sum(x)**2
# # the region
# low = np.zeros(10)
# high = np.ones(10)
# # values of N to test with
# N = [10**i for i in range(3, 8)]
# Nlen = len(N)
# sintegral = np.zeros(Nlen)
# serror = np.zeros(Nlen)
# for i in range(Nlen):
#     print('N = {:.4g}'.format(N[i]))
#     sintegral[i], serror[i], pts = mc_integrate(func, low, high, N[i])

# R = 1
# pi, pts = estimate_pi(R, 10**6)
# fig, ax = plt.subplots(figsize=(5,5))
# plt.xlim([-R-0.00, R+0.00])
# plt.ylim([-R-0.00, R+0.00])
# x1 = np.linspace(-R, R, 1000)
# y1 = np.sqrt(R**2 - x1**2)
# y2 = -1*y1
# plt.plot(x1, y1, "k-", lw=2)
# plt.plot(x1, y2, "k-", lw=2)
# inside_pts = list(filter(lambda x: x[0]**2 + x[1]**2 < R**2, pts))
# inx = []
# iny = []
# for pt in inside_pts:
#     inx.append(pt[0])
#     iny.append(pt[1])
# plt.plot(pts[:,0], pts[:,1], 'b,')
# plt.plot(inx, iny, 'r,')
# plt.xlabel('x')
# plt.ylabel('y')
#
# plt.show()

# save, svar = mc_integrate(hyperbatman, [-7, -3, -1], [7, 3, 1], 10**5)

# ns = [i for i in range(3, 11)]
# sactual = []
# scomputed = []
# serror = []
# srel = []
# j = 0
# for n in ns:
#     print('n =', n)
#     low = [-7, -3]
#     high = [7, 3]
#     for i in range(2, n):
#         low.append(-1)
#         high.append(1)
#     ave, var = mc_integrate(hyperbatman, low, high, 10**6)
#     scomputed.append(ave)
#     sactual.append(batman_area*2**(n-2))
#     srel.append(abs(sactual[j]-scomputed[j])/abs(sactual[j]))
#     serror.append(var)
#     j += 1

# Ns = [10**i for i in range(3, 8)]
# V = 14*7*2**6
# actual = batman_area*2**6
# low = [-7, -3, -1, -1, -1, -1, -1, -1]
# high = [7, 3, 1, 1, 1, 1, 1, 1]
# scomputed = []
# serror = []
# srel = []
# for N in Ns:
#     print('N =', N)
#     ave, var = mc_integrate(hyperbatman, low, high, N)
#     scomputed.append(ave)
#     serror.append(var)
#     srel.append(abs(actual-ave)/abs(actual))
"""
