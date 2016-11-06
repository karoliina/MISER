import numpy as np
import matplotlib.pyplot as plt
from RandomNumberGenerators import Ran
import pdb

def circle(point):
    x = point[0]
    y = point[1]
    if x**2 + y**2 < 1:
        return 1
    else:
        return 0

def paraboloid(point):
    x = point[0]
    y = point[1]
    return x**2 + y**2

def y1(x):
    return 3*np.sqrt(1 - x**2/49)

def y2(x):
    return (6*np.sqrt(10)/7 + (1.5 - 0.5*abs(x))) - \
           (6*np.sqrt(10)/14)*np.sqrt(4 - (abs(x) - 1)**2)

def y3(x):
    return 9 - 8*abs(x)

def y4(x):
    return 3*abs(x) + 0.75

def y5(x):
    return (abs(x/2) - (3*np.sqrt(33) - 7)/112*x**2 - 3) + \
            np.sqrt(1 - (abs(abs(x) - 2) - 1)**2)

def batman(point):
    """Checks if the point (x,y) is within the Batman curve.

    :param point: array-like - contains the x and y coordinates of the point
    :return: 1 if point is inside the Batman, 0 if it isn't
    """
    x = point[0]
    y = point[1]
    absx = abs(x)
    absy = abs(y)
    fxl = fxu = 0
    if absx > 7 or absy > 3:
        return 0
    elif absx > 4 and absx <= 7:
        fxu = y1(x)
        fxl = -1*fxu
    elif absx > 3 and absx <= 4:
        fxu = y1(x)
        fxl = y5(x)
    elif absx > 1 and absx <= 3:
        fxu = y2(x)
        fxl = y5(x)
    elif absx > 0.75 and absx <= 1:
        fxu = y3(x)
        fxl = y5(x)
    elif absx > 0.5 and absx <= 0.75:
        fxu = y4(x)
        fxl = y5(x)
    elif absx <= 0.5:
        fxu = 2.25
        fxl = y5(x)
    if y >= fxl and y <= fxu:
        return 1
    else:
        return 0


def hyperbatman(point):
    """Checks if the point (x0, x1, ..., xn) is within a n-dimensional Batman-shaped cylinder, where
    the cylinder extends from -1 to 1 in the dimensions 2..n.

    :param point: array-like - contains the x and y coordinates of the point
    :return: 1 if point is inside the Batman, 0 if it isn't
    """
    x = point[0]
    y = point[1]
    absx = abs(x)
    absy = abs(y)
    fxl = fxu = 0
    if absx > 7 or absy > 3:
        return 0
    elif absx > 4 and absx <= 7:
        fxu = y1(x)
        fxl = -1*fxu
    elif absx > 3 and absx <= 4:
        fxu = y1(x)
        fxl = y5(x)
    elif absx > 1 and absx <= 3:
        fxu = y2(x)
        fxl = y5(x)
    elif absx > 0.75 and absx <= 1:
        fxu = y3(x)
        fxl = y5(x)
    elif absx > 0.5 and absx <= 0.75:
        fxu = y4(x)
        fxl = y5(x)
    elif absx <= 0.5:
        fxu = 2.25
        fxl = y5(x)
    if y >= fxl and y <= fxu:
        for i in range(2, len(point)):
            if point[i] < -1 or point[i] > 1:
                return 0
        return 1
    else:
        return 0


batman_area = 955/48 - 2/7*(2*np.sqrt(33) + 7*np.pi + 3*np.sqrt(10)*(np.pi - 1)) \
              + 21*(np.arccos(3/7) + np.arccos(4/7))

# x = np.linspace(-7, 7, 500)
# yu = np.zeros(500)
# yl = np.zeros(500)
# for i in range(500):
#     yu[i] = batman_upper(x[i])
#     yl[i] = batman_lower(x[i])

def batman_upper(x):
    if x >= -7 and x < -3:
        return y1(x)
    elif x >= -3 and x < -1:
        return y2(x)
    elif x >= -1 and x < -0.75:
        return y3(x)
    elif x >= -0.75 and x < -0.5:
        return y4(x)
    elif x >= -0.5 and x <= 0.5:
        return 2.25
    elif x > 0.5 and x <= 0.75:
        return y4(x)
    elif x > 0.75 and x <= 1:
        return y3(x)
    elif x > 1 and x <= 3:
        return y2(x)
    elif x > 3 and x <= 7:
        return y1(x)
    else:
        return float('NaN')

def batman_lower(x):
    if x >= -7 and x < -4:
        return -1*y1(x)
    elif x >= -4 and x <= 4:
        return y5(x)
    elif x > 4 and x <= 7:
        return -1*y1(x)
    else:
        return float('NaN')

def batman_dy(x):
    return batman_upper(x) - batman_lower(x)


# ran = Ran(5331)
#
# def ranpt(low, high):
#     """Returns a uniformly random point in an n-dimensional rectangular region.
#
#     :param low: array-like - the "lower left" coordinates
#     :param high: array-like - the "upper right" coordinates
#     :return: np.array - the n-dimensional random point
#     """
#     n = len(low)
#     pt = np.zeros(n)
#     for i in range(n):
#         pt[i] = ran.rand(low=low[i], high=high[i])
#     return pt
#
# fig, ax = plt.subplots()
# plt.xlim([-7, 7])
# plt.ylim([-3, 3])
# x = np.linspace(-7, 7, 500)
# batman_up = [batman_upper(pt) for pt in x]
# batman_low = [batman_lower(pt) for pt in x]
# plt.plot(x, batman_up, 'r-')
# plt.plot(x, batman_low, 'r-')
#
# for i in range(100):
#     pt = ranpt([-7, -3], [7, 3])
#     print('[{:.4g}, {:.4g}]: f = {}'.format(pt[0], pt[1], batman(pt)))
#     ax.plot(pt[0], pt[1], 'b.')
#     fig.canvas.draw()
#     plt.show()
#     pdb.set_trace()
