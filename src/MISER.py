from RandomNumberGenerators import Ran
from test_functions import batman, circle, batman_upper, batman_lower, hyperbatman, batman_area
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import pdb

class MISER:

    def __init__(self):
        """Initialize the random number generators used with the MISER function."""
        self.ran = Ran(5331)
        self.points = []
        self.terminals = []

    def ranpt(self, low, high):
        """Returns a uniformly random point in an n-dimensional rectangular region.

        :param low: array-like - the "lower left" coordinates
        :param high: array-like - the "upper right" coordinates
        :return: np.array - the n-dimensional random point
        """
        n = len(low)
        pt = np.zeros(n)
        for i in range(n):
            pt[i] = self.ran.rand(low=low[i], high=high[i])
        return pt

    def mean_variance(self, a):
        """Computes and returns estimates of the mean and variance of the array a.

        :param a: array-like - an 1D array of numbers
        :return: mean - the estimated mean of a
                 var - the estimated variance of a
        """
        n = len(a)
        if n == 0:
            return 0, 0
        a2 = list(map(lambda x: x**2, a))
        mean = sum(a) / n
        return mean, sum(a2) / n - mean**2

    def MISER(self, func, low, high, N, dith):
        """Implementation of the MISER algorithm of Press and Farrar. Monte Carlo samples the function func in the
        n-dimensional rectangular volume region using recursive stratified sampling.

        :param func: function from R^n to R
        :param low: array-like - the n "lower-left" coordinates
        :param high: array-like - the n "upper-right" coordinates
        :param N: integer - total number of points to sample
        :param dith: float - the dithering parameter, a nonzero value will cause the subvolumes to be divided at a
                     point 0.5 +- dith instead of in the middle
        :return: ave - estimate of the mean of func within the region
                 err - error estimate of the integral
        """
        MNBS = 60 # minimum number of evaluations left so that a subregion will be further bisected
        MNPT = 15 # minimum number of evaluations in a terminal subregion
        PFAC = 0.1 # fraction of evaluations to use for exploring the variance
        BIG = 1e30
        n = len(low)
        V = high[0] - low[0]
        for i in range(1, n):
            V *= high[i] - low[i]
        N_pre = max(MNPT, int(N*PFAC))
        # too few overall points or preliminary sampling points left,
        # do straight Monte Carlo
        if N < MNBS:
            self.terminals.append(np.concatenate((low, high)))
            fvals = np.zeros(N)
            for i in range(N):
                pt = self.ranpt(low, high)
                self.points.append(pt)
                fvals[i] = func(pt)
            avg, var = self.mean_variance(fvals)
            # print('avg = {:.4g}'.format(avg))
            return avg, var/N, N
        else: # perform preliminary sampling
            mid = np.zeros(n)
            # initialize the midpoints for each dimension
            for j in range(n):
                # if dith != 0, select midpoint s randomly
                s = np.sign(self.ran.rand(-1, 1)) * dith
                mid[j] = (0.5 + s)*low[j] + (0.5 - s)*high[j]
            left_fvals = []
            right_fvals = []
            for i in range(n):
                left_fvals.append([])
                right_fvals.append([])
            for i in range(N_pre):
                pt = self.ranpt(low, high)
                fval = func(pt)
                # for each dimension, accumulate sums for the subregion pt fell into
                for j in range(n):
                    if pt[j] <= mid[j]: # point is in left subregion
                        left_fvals[j].append(fval)
                    else: # point is in right subregion
                        right_fvals[j].append(fval)
            # choose which dimension jb to bisect
            sumb = BIG
            jb = -1
            ssums = np.zeros(n) # variance estimators for each dimension
            left_variances = np.zeros(n)
            right_variances = np.zeros(n)
            for j in range(n):
                mean_left, left_variances[j] = self.mean_variance(left_fvals[j])
                mean_right, right_variances[j] = self.mean_variance(right_fvals[j])
                sigma_left = left_variances[j]**(1/3)
                sigma_right = right_variances[j]**(1/3)
                ssums[j] = sigma_left + sigma_right
                if ssums[j] <= sumb:
                    sumb = ssums[j]
                    jb = j
                    sigma_left_b = sigma_left
                    sigma_right_b = sigma_right
            # if the samples N_pre were in the same half for all dimensions,
            # choose jb randomly (this can happen if MNPT is too small).
            # also choose jb randomly if the variances are equal in all dimensions
            if jb == -1 or len(set(ssums)) == 1:
                jb = int(self.ran.rand(0, n))
            reg_left = low[jb]
            reg_mid = mid[jb]
            reg_right = high[jb]
            frac_left = abs((reg_mid - reg_left)/(reg_right - reg_left))
            frac_right = 1 - frac_left
            # number of points allocated to left subregion - equation (20)
            if sigma_left_b > 0 or sigma_right_b > 0:
                N_left = MNPT + int((N - N_pre - 2*MNPT)*(frac_left*sigma_left_b / \
                    (frac_left*sigma_left_b + frac_right*sigma_right_b)))
            else:
                N_left = MNPT + int((N - N_pre - 2*MNPT)/2)
            # number of points allocated to right subregion
            N_right = N - N_pre - N_left
            # gather limits of subregions
            low_tmp = np.zeros(n)
            high_tmp = np.zeros(n)
            for j in range(n):
                low_tmp[j] = low[j]
                high_tmp[j] = high[j]
            high_tmp[jb] = mid[jb]
            # recursion in left subregion
            ave_left, var_left, Nl = self.MISER(func, low_tmp, high_tmp, N_left, dith)
            low_tmp[jb] = mid[jb]
            high_tmp[jb] = high[jb]
            # recursion in right subregion
            ave_right, var_right, Nr = self.MISER(func, low_tmp, high_tmp, N_right, dith)
            # combine regions
            ave = frac_left*ave_left + frac_right*ave_right
            var = var_left*frac_left**2 + var_right*frac_right**2
            return ave, var, Nl + Nr

    def get_points(self):
        """Return a list of points for each dimension of the most recent
        calculation."""
        N = len(self.points)
        if N == 0:
            return None
        n = len(self.points[0])
        pts = []
        for j in range(n):
            tmp = np.zeros(N)
            for i in range(N):
                tmp[i] = self.points[i][j]
            pts.append(tmp.copy())
        return pts


