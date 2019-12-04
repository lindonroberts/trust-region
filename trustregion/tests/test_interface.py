# Ensure compatibility with Python 2
from __future__ import absolute_import, division, print_function, unicode_literals

from math import sqrt
import numpy as np
import unittest

import trustregion


BALL_EPS=1e-12


def model_value(g, H, s):
    return np.dot(g, s) + 0.5*np.dot(s, H.dot(s))


def cauchy_pt(g, hess, delta):
    # General expression for the Cauchy point
    crv = np.dot(g, hess.dot(g))
    gnorm = np.linalg.norm(g)
    if crv <= 0.0:
        alpha = delta / gnorm
    else:
        alpha = min(delta / gnorm, gnorm**2 / crv)
    s = -alpha * g
    red = model_value(g, hess, s)
    crvmin = np.dot(s, hess.dot(s)) / np.dot(s, s)
    if crvmin < 0.0:
        crvmin = -1.0
    return s, red, crvmin


def cauchy_pt_box(g, hess, delta, lower, upper):
    # General expression for the Cauchy point, lower <= s <= upper
    crv = np.dot(g, hess.dot(g))
    gnorm = np.linalg.norm(g)
    if crv <= 0.0:
        alpha = delta / gnorm
    else:
        alpha = min(delta / gnorm, gnorm**2 / crv)
    # print("alpha = %g" % alpha)
    # Then cap with bounds:
    for i in range(len(g)):
        if g[i] > 0:  # s[i] negative, will hit lower
            alpha = min(alpha, -lower[i] / g[i])
        elif g[i] < 0:  # s[i] positive, will hit upper
            alpha = min(alpha, -upper[i] / g[i])
        # print("alpha = %g after i=%g" % (alpha, i))
    s = -alpha * g
    red = model_value(g, hess, s)
    crvmin = np.dot(s, hess.dot(s)) / np.dot(s, s)
    if crvmin < 0.0:
        crvmin = -1.0
    return s, red, crvmin


class TestUncInternal(unittest.TestCase):
    def runTest(self):
        n = 3
        g = np.array([1.0, 0.0, 1.0])
        H = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
        Delta = 2.0
        d, gnew, crvmin = trustregion.solve(g, H, Delta, verbose_output=True)
        true_d = np.array([-1.0, 0.0, -0.5])
        est_min = model_value(g, H, d)
        true_min = model_value(g, H, true_d)
        # Hope to get actual correct answer for internal minimum?
        # self.assertTrue(np.all(d == true_d), msg='Wrong answer')
        # self.assertAlmostEqual(est_min, true_min, msg='Wrong min value')
        s_cauchy, red_cauchy, crvmin_cauchy = cauchy_pt(g, H, Delta)
        self.assertTrue(est_min <= red_cauchy, msg='Cauchy reduction not achieved')
        self.assertTrue(np.all(gnew == g + H.dot(d)), msg='Wrong gnew')
        # print(crvmin)
        self.assertAlmostEqual(crvmin, 1.2, msg='Wrong crvmin')
        self.assertTrue(np.linalg.norm(d) <= Delta + BALL_EPS, msg='Ball constraint violated')


class TestUncBdry(unittest.TestCase):
    def runTest(self):
        n = 3
        g = np.array([1.0, 0.0, 1.0])
        H = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
        Delta = 5.0 / 12.0
        d, gnew, crvmin = trustregion.solve(g, H, Delta, verbose_output=True)
        true_d = np.array([-1.0 / 3.0, 0.0, -0.25])
        est_min = model_value(g, H, d)
        true_min = model_value(g, H, true_d)
        # Hope to get actual correct answer
        # self.assertTrue(np.all(d == true_d), msg='Wrong answer')
        # self.assertAlmostEqual(est_min, true_min, msg='Wrong min value')
        s_cauchy, red_cauchy, crvmin_cauchy = cauchy_pt(g, H, Delta)
        self.assertTrue(est_min <= red_cauchy, msg='Cauchy reduction not achieved')
        self.assertTrue(np.all(gnew == g + H.dot(d)), msg='Wrong gnew')
        self.assertAlmostEqual(crvmin, 0.0, msg='Wrong crvmin')
        self.assertTrue(np.linalg.norm(d) <= Delta+BALL_EPS, msg='Ball constraint violated')


class TestUncBdry2(unittest.TestCase):
    def runTest(self):
        n = 3
        g = np.array([1.0, 0.0, 1.0])
        H = np.array([[-2.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
        Delta = 5.0 / 12.0
        d, gnew, crvmin = trustregion.solve(g, H, Delta, verbose_output=True)
        true_d = np.array([-1.0 / 3.0, 0.0, -0.25])
        est_min = model_value(g, H, d)
        true_min = model_value(g, H, true_d)
        # Hope to get actual correct answer
        # self.assertTrue(np.all(d == true_d), msg='Wrong answer')
        # self.assertAlmostEqual(est_min, true_min, msg='Wrong min value')
        s_cauchy, red_cauchy, crvmin_cauchy = cauchy_pt(g, H, Delta)
        self.assertTrue(est_min <= red_cauchy, msg='Cauchy reduction not achieved')
        self.assertTrue(np.all(gnew == g + H.dot(d)), msg='Wrong gnew')
        self.assertAlmostEqual(crvmin, 0.0, msg='Wrong crvmin')
        self.assertTrue(np.linalg.norm(d) <= Delta+BALL_EPS, msg='Ball constraint violated')


class TestUncBdry3(unittest.TestCase):
    def runTest(self):
        n = 3
        g = np.array([0.0, 0.0, 1.0])
        H = np.array([[-2.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
        Delta = 0.5
        d, gnew, crvmin = trustregion.solve(g, H, Delta, verbose_output=True)
        true_d = np.array([0.0, 0.0, -0.5])
        est_min = model_value(g, H, d)
        true_min = model_value(g, H, true_d)
        # Hope to get actual correct answer
        # self.assertTrue(np.all(d == true_d), msg='Wrong answer')
        # self.assertAlmostEqual(est_min, true_min, msg='Wrong min value')
        s_cauchy, red_cauchy, crvmin_cauchy = cauchy_pt(g, H, Delta)
        self.assertTrue(est_min <= red_cauchy, msg='Cauchy reduction not achieved')
        self.assertTrue(np.all(gnew == g + H.dot(d)), msg='Wrong gnew')
        self.assertAlmostEqual(crvmin, 0.0, msg='Wrong crvmin')
        # self.assertAlmostEqual(crvmin, crvmin_cauchy, msg='Wrong crvmin')
        self.assertTrue(np.linalg.norm(d) <= Delta+BALL_EPS, msg='Ball constraint violated')


class TestUncHard(unittest.TestCase):
    def runTest(self):
        n = 3
        g = np.array([0.0, 0.0, 1.0])
        H = np.array([[-2.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]])
        Delta = sqrt(2.0)
        d, gnew, crvmin = trustregion.solve(g, H, Delta, verbose_output=True)
        true_d = np.array([1.0, 0.0, -1.0])  # non-unique solution
        est_min = model_value(g, H, d)
        true_min = model_value(g, H, true_d)
        # Hope to get actual correct answer
        # self.assertTrue(np.all(d == true_d), msg='Wrong answer')
        # self.assertAlmostEqual(est_min, true_min, msg='Wrong min value')
        s_cauchy, red_cauchy, crvmin_cauchy = cauchy_pt(g, H, Delta)
        self.assertTrue(est_min <= red_cauchy, msg='Cauchy reduction not achieved')
        self.assertTrue(np.all(gnew == g + H.dot(d)), msg='Wrong gnew')
        self.assertAlmostEqual(crvmin, 0.0, msg='Wrong crvmin')
        self.assertTrue(np.linalg.norm(d) <= Delta+BALL_EPS, msg='Ball constraint violated')


class TestConInternal(unittest.TestCase):
    def runTest(self):
        n = 3
        g = np.array([1.0, 0.0, 1.0])
        H = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
        Delta = 2.0
        sl = np.array([-0.5, -10.0, -10.0])
        su = np.array([10.0, 10.0, 10.0])
        d, gnew, crvmin = trustregion.solve(g, H, Delta, sl=sl, su=su, verbose_output=True)
        true_d = np.array([-1.0, 0.0, -0.5])
        est_min = model_value(g, H, d)
        true_min = model_value(g, H, true_d)
        # Hope to get actual correct answer for internal minimum?
        # self.assertTrue(np.all(d == true_d), msg='Wrong answer')
        # self.assertAlmostEqual(est_min, true_min, msg='Wrong min value')
        s_cauchy, red_cauchy, crvmin_cauchy = cauchy_pt_box(g, H, Delta, sl, su)
        # print(s_cauchy)
        # print(d)
        self.assertTrue(est_min <= red_cauchy, msg='Cauchy reduction not achieved')
        self.assertTrue(np.all(gnew == g + H.dot(d)), msg='Wrong gnew')
        # print(crvmin)
        self.assertAlmostEqual(crvmin, -1.0, msg='Wrong crvmin')
        self.assertTrue(np.linalg.norm(d) <= Delta+BALL_EPS, msg='Ball constraint violated')
        self.assertTrue(np.max(d - sl) >= 0.0, msg='Lower bound violated')
        self.assertTrue(np.max(d - su) <= 0.0, msg='Upper bound violated')


class TestConBdry(unittest.TestCase):
    def runTest(self):
        n = 3
        g = np.array([1.0, 0.0, 1.0])
        H = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
        Delta = 5.0 / 12.0
        sl = np.array([-0.3, -0.01, -0.1])
        su = np.array([10.0, 1.0, 10.0])
        d, gnew, crvmin = trustregion.solve(g, H, Delta, sl=sl, su=su, verbose_output=True)
        true_d = np.array([-1.0 / 3.0, 0.0, -0.25])
        est_min = model_value(g, H, d)
        true_min = model_value(g, H, true_d)
        # Hope to get actual correct answer
        # self.assertTrue(np.all(d == true_d), msg='Wrong answer')
        # self.assertAlmostEqual(est_min, true_min, msg='Wrong min value')
        s_cauchy, red_cauchy, crvmin_cauchy = cauchy_pt_box(g, H, Delta, sl, su)
        self.assertTrue(est_min <= red_cauchy, msg='Cauchy reduction not achieved')
        self.assertTrue(np.max(np.abs(gnew - g - H.dot(d))) < 1e-10, msg='Wrong gnew')
        # print(crvmin)
        self.assertAlmostEqual(crvmin, -1.0, msg='Wrong crvmin')
        # self.assertAlmostEqual(crvmin, crvmin_cauchy, msg='Wrong crvmin')
        self.assertTrue(np.linalg.norm(d) <= Delta+BALL_EPS, msg='Ball constraint violated')
        self.assertTrue(np.max(d - sl) >= 0.0, msg='Lower bound violated')
        self.assertTrue(np.max(d - su) <= 0.0, msg='Upper bound violated')


class TestLinear(unittest.TestCase):
    def runTest(self):
        g = np.array([1.0, -1.0])
        a = np.array([-2.0, -2.0])
        b = np.array([1.0, 2.0])
        delta = 2.0
        for H in [None, np.zeros((len(g), len(g)))]:
            x = trustregion.solve(g, H, delta, sl=a, su=b)
            xtrue = np.array([-sqrt(2.0), sqrt(2.0)])
            self.assertTrue(np.max(np.abs(x - xtrue)) < 1e-10, msg='Wrong step')
            self.assertTrue(np.linalg.norm(x) <= delta+BALL_EPS, msg='Ball constraint violated')
            self.assertTrue(np.max(x - a) >= 0.0, msg='Lower bound violated')
            self.assertTrue(np.max(x - b) <= 0.0, msg='Upper bound violated')


class TestLinear2(unittest.TestCase):
    def runTest(self):
        g = np.array([1.0, -1.0])
        a = np.array([-2.0, -2.0])
        b = np.array([1.0, 2.0])
        delta = 5.0
        for H in [None, np.zeros((len(g), len(g)))]:
            x = trustregion.solve(g, H, delta, sl=a, su=b)
            xtrue = np.array([-2.0, 2.0])
            self.assertTrue(np.max(np.abs(x - xtrue)) < 1e-10, msg='Wrong step')
            self.assertTrue(np.linalg.norm(x) <= delta+BALL_EPS, msg='Ball constraint violated')
            self.assertTrue(np.max(x - a) >= 0.0, msg='Lower bound violated')
            self.assertTrue(np.max(x - b) <= 0.0, msg='Upper bound violated')


class TestLinearOldBug(unittest.TestCase):
    def runTest(self):
        g = np.array([-1.0, -1.0])
        a = np.array([-2.0, -2.0])
        b = np.array([0.1, 0.9])
        delta = sqrt(2.0)
        for H in [None, np.zeros((len(g), len(g)))]:
            x = trustregion.solve(g, H, delta, sl=a, su=b)
            xtrue = b
            # print(x)
            self.assertTrue(np.max(np.abs(x - xtrue)) < 1e-10, msg='Wrong step')
            self.assertTrue(np.linalg.norm(x) <= delta+BALL_EPS, msg='Ball constraint violated')
            self.assertTrue(np.max(x - a) >= 0.0, msg='Lower bound violated')
            self.assertTrue(np.max(x - b) <= 0.0, msg='Upper bound violated')


class TestLinearOldBug2(unittest.TestCase):
    def runTest(self):
        g = np.array([-1.0, -1.0, -1.0])
        a = np.array([-2.0, -2.0, -2.0])
        b = np.array([0.9, 0.1, 5.0])
        delta = sqrt(3.0)
        for H in [None, np.zeros((len(g), len(g)))]:
            x = trustregion.solve(g, H, delta, sl=a, su=b)
            xtrue = np.array([0.9, 0.1, sqrt(3.0 - 0.81 - 0.01)])
            # print(x)
            self.assertTrue(np.max(np.abs(x - xtrue)) < 1e-10, msg='Wrong step')
            self.assertTrue(np.linalg.norm(x) <= delta+BALL_EPS, msg='Ball constraint violated')
            self.assertTrue(np.max(x - a) >= 0.0, msg='Lower bound violated')
            self.assertTrue(np.max(x - b) <= 0.0, msg='Upper bound violated')


class TestLinear2WithZeros(unittest.TestCase):
    def runTest(self):
        g = np.array([0.0, -1.0])
        a = np.array([-2.0, -2.0])
        b = np.array([1.0, 2.0])
        delta = 5.0
        for H in [None, np.zeros((len(g), len(g)))]:
            x = trustregion.solve(g, H, delta, sl=a, su=b)
            xtrue = np.array([0.0, 2.0])
            self.assertTrue(np.max(np.abs(x - xtrue)) < 1e-10, msg='Wrong step')
            self.assertTrue(np.linalg.norm(x) <= delta+BALL_EPS, msg='Ball constraint violated')
            self.assertTrue(np.max(x - a) >= 0.0, msg='Lower bound violated')
            self.assertTrue(np.max(x - b) <= 0.0, msg='Upper bound violated')


class TestLinear2WithAlmostZeros(unittest.TestCase):
    def runTest(self):
        g = np.array([1e-15, -1.0])
        a = np.array([-2.0, -2.0])
        b = np.array([1.0, 2.0])
        delta = 5.0
        x = trustregion.solve(g, None, delta, sl=a, su=b)
        xtrue = np.array([0.0, 2.0])
        self.assertTrue(np.max(np.abs(x - xtrue)) < 1e-10, msg='Wrong step')
        self.assertTrue(np.linalg.norm(x) <= delta+BALL_EPS, msg='Ball constraint violated')
        self.assertTrue(np.max(x - a) >= 0.0, msg='Lower bound violated')
        self.assertTrue(np.max(x - b) <= 0.0, msg='Upper bound violated')


class TestLinear2WithAlmostZeros2(unittest.TestCase):
    def runTest(self):
        g = np.array([1e-15, 0.0])
        a = np.array([-2.0, -2.0])
        b = np.array([1.0, 2.0])
        delta = 5.0
        for H in [None, np.zeros((len(g), len(g)))]:
            x = trustregion.solve(g, H, delta, sl=a, su=b)
            # Since objective is essentially zero, will accept any x within the defined region
            self.assertTrue(np.linalg.norm(x) <= delta+BALL_EPS, msg='Ball constraint violated')
            self.assertTrue(np.max(x - a) >= 0.0, msg='Lower bound violated')
            self.assertTrue(np.max(x - b) <= 0.0, msg='Upper bound violated')


class TestZero(unittest.TestCase):
    def runTest(self):
        n = 5
        g = np.ones((n,))
        sl = np.zeros((n,))
        su = np.zeros((n,))
        delta = 1.0
        for H in [None, np.zeros((len(g), len(g)))]:
            x = trustregion.solve(g, H, delta, sl=sl, su=su)
            self.assertAlmostEqual(np.linalg.norm(x), 0.0, msg='Nonzero step')


class TestZero2(unittest.TestCase):
    def runTest(self):
        n = 5
        g = np.ones((n,))
        sl = -np.ones((n,))
        su = np.ones((n,))
        delta = 0.0
        for H in [None, np.zeros((len(g), len(g)))]:
            x = trustregion.solve(g, H, delta, sl=sl, su=su)
            self.assertAlmostEqual(np.linalg.norm(x), 0.0, msg='Nonzero step')


class TestOneSided(unittest.TestCase):
    def runTest(self):
        g = np.array([1.0, 1.0])
        a = np.array([0.0, 0.0])
        b = np.array([1.0, 2.0])
        delta = 5.0
        for H in [None, np.zeros((len(g), len(g)))]:
            x = trustregion.solve(g, None, delta, sl=a, su=b)
            xtrue = np.array([0.0, 0.0])
            self.assertTrue(np.max(np.abs(x - xtrue)) < 1e-10, msg='Wrong step')
            self.assertTrue(np.linalg.norm(x) <= delta + BALL_EPS, msg='Ball constraint violated')
            self.assertTrue(np.max(x - a) >= 0.0, msg='Lower bound violated')
            self.assertTrue(np.max(x - b) <= 0.0, msg='Upper bound violated')


class TestOneSided2(unittest.TestCase):
    def runTest(self):
        g = np.array([1.0, -1.0])
        a = np.array([0.0, -1.0])
        b = np.array([1.0, 0.0])
        delta = 5.0
        for H in [None, np.zeros((len(g), len(g)))]:
            x = trustregion.solve(g, None, delta, sl=a, su=b)
            xtrue = np.array([0.0, 0.0])
            self.assertTrue(np.max(np.abs(x - xtrue)) < 1e-10, msg='Wrong step')
            self.assertTrue(np.linalg.norm(x) <= delta + BALL_EPS, msg='Ball constraint violated')
            self.assertTrue(np.max(x - a) >= 0.0, msg='Lower bound violated')
            self.assertTrue(np.max(x - b) <= 0.0, msg='Upper bound violated')


class TestAnalytic(unittest.TestCase):
    def runTest(self):
        g = np.array([1.0, -1.0])
        delta = 5.0
        for H in [None, np.zeros((len(g), len(g)))]:
            x = trustregion.solve(g, H, delta)
            xtrue = delta * g / np.linalg.norm(g)
            self.assertTrue(np.max(np.abs(x - xtrue)) < 1e-10, msg='Wrong step')
            self.assertTrue(np.linalg.norm(x) <= delta + BALL_EPS, msg='Ball constraint violated')


class TestTypeErrors(unittest.TestCase):
    def runTest(self):
        g = np.array([1.0, 0.0, 1.0])
        H = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
        Delta = 2.0
        sl = np.array([-0.5, -10.0, -10.0])
        su = np.array([10.0, 10.0, 10.0])
        self.assertRaises(ValueError, trustregion.solve, 'abc', H, Delta)
        self.assertRaises(ValueError, trustregion.solve, g, 'abc', Delta)
        self.assertRaises(ValueError, trustregion.solve, g, H, None)
        self.assertRaises(ValueError, trustregion.solve, g, H, Delta, sl='abc', su=su)
        self.assertRaises(ValueError, trustregion.solve, g, H, Delta, sl=sl, su='abc')
        self.assertRaises(ValueError, trustregion.solve, g, H, Delta, verbose_output='abc')
        # Need both or neither box bounds
        self.assertRaises(ValueError, trustregion.solve, g, H, Delta, sl=None, su=su)
        self.assertRaises(ValueError, trustregion.solve, g, H, Delta, sl=sl, su=None)


class TestDimensionErrors(unittest.TestCase):
    def runTest(self):
        # Where all vector/matrix inputs have wrong dimension or shape
        g = np.array([1.0, 0.0])
        H1 = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
        H2 = np.eye(2)
        H3 = np.zeros((2,3))
        H4 = np.zeros((3, 2))
        sl1 = -np.arange(5)
        sl2 = -np.ones((2,))
        su1 = np.arange(5)
        su2 = np.ones((2,))
        delta = 1.0
        self.assertRaises(ValueError, trustregion.solve, g, H1, delta)
        self.assertRaises(ValueError, trustregion.solve, g, H3, delta)
        self.assertRaises(ValueError, trustregion.solve, g, H4, delta)
        self.assertRaises(ValueError, trustregion.solve, g, H2, delta, sl=sl1, su=su2)
        self.assertRaises(ValueError, trustregion.solve, g, H2, delta, sl=sl2, su=su1)


class TestNumberErrors(unittest.TestCase):
    def runTest(self):
        # Where values are dodgy
        g = np.array([1.0, 0.0])
        H1 = np.eye(2)
        H2 = np.eye(2); H2[1,0] = -1.0  # non-symmetric
        H3 = np.eye(2, dtype=np.complex); H3[0,0] = 1.0+1j  # complex
        delta1 = -1.0
        delta2 = 1.0
        sl1 = np.array([-1.0, 0.01])
        sl2 = -np.ones((2,))
        su1 = np.array([-0.001, 1.0])
        su2 = np.ones((2,))
        self.assertRaises(ValueError, trustregion.solve, g, H1, delta1)
        self.assertRaises(ValueError, trustregion.solve, g, H2, delta2)
        self.assertRaises(ValueError, trustregion.solve, g, H3, delta2)
        self.assertRaises(ValueError, trustregion.solve, g, H1, delta2, sl=sl1, su=su2)
        self.assertRaises(ValueError, trustregion.solve, g, H1, delta2, sl=sl2, su=su1)