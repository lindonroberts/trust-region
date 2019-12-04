"""
Examples of usage of trustregion.solve
"""
# Ensure compatibility with Python 2
from __future__ import absolute_import, division, print_function, unicode_literals


import numpy as np
import trustregion

# Example 1: Linear objective, no box constraints
# min_{s} g*s  subject to ||s||_2 <= delta
g = np.array([1.0, 2.0, 3.0])
delta = 1.0
s = trustregion.solve(g, None, delta)  # H=None or H=np.zeros(...) both valid
print("Example 1, s =", s)

# Example 2: Linear objective, with box constraints
# min_{s} g*s  subject to ||s||_2 <= delta and sl <= s <= su
g = np.array([-1.0, 2.0, 3.0])
delta = 0.5
sl = np.array([-10.0, -0.5, -0.3])  # lower bounds, need sl <= 0
su = np.array([0.5, 10.0, 1.0])  # upper bounds, need su >= 0
s = trustregion.solve(g, None, delta, sl=sl, su=su)  # H=None or H=np.zeros(...) both valid
print("Example 2, s =", s)

# Example 3: Quadratic objective, no box constraints
# min_{s} g*s + 0.5*s*H*s  subject to ||s||_2 <= delta
g = np.array([1.0, 2.0, 3.0])
H = np.eye(3)  # must be real, symmetric
delta = 1.0
s = trustregion.solve(g, H, delta)
print("Example 3, s =", s)

# Example 4: Quadratic objective, with box constraints
# min_{s} g*s + 0.5*s*H*s  subject to ||s||_2 <= delta and sl <= s <= su
g = np.array([1.0, 2.0, 3.0])
H = np.eye(3)  # must be real, symmetric
delta = 1.0
sl = np.array([-10.0, -0.5, -0.3])  # lower bounds, need sl <= 0
su = np.array([0.5, 10.0, 1.0])  # upper bounds, need su >= 0
s = trustregion.solve(g, H, delta, sl=sl, su=su)
print("Example 4, s =", s)

# Verbose output (use example 4)
# s = solution to trust-region subproblem
# gnew = gradient of objective at solution s (i.e. gnew = g + H*s)
# crvmin = smallest curvature of H found, except for:
#  - crvmin = 0 if ||s||=delta
#  - crvmin = -1 if box constraints active in all coordinates
s, gnew, crvmin = trustregion.solve(g, H, delta, sl=sl, su=su, verbose_output=True)
print("Example 4 (verbose output):")
print("s =", s)
print("gnew =", gnew)
print("crvmin =", crvmin)
