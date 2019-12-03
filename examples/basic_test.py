import numpy as np
import trustregion


n = 3
g = np.array([1.0, 0.0, 1.0])
H = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]])
Delta = 2.0
xopt = np.ones((n,))  # trying nonzero (since bounds inactive)
sl = -1e20 * np.ones((n,))
su = 1e20 * np.ones((n,))
d = trustregion.solve(g, H, Delta, sl=sl, su=su)
true_d = np.array([-1.0, 0.0, -0.5])
print("s =", d)
print("True soln =", true_d)

Delta = 5.0 / 12.0
d = trustregion.solve(g, H, Delta, sl=sl, su=su)
true_d = np.array([-1.0 / 3.0, 0.0, -0.25])
print("s =", d)
print("True soln =", true_d)
print("Done")

g = np.array([1.0, -1.1])
a = np.array([-2.0, -2.0])
b = np.array([1.0, 2.0])
delta = 2.8
d = trustregion.solve(g, None, delta, sl=a, su=b)
print("trustregion =", d)
print("Done")