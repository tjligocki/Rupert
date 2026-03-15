#!/usr/bin/env python3

import numpy as np
from scipy.optimize import minimize
import time
import argparse

# Argument parser for method and tolerance
parser = argparse.ArgumentParser(description="Optimize r(theta, n)")
parser.add_argument('--method', type=str, default='Powell', choices=['Nelder-Mead', 'Powell', 'COBYLA'], help='Optimization method')
parser.add_argument('--ftol', type=float, default=1e-10, help='Function tolerance for convergence')
args = parser.parse_args()

# Construct continuous orthonormal basis given normal n
def construct_basis(n):
  n_x, n_y, n_z = n
  if abs(n_x) <= abs(n_y) and abs(n_x) <= abs(n_z):
    u0 = np.array([0, -n_z, n_y])
  elif abs(n_y) <= abs(n_x) and abs(n_y) <= abs(n_z):
    u0 = np.array([-n_z, 0, n_x])
  else:
    u0 = np.array([-n_y, n_x, 0])
  u0 = u0 / np.linalg.norm(u0)
  v0 = np.cross(n, u0)
  v0 = v0 / np.linalg.norm(v0)
  return u0, v0

# Compute r(theta, n)
def r_theta_n(theta, n):
  u0, v0 = construct_basis(n)
  cos_t = np.cos(theta)
  sin_t = np.sin(theta)

  u = cos_t * u0 + sin_t * v0
  v = -sin_t * u0 + cos_t * v0

  vals = []
  for j in range(3):
    u_j = u[j]
    v_j = v[j]
    vals.append(abs(u_j + v_j))
    vals.append(abs(u_j - v_j))

  return 1 / max(vals)

# Wrapper function to optimize (negative r because we minimize)
def objective(x):
  theta = x[0]
  n = x[1:]
  n = n / np.linalg.norm(n)  # ensure unit vector
  return -r_theta_n(theta, n)

# Use a different random seed each time
seed = int(time.time()) % (2**32 - 1)
np.random.seed(seed)

# Generate random initial theta and unit normal
theta_init = np.random.uniform(0, 2*np.pi)
n_init = np.random.randn(3)
n_init = n_init / np.linalg.norm(n_init)
x0 = np.concatenate(([theta_init], n_init))

print()

# Print the seed and initial values
print(f"Random seed: {seed}")
print(f"Initial theta: {theta_init:.10f} radians")
print(f"Initial n: {n_init}")

# Perform optimization
res = minimize(objective, x0, method=args.method, options={'maxiter': 5000, 'disp': True, 'ftol': args.ftol})

# Extract results
theta_opt = res.x[0]
n_opt = res.x[1:]
n_opt = n_opt / np.linalg.norm(n_opt)
r_opt = -res.fun

# Print result
print(f"\nOptimization method: {args.method}")
print(f"Function tolerance: {args.ftol}")
print(f"Maximum r ≈ {r_opt:.10f}")
print(f"Optimal theta ≈ {theta_opt:.10f} radians")
print(f"Optimal n ≈ {n_opt}")

print()
