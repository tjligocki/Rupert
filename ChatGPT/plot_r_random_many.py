#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Define function to construct continuous orthonormal basis given normal n
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

# Function to compute f(theta) given n
def f_theta(theta, u0, v0):
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

  return max(vals)

# Random seed based on current time
import time
np.random.seed(int(time.time()) % (2**32 - 1))

# Function to plot multiple random normals
def plot_multiple_random_normals(num_normals=3):
  theta_vals = np.linspace(0, 2*np.pi, 500)
  plt.figure(figsize=(12,8))

  for i in range(num_normals):
    n_rand = np.random.randn(3)
    n_rand = n_rand / np.linalg.norm(n_rand)
    u0_rand, v0_rand = construct_basis(n_rand)

    f_vals = np.array([f_theta(theta, u0_rand, v0_rand) for theta in theta_vals])
    r_vals = 1 / f_vals

    max_r = np.max(r_vals)
    theta_max_r = theta_vals[np.argmax(r_vals)]

    print(f"Normal {i+1}: n = {n_rand}")
    print(f"  Maximum r(θ) ≈ {max_r:.10f} at θ ≈ {theta_max_r:.10f} radians\n")

    plt.plot(theta_vals, r_vals, label=f'Normal {i+1} max r ≈ {max_r:.5f}')

  plt.title(f'r(θ) for {num_normals} Random Normal Vectors')
  plt.xlabel('θ (radians)')
  plt.ylabel('r(θ)')
  plt.ylim(0, 1.5)
  plt.grid(True)
  plt.legend(loc='upper right')
  plt.tight_layout()
  plt.show()

# Example usage
plot_multiple_random_normals(num_normals=20)

