#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Define reference orthonormal basis in the plane orthogonal to n = (1,1,1)
u0 = np.array([1/np.sqrt(2), -1/np.sqrt(2), 0])
v0 = np.array([1/np.sqrt(6), 1/np.sqrt(6), -2/np.sqrt(6)])

# Function to compute f(theta)
def f_theta(theta):
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

# Generate theta values and evaluate f(theta) and r(theta)
theta_vals = np.linspace(0, 2*np.pi, 500)
f_vals = np.array([f_theta(theta) for theta in theta_vals])
r_vals = 1 / f_vals

# Find maximum r(theta)
max_r = np.max(r_vals)
theta_max_r = theta_vals[np.argmax(r_vals)]

# Plot the result
plt.figure(figsize=(10,6))
plt.plot(theta_vals, r_vals, label='r(θ) = 1 / f(θ)')
plt.axhline(max_r, color='red', linestyle='--', label=f'Max r(θ) ≈ {max_r:.6f}')
plt.axvline(theta_max_r, color='green', linestyle='--', label=f'θ max ≈ {theta_max_r:.4f} rad')
plt.title('Maximum r(θ) for Square in Plane Perpendicular to (1,1,1)')
plt.xlabel('θ (radians)')
plt.ylabel('r(θ)')
plt.grid(True)
plt.legend(loc='center right')
plt.tight_layout()
plt.show()

# Print the result
print(f"Maximum r(θ) ≈ {max_r:.10f} at θ ≈ {theta_max_r:.10f} radians")

