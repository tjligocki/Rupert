#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Define the cube (edges from -1 to 1)
cube_vertices = np.array([
  [-1, -1, -1],
  [-1, -1,  1],
  [-1,  1, -1],
  [-1,  1,  1],
  [ 1, -1, -1],
  [ 1, -1,  1],
  [ 1,  1, -1],
  [ 1,  1,  1]
])

cube_faces = [
  [0, 1, 3, 2],
  [4, 5, 7, 6],
  [0, 1, 5, 4],
  [2, 3, 7, 6],
  [0, 2, 6, 4],
  [1, 3, 7, 5]
]

# Define the trial u and v vectors
u = np.array([1/np.sqrt(2), -1/np.sqrt(2), 0])
v = np.array([1/np.sqrt(6), 1/np.sqrt(6), -2/np.sqrt(6)])

# Define function to generate square vertices in proper order
def generate_square(center, s, u, v):
  vertices = []
  for sign1, sign2 in [(-1, -1), (1, -1), (1, 1), (-1, 1)]:
      vertex = center + (s/2)*(sign1*u + sign2*v)
      vertices.append(vertex)
  return np.array(vertices)

# Center at origin
center = np.array([0,0,0])

# Unscaled square (s=1)
s_unscaled = 1
square_unscaled = generate_square(center, s_unscaled, u, v)

# Scaled square (s computed)
s_scaled = 2 / np.max([
  np.max([np.abs(u[0] + v[0]), np.abs(u[0] - v[0])]),
  np.max([np.abs(u[1] + v[1]), np.abs(u[1] - v[1])]),
  np.max([np.abs(u[2] + v[2]), np.abs(u[2] - v[2])])
])
square_scaled = generate_square(center, s_scaled, u, v)

# Create plot
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d', proj_type='ortho')

# Plot cube faces
for face in cube_faces:
  face_vertices = cube_vertices[face]
  poly3d = [[tuple(point) for point in face_vertices]]
  ax.add_collection3d(Poly3DCollection(poly3d, facecolors='lightgray', linewidths=1, edgecolors='black', alpha=0.3))

# Plot unscaled square
unscaled_faces = [[tuple(point) for point in square_unscaled]]
ax.add_collection3d(Poly3DCollection(unscaled_faces, facecolors='blue', alpha=0.5))
ax.scatter(square_unscaled[:,0], square_unscaled[:,1], square_unscaled[:,2], color='blue', s=20)

# Plot scaled square
scaled_faces = [[tuple(point) for point in square_scaled]]
ax.add_collection3d(Poly3DCollection(scaled_faces, facecolors='orange', alpha=0.6))
ax.scatter(square_scaled[:,0], square_scaled[:,1], square_scaled[:,2], color='orange', s=20)

# Plot cube corners
ax.scatter(cube_vertices[:,0], cube_vertices[:,1], cube_vertices[:,2], color='gray', s=20)

# Set limits and labels
ax.set_xlim([-1.5, 1.5])
ax.set_ylim([-1.5, 1.5])
ax.set_zlim([-1.5, 1.5])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Unscaled and Scaled Square Inside Cube (Orthographic Projection)')
ax.view_init(elev=20, azim=30)

plt.tight_layout()
plt.show()

