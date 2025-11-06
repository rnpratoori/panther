# -*- coding: utf-8 -*-
"""
Standalone Python script to calculate the maximum number of non-overlapping
spheres that can be placed in a 3D grid within a rectangular prism with a
defined buffer.

This script calculates and prints the (x, y, z) coordinates for the center
of each possible sphere for both gridded and staggered arrangements.
"""
import math

# -----------------------------------------------------------------------------
#                               INPUT PARAMETERS
# -----------------------------------------------------------------------------
# Arrangement type: 'grid' or 'staggered'
ARRANGEMENT = 'grid'  # Change to 'staggered' for that arrangement

# Dimensions of the containing rectangular prism
RECT_WIDTH = 2.0
RECT_HEIGHT = 1.0
RECT_DEPTH = 1.0

# Buffer space from the edge of the prism
BUFFER = 0.1

# Properties of the spheres
SPHERE_RADIUS = 0.1
EDGE_SPACING = 0.4

# -----------------------------------------------------------------------------
#                                 CALCULATIONS
# -----------------------------------------------------------------------------

# Calculate the center-to-center distance from the edge spacing
center_to_center_spacing = EDGE_SPACING + 2 * SPHERE_RADIUS

# Step 1: Define the Effective Placement Area
effective_width = RECT_WIDTH - 2 * BUFFER
effective_height = RECT_HEIGHT - 2 * BUFFER
effective_depth = RECT_DEPTH - 2 * BUFFER

# Step 2: Determine the Starting Point
x_start = BUFFER + SPHERE_RADIUS
y_start = RECT_HEIGHT - BUFFER - SPHERE_RADIUS
z_start = RECT_DEPTH - BUFFER - SPHERE_RADIUS

# Step 3: Calculate the Number of Spheres that Fit
span_x = effective_width - 2 * SPHERE_RADIUS
span_y = effective_height - 2 * SPHERE_RADIUS
span_z = effective_depth - 2 * SPHERE_RADIUS

if span_x < 0 or span_y < 0 or span_z < 0:
    print("Error: The sphere radius is too large to fit within the effective area.")
    exit()

num_x = 1 + math.floor(span_x / center_to_center_spacing)
num_y = 1 + math.floor(span_y / center_to_center_spacing)
num_z = 1 + math.floor(span_z / center_to_center_spacing)

# Step 4: Generate and Store the Center Coordinates
all_centers = []
horizontal_offset = center_to_center_spacing / 2.0
right_boundary = RECT_WIDTH - BUFFER
top_boundary = RECT_HEIGHT - BUFFER 

if ARRANGEMENT == 'grid':
    for i in range(num_z):
        current_z = z_start - i * center_to_center_spacing
        for j in range(num_y):
            current_y = y_start - j * center_to_center_spacing
            for k in range(num_x):
                current_x = x_start + k * center_to_center_spacing
                all_centers.append((current_x, current_y, current_z))

elif ARRANGEMENT == 'staggered':
    for i in range(num_z):
        current_z = z_start - i * center_to_center_spacing
        start_x_for_layer = x_start
        start_y_for_layer = y_start

        # Offset every other layer (in the z-direction)
        if i % 2 == 1:
            start_x_for_layer += horizontal_offset
            start_y_for_layer -= horizontal_offset

        for j in range(num_y):
            current_y = start_y_for_layer - j * center_to_center_spacing
            for k in range(num_x):
                current_x = start_x_for_layer + k * center_to_center_spacing
                
                # Ensure the staggered sphere does not extend beyond the boundaries
                if (current_x + SPHERE_RADIUS) <= right_boundary and \
                   (current_y + SPHERE_RADIUS) <= top_boundary:
                    all_centers.append((current_x, current_y, current_z))

# -----------------------------------------------------------------------------
#                                   OUTPUT
# -----------------------------------------------------------------------------
print(f"--- Sphere Placement Calculation Results ({ARRANGEMENT.title()} Arrangement) ---")
print(f"Prism Dimensions: {RECT_WIDTH} x {RECT_HEIGHT} x {RECT_DEPTH}")
print(f"Sphere Radius: {SPHERE_RADIUS}, Edge-to-Edge Spacing: {EDGE_SPACING}")
print(f"Buffer: {BUFFER}\n")

if all_centers:
    print(f"Found a total of {len(all_centers)} possible sphere centers.")
    print("Coordinates (x, y, z):")
    for index, center in enumerate(all_centers):
        print(f"  Sphere {index + 1:>3}: ({center[0]:.4f}, {center[1]:.4f}, {center[2]:.4f})")
else:
    print("No complete spheres could be placed with the given parameters.")

print("\n--- End of Script ---")
