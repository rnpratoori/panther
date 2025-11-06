# -*- coding: utf-8 -*-
"""
Standalone Python script to calculate the maximum number of non-overlapping
circles that can be placed in a grid within a rectangle with a defined buffer.

This script calculates and prints the (x, y) coordinates for the center
of each possible circle.
"""
import math

# -----------------------------------------------------------------------------
#                               INPUT PARAMETERS
# -----------------------------------------------------------------------------
# You can change these values to see how the results change.

# Dimensions of the containing rectangle
RECT_WIDTH = 2.0
RECT_HEIGHT = 1.0

# Buffer space from the edge of the rectangle
BUFFER = 0.1

# Properties of the circles
CIRCLE_RADIUS = 0.05
# CHANGED: This is now the spacing between the *edges* of the circles.
EDGE_SPACING = 0.2

# --- ASCII Diagram of the Setup ---
#
#    +----------------------------------------------------+ RECT_WIDTH (2.0)
#    |  BUFFER (0.1)                                      |
#    |  +----------------------------------------------+  |
#    |  |               EFFECTIVE AREA               |  | RECT_HEIGHT (1.0)
#    |  |                                              |  |
#    |  |      (o) spacing (o)----(o)                  |  |
#    |  |       |           |      |                   |  |
#    |  |      (o)----(o)----(o)                  |  |
#    |  |       .           .      .                   |  |
#    |  |       .           .      .                   |  |
#    |  +----------------------------------------------+  |
#    |                                                    |
#    +----------------------------------------------------+
#

# -----------------------------------------------------------------------------
#                                 CALCULATIONS
# -----------------------------------------------------------------------------

# NEW: Calculate the center-to-center distance from the edge spacing
center_to_center_spacing = EDGE_SPACING + 2 * CIRCLE_RADIUS

# Step 1: Define the Effective Placement Area
# This is the inner rectangle where circles can actually be placed.
effective_width = RECT_WIDTH - 2 * BUFFER
effective_height = RECT_HEIGHT - 2 * BUFFER

# Step 2: Determine the Starting Point
# The center of the first row of circles is in the top-left of the effective area,
# offset downwards by its own radius.
x_start = BUFFER + CIRCLE_RADIUS
y_start = RECT_HEIGHT - BUFFER - CIRCLE_RADIUS

# Step 3: Calculate the Number of Circles that Fit
# First, find the total distance the centers can span.
span_x = effective_width - 2 * CIRCLE_RADIUS
span_y = effective_height - 2 * CIRCLE_RADIUS

# Check if any circles can fit at all.
if span_x < 0 or span_y < 0:
    print("Error: The circle radius is too large to fit within the effective area.")
    # Exit the script cleanly if no circles fit
    exit()

# Calculate the maximum number of circles that could fit in any given row.
# The actual number in staggered rows might be one less.
num_x = 1 + math.floor(span_x / center_to_center_spacing)
num_y = 1 + math.floor(span_y / center_to_center_spacing)

# Step 4: Generate and Store the Center Coordinates for a STAGGERED grid
all_centers = []
horizontal_offset = center_to_center_spacing / 2.0
right_boundary = RECT_WIDTH - BUFFER

for i in range(num_y):
    # Subtract spacing to place circles downwards from the top
    current_y = y_start - i * center_to_center_spacing
    
    start_x_for_row = x_start
    # Check if this is an even row (2nd, 4th, etc., counting from the top).
    # In a 0-indexed loop, this corresponds to odd indices (1, 3, 5...).
    if i % 2 == 1:
        start_x_for_row += horizontal_offset

    for j in range(num_x):
        current_x = start_x_for_row + j * center_to_center_spacing
        
        # Ensure the staggered circle does not extend beyond the right buffer
        if (current_x + CIRCLE_RADIUS) <= right_boundary:
            all_centers.append((current_x, current_y))


# -----------------------------------------------------------------------------
#                                   OUTPUT
# -----------------------------------------------------------------------------
print("--- Circle Placement Calculation Results (Staggered Grid) ---")
print(f"Rectangle Dimensions: {RECT_WIDTH} x {RECT_HEIGHT}")
print(f"Circle Radius: {CIRCLE_RADIUS}, Edge-to-Edge Spacing: {EDGE_SPACING}")
print(f"Buffer: {BUFFER}\n")

if all_centers:
    print(f"Found a total of {len(all_centers)} possible circle centers.")
    print("Coordinates (x, y):")
    for index, center in enumerate(all_centers):
        # Format the output for readability, rounding to 4 decimal places
        print(f"  Circle {index + 1:>2}:  ({center[0]:.4f}, {center[1]:.4f})")
else:
    # This message is a fallback, though the earlier check should catch this.
    print("No complete circles could be placed with the given parameters.")

print("\n--- End of Script ---")



