"""Calculate circle layouts and area coverage for 2D void-spacing studies."""

import math
from itertools import product

# ----------------------------------------------------------------------
# Domain dimensions and buffer
RECT_WIDTH = 2.0
RECT_HEIGHT = 1.0
BUFFER = 0.1

# Radii and spacings to test
RADII = [0.05, 0.10, 0.15]
EDGE_SPACINGS = [0.2, 0.4]

# ----------------------------------------------------------------------
def compute_circle_centers(radius, edge_spacing):
    """Return list of circle centers for given radius and edge spacing."""
    center_to_center_spacing = edge_spacing + 2 * radius

    # Effective placement area
    effective_width = RECT_WIDTH - 2 * BUFFER
    effective_height = RECT_HEIGHT - 2 * BUFFER

    # Starting position (top-left circle)
    x_start = BUFFER + radius
    y_start = RECT_HEIGHT - BUFFER - radius

    # Span available for centers
    span_x = effective_width - 2 * radius
    span_y = effective_height - 2 * radius

    # How many circles in x and y directions
    if span_x < 0 or span_y < 0:
        return []  # no circles fit
    num_x = 1 + math.floor(span_x / center_to_center_spacing)
    num_y = 1 + math.floor(span_y / center_to_center_spacing)

    centers = []
    for i in range(num_y):
        current_y = y_start - i * center_to_center_spacing
        for j in range(num_x):
            current_x = x_start + j * center_to_center_spacing
            centers.append((current_x, current_y))
    return centers

# ----------------------------------------------------------------------
# Run all combinations and compute percentages
rect_area = RECT_WIDTH * RECT_HEIGHT

print("--- Circle Coverage Results ---")
for radius, spacing in product(RADII, EDGE_SPACINGS):
    centers = compute_circle_centers(radius, spacing)
    circle_area = len(centers) * math.pi * radius**2
    percentage = 100.0 * circle_area / rect_area
    print(f"Radius={radius:.2f}, Spacing={spacing:.2f} → "
          f"{len(centers)} circles, "
          f"Coverage={percentage:.2f}%")
