"""Generate an RGB Maxwell triangle legend for ternary composition plots."""

import numpy as np
import matplotlib.pyplot as plt

def plot_maxwell_triangle_legend(filename='maxwell_triangle_legend.png', resolution=400):
    # Create a grid of barycentric coordinates
    n = resolution
    img = np.ones((n, n, 3))
    for i in range(n):
        for j in range(n):
            # Barycentric coordinates
            c1 = i / (n - 1)
            c2 = j / (n - 1)
            c3 = 1 - c1 - c2
            if c3 < 0 or c3 > 1:
                img[j, i, :] = 1  # white outside triangle
            else:
                img[j, i, :] = [c1, c2, c3]  # RGB mapping

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(img, origin='lower', extent=(0, 1, 0, 1))

    # Draw triangle borders
    triangle = np.array([[0, 0], [1, 0], [0, 1], [0, 0]])
    ax.plot(triangle[:, 0], triangle[:, 1], color='k', lw=2)

    # Annotate corners
    ax.text(-0.05, -0.05, 'c3 (Blue)', color='blue', fontsize=12, ha='right', va='top')
    ax.text(1.05, -0.05, 'c1 (Red)', color='red', fontsize=12, ha='left', va='top')
    ax.text(-0.05, 1.05, 'c2 (Green)', color='green', fontsize=12, ha='right', va='bottom')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, 1.1)
    ax.set_aspect('equal')
    ax.set_title('Maxwell Triangle (RGB Legend)')
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"Saved Maxwell triangle legend as {filename}")

if __name__ == "__main__":
    plot_maxwell_triangle_legend()
