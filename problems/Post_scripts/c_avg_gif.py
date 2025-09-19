import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter

def animate_averages(df, title, filename, fps=20, as_gif=False):
    fig, ax = plt.subplots()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Time')
    ax.set_ylabel('Average value')
    ax.set_title(title)

    line1, = ax.plot([], [], 'r-', label='avg_c1')
    line2, = ax.plot([], [], 'g-', label='avg_c2')
    line3, = ax.plot([], [], 'b-', label='avg_c3')
    ax.legend()

    def init():
        line1.set_data([], [])
        line2.set_data([], [])
        line3.set_data([], [])
        return line1, line2, line3

    def update(frame):
        x = df['time'][:frame+1]
        y1 = df['avg_c1'][:frame+1]
        y2 = df['avg_c2'][:frame+1]
        y3 = df['avg_c3'][:frame+1]
        line1.set_data(x, y1)
        line2.set_data(x, y2)
        line3.set_data(x, y3)
        return line1, line2, line3

    anim = FuncAnimation(fig, update, frames=len(df), init_func=init, blit=True, interval=1000/fps)
    if as_gif:
        anim.save(filename, writer=PillowWriter(fps=fps))
    else:
        anim.save(filename, writer=FFMpegWriter(fps=fps))
    plt.close(fig)
    print(f"Saved animation to {filename}")

# Example usage after your Excel/plot saving loop:
# animate_averages(df, "Averages vs Time (block0)", "block0_anim.mp4")