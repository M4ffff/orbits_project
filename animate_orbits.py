import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import argparse


parser = argparse.ArgumentParser(description='Produces an animated plot of the planets in motion.')
parser.add_argument('--filename', default='orbit_positions.csv', help='Name of the csv file to animate')
parser.add_argument('--frames', default=100, help='Threshold number of frames')

args = parser.parse_args()

threshold_frames = int(args.frames)

# Define our plotting area
fig, ax = plt.subplots(1, 2, figsize=(8, 4))
ax[0].set_aspect("equal")
ax[1].set_aspect("equal")

data_file = args.filename

positions = pd.read_csv(data_file)

num_positions = len(positions)

if num_positions > threshold_frames:
    print(f"WARNING: More than threshold number of frames: {threshold_frames}. Reducing and averaging for plotting.")
    positions["bins"] = pd.cut(positions.index.values, bins=threshold_frames)
    positions = positions.groupby("bins").mean()

num_bodies = int(len(positions.columns)/2)
num_frames = len(positions)

print(f"Number of frames: {num_frames}")

# Define time for each frame in the animation
# frames_per_second = 20
# interval = 1000/frames_per_second # time for each frame in milliseconds
interval = 100 - 25/num_frames  # time for each frame in milliseconds

def animate(i):
    ax[0].clear() 
    ax[1].clear() 

    print(f"Frame: {i}")

    nn = 3
    if num_bodies < 3:
        nn = 2
    for j in range(1, nn+1):
        x_label = f"x{j}"
        y_label = f"y{j}"
        x = positions[x_label]
        y = positions[y_label]

        if j == 1:
            x_min = np.min(x)
            y_min = np.min(y)
            x_max = np.max(x)
            y_max = np.max(y)
        else:
            if np.min(x) < x_min:
                x_min = np.min(x)
            if np.min(y) < y_min:
                y_min = np.min(y)
            if np.max(x) > x_max:
                x_max = np.max(x)
            if np.max(y) > y_max:
                y_max = np.max(y)

        x_current = x.iloc[:i + 1]
        y_current = y.iloc[:i + 1]

        ax[0].plot(x_current, y_current, color=f"C{j}")
        ax[0].scatter(x.iloc[0], y.iloc[0], marker='x', color=f"C{j}")
        ax[0].scatter(x_current.iloc[-1], y_current.iloc[-1], marker='o', color=f"C{j}")

    x_min = x_min - np.abs(x_max) * 0.1
    y_min = y_min - np.abs(y_max) * 0.1
    x_max = x_max + np.abs(x_max) * 0.1
    y_max = y_max + np.abs(y_max) * 0.1

    ax[0].set_xlim(x_min, x_max)
    ax[0].set_ylim(y_min, y_max)

    for j in range(2, num_bodies + 1):
        x_label = f"x{j}"
        y_label = f"y{j}"
        x = positions[x_label]
        y = positions[y_label]

        x_current = x.iloc[:i + 1]
        y_current = y.iloc[:i + 1]
        if j == 2:
            x_lim = x_current.iloc[-1]
            y_lim = y_current.iloc[-1]

        ax[1].plot(x_current, y_current, color=f"C{j}")
        ax[1].scatter(x_current.iloc[-1], y_current.iloc[-1], marker='o', color=f"C{j}")

    ax[1].set_xlim(x_lim - 1e9, x_lim + 1e9)
    ax[1].set_ylim(y_lim - 1e9, y_lim + 1e9)

anim = animation.FuncAnimation(fig, animate, frames=num_frames, interval=interval)
anim.save("orbits.gif")