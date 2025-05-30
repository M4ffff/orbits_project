import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import argparse


parser = argparse.ArgumentParser(description='Produces an animated plot of the planets in motion.')
parser.add_argument('--filename', default='1.00_sun__1.0_moon.csv', help='Name of the standard csv file to animate')
parser.add_argument('--extra_filename', default='1.00_sun__1.0_moon.csv', help='Name of the csv file to animate')
parser.add_argument('--frames', default=100, help='Threshold number of frames')

args = parser.parse_args()

threshold_frames = int(args.frames)

# Define our plotting area
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].set_aspect("equal")
ax[1].set_aspect("equal")

standard_data_file = args.filename
extra_data_file = args.extra_filename
print(extra_data_file)

standard_positions = pd.read_csv(standard_data_file)
extra_positions = pd.read_csv(extra_data_file)

num_standard_positions = len(standard_positions)

if num_standard_positions > threshold_frames:
    print(f"WARNING: More than threshold number of frames: {threshold_frames}. Reducing and averaging for plotting.")
    standard_positions["bins"] = pd.cut(standard_positions.index.values, bins=threshold_frames)
    standard_positions = standard_positions.groupby("bins").mean()
    extra_positions["bins"] = pd.cut(extra_positions.index.values, bins=threshold_frames)
    extra_positions = extra_positions.groupby("bins").mean()

num_bodies = int(len(standard_positions.columns)/2)
num_frames = len(standard_positions)

print(f"Number of frames: {num_frames}")

# Define time for each frame in the animation
# frames_per_second = 20
# interval = 1000/frames_per_second # time for each frame in milliseconds
interval = 100 - 25/num_frames  # time for each frame in milliseconds

positions = [extra_positions, standard_positions]

def animate(i):
    # clear previous plots
    ax[0].clear() 
    ax[1].clear() 
    
    if (i % 10) == 0:
        print(f"Frame: {i}")
    
    # make into a for loop
    local_minima = [[], [], [], []]
    x_rel_lim = []
    y_rel_lim = []
    
    for k, position in enumerate(positions):
        nn = 3
        if num_bodies < 3:
            nn = 2
        for j in range(1, nn+1):
            x_label = f"x{j}"
            y_label = f"y{j}"
            # all x positions for sun/earth/moon
            x = position[x_label]
            y = position[y_label]

            # find initial minimum
            if j == 1:
                x_min = np.min(x)
                y_min = np.min(y)
                x_max = np.max(x)
                y_max = np.max(y)
            else:
                if np.min(x) < x_min:
                    local_x_min = np.min(x)
                    local_minima[0].append(local_x_min)
                if np.min(y) < y_min:
                    local_y_min = np.min(y)
                    local_minima[1].append(local_y_min)
                if np.max(x) > x_max:
                    local_x_max = np.max(x)
                    local_minima[2].append(local_x_max)
                if np.max(y) > y_max:
                    local_y_max = np.max(y)
                    local_minima[3].append(local_y_max)

            # x and y values up to current iteration
            x_current = x.iloc[:i + 1]
            y_current = y.iloc[:i + 1]

            # Plot all positions
            ax[0].plot(x_current, y_current, color=f"C{j+k}")
            # Marks starting position
            ax[0].scatter(x.iloc[0], y.iloc[0], marker='x', color=f"C{j+k}")
            # Marks current position
            ax[0].scatter(x_current.iloc[-1], y_current.iloc[-1], marker='o', color=f"C{j+k}")

        x_min = 1.1*np.min(local_minima[0]) #- np.abs(x_max) * 0.1
        y_min = 1.1*np.min(local_minima[1]) #- np.abs(y_max) * 0.1
        x_max = 1.1*np.max(local_minima[2]) #+ np.abs(x_max) * 0.1
        y_max = 1.1*np.max(local_minima[3]) # + np.abs(y_max) * 0.1

        ax[0].set_xlim(x_min, x_max)
        ax[0].set_ylim(y_min, y_max)

        # Loop through each body
        # Plot path and current position of each body
        
        for j in range(2, num_bodies + 1):
            colour = f"C{j+k}"
            if j==2:
                colour='C0'
                           
            x_label = f"x{j}"
            earth_x_label = "x2"
            y_label = f"y{j}"
            earth_y_label = "y2"
            # get all x posiitons for earth/moon
            x = position[x_label] 
            y = position[y_label] 
            # get all x positions relative to earth
            x_rel = position[x_label] - position[earth_x_label]
            y_rel = position[y_label] - position[earth_y_label]

            # x and y values up to current iteration
            x_current = x.iloc[:i + 1]
            y_current = y.iloc[:i + 1]
            
            # x and y values of positions relative to earth up to current iteration
            x_rel_current = x_rel.iloc[:i + 1]
            y_rel_current = y_rel.iloc[:i + 1]
            
            if j == 3:
                x_rel_lim.append( np.max([np.abs(np.min(x_rel)), np.max(x_rel)]))
                y_rel_lim.append( np.max([np.abs(np.min(y_rel)), np.max(y_rel)]))

            # Plot current position
            ax[1].scatter(x_rel_current, y_rel_current, marker='.', color=colour)
            # Plot the earth:
            ax[1].scatter(x_rel_current.iloc[-1], y_rel_current.iloc[-1], marker='.', color=colour)
        
    ax[1].set_xlim(-np.max(x_rel_lim)*1.1, np.max(x_rel_lim)*1.1)
    ax[1].set_ylim(-np.max(y_rel_lim)*1.1, np.max(y_rel_lim)*1.1)
        # ax[3].legend()

anim = animation.FuncAnimation(fig, animate, frames=num_frames, interval=interval)
gif_name = f"{extra_data_file[:-4]}_{args.frames}_frames"
anim.save(f"{gif_name}.gif")
print(f"gif saved as {gif_name}.gif")