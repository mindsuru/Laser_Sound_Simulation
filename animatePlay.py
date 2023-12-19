import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import json


def import_data(filename):
    # Daten aus der JSON-Datei lesen
    with open(filename, 'r') as file:
        data = json.load(file)
    return data

def create_grid(radius):
    x = np.linspace(-radius, radius, 50)
    y = np.linspace(-radius, radius, 50)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    return X, Y, Z

def animatePlay(frames, radius, fps):
    X, Y, Z = create_grid(radius)  # Koordinatenraster erstellen
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim([-radius, radius])
    ax.set_ylim([-radius, radius])
    ax.set_zlim([-radius/10000000, radius/10000000])  # Z-Achse angepasst
    surf = None

    def animate(i):
        nonlocal surf
        if surf:
            surf.remove()

        frame_data = frames[i]
        for x, y, z in frame_data:
            # Umrechnen der Koordinaten in Indizes
            ix = np.searchsorted(X[0], x)
            iy = np.searchsorted(Y[:, 0], y)
            Z[iy, ix] = z

        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
        return fig,

    anim = FuncAnimation(fig, animate, frames=len(frames), interval=1000/(fps), repeat=True)
    plt.show()

filename = "/home/alex/Dokumente/VsCode/Skripte/Python_Skript/Laser_Sound_Simulation/animation_data.json"
frames = import_data(filename)
animatePlay(frames, 40, 200)