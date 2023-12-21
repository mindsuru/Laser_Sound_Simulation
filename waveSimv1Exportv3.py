import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys
import cProfile

def de(String):
    print(String)
    sys.exit()

def waveFunctionAir(A, f, r, t):
    c = 343000                   # Schallgeschwindigkeit in mm/s
    lambda_ = c / f              # Wellenlänge
    k = 2 * np.pi / lambda_      # Wellenzahl
    omega = 2 * np.pi * f        # Kreisfrequenz
    epsilon = 1e-6               # Schwellenwert, um zu bestimmen, wann r als Null behandelt wird
    
    phi = -k * r
    # Vermeide Division durch sehr kleine Werte von r
    r_safe = np.where(np.abs(r) < epsilon, epsilon, r)

    # Berechne die Wellenfunktion
    waveFunction = A * np.cos(k * r_safe - omega * (t) + phi) / r_safe

    # Für sehr kleine r, verwende die alternative Berechnung
    waveFunction = np.where(np.abs(r) < epsilon, A * np.cos(-omega * (t)), waveFunction)

    
    return waveFunction

class MembraneSurface:
    def __init__(self, center, radius, frequency, amplitude, source):
        self.center = center
        self.radius = radius
        self.source = source
        self.amplitude = amplitude
        
        resolution = radius * 2 + 1
        x = np.linspace(-radius, radius, resolution)
        y = np.linspace(-radius, radius, resolution)
        X, Y = np.meshgrid(x, y)

        R = np.sqrt(X**2 + Y**2)
        self.mask = R <= radius - 1

        self.X, self.Y = X + center[0], Y + center[1]
        self.Z = np.zeros_like(self.X)

        self.velocity = np.zeros_like(self.X)
        self.F_excitation = np.zeros_like(self.X)
        self.F_elasticity = np.zeros_like(self.X)
        self.F_damping = np.zeros_like(self.X)
        self.frequency = frequency * 2 * np.pi
        self.dx = (2 * radius) / (resolution - 1)
        self.dy = (2 * radius) / (resolution - 1)

        self.immovable = np.zeros_like(self.X, dtype=bool)

        # Definiere einen erweiterten Randbereich
        extended_radius = radius - 1

        # Setze alle Punkte außerhalb des erweiterten Radius als immovable
        R = np.sqrt((self.X - center[0])**2 + (self.Y - center[1])**2)
        self.immovable[R > extended_radius] = True


    def update_Z(self, current_time, dt, elasticity, damping, c):
        self.Z[self.immovable] = 0
        self.velocity[self.immovable] = 0

        Z_xx = np.zeros_like(self.Z)
        Z_yy = np.zeros_like(self.Z)

        inner_area = ~self.immovable[1:-1, 1:-1]

        # Berechnung der zweiten räumlichen Ableitungen von Z
        Z_xx[1:-1, 1:-1][inner_area] = (np.roll(self.Z, -1, axis=1) - 2 * self.Z + np.roll(self.Z, 1, axis=1))[1:-1, 1:-1][inner_area] / self.dx**2
        Z_yy[1:-1, 1:-1][inner_area] = (np.roll(self.Z, -1, axis=0) - 2 * self.Z + np.roll(self.Z, 1, axis=0))[1:-1, 1:-1][inner_area] / self.dy**2
        
        # Berechnung der elastischen und Dämpfungskräfte
        self.F_elasticity = elasticity * (Z_xx + Z_yy)
        self.F_damping = damping * self.velocity

        # Integration der Schallquelle als Anregungsfunktion
        R = np.sqrt((self.X - self.source[0])**2 + (self.Y - self.source[1])**2 + (self.Z - self.source[2])**2)
        self.F_excitation[self.mask & ~self.immovable] = waveFunctionAir(self.amplitude, self.frequency / (2 * np.pi), R[self.mask & ~self.immovable], current_time)

        #print('elasticity: ', self.F_elasticity)
        #print('damping: ', self.F_damping)
        #print('excitation: ', self.F_excitation)

        contains_nan = np.isnan(self.F_excitation).any()
        if(contains_nan):
            sys.exit()
        F_totalPower = self.F_elasticity + self.F_excitation + self.F_damping

        
        # Hinzufügen des Gewichts in kg eines Punktes 
        mass = 0.000000001
        F_totalPower /= mass
        # Aktualisierung der Beschleunigung und der Geschwindigkeit
        a = c**2 * (Z_xx + Z_yy) + F_totalPower
        self.velocity += a * dt

        #print(self.velocity)

        # Aktualisierung der Auslenkung
        self.Z += self.velocity * dt

        #print(self.Z)     

    def get_XYZ(self):
        return self.X, self.Y, self.Z



def export_data_for_blender(frames, filename="animation_data.json"):
    with open(filename, "w") as file:
        json.dump(frames, file)
    print(f"Daten wurden erfolgreich nach {filename} exportiert.")

def animatePlay(frames, membrane, fps):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim([-membrane.radius, membrane.radius])
    ax.set_ylim([-membrane.radius, membrane.radius])
    ax.set_zlim([-membrane.radius, membrane.radius])
    surf = None

    def animate(i):
        nonlocal surf
        if surf:
            surf.remove()

        frame_data = frames[i]
        Z = np.zeros_like(membrane.X)
        for x, y, z in frame_data:
            # Umrechnen der Koordinaten in Indizes
            ix = np.searchsorted(membrane.X[0], x)
            iy = np.searchsorted(membrane.Y[:, 0], y)
            Z[iy, ix] = z

        surf = ax.plot_surface(membrane.X, membrane.Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
        return fig,

    anim = FuncAnimation(fig, animate, frames=len(frames), interval=1000/(fps), repeat=False)
    plt.show()


def selectFrame(simulation_steps, desired_frames):

    steps_per_frame = simulation_steps / desired_frames

    # Auswahl der Schritte für jeden Frame
    selected_steps = [int(round(steps_per_frame * i)) for i in range(desired_frames)]

    return selected_steps



def main():
    profiler = cProfile.Profile()
    frames = []                                     # Container für die Frames
    current_time = 0                                # Initalisierung des Startzeitpunktes

    ## Eigenschaften der Schallwelle
    frequency_hz = 880
    amplitude = 3
    #  Dauer einer Schwingungsperiode in Sekunden
    T_ = 1/frequency_hz

    ## Eigenschaften der Membran
    #  Elastizitätsmodul
    elasticity = 1e-6
    #  Dämpfung der Membranschwingung
    damping = 0.000005
    radius = 30
    #  Wellengeschwindigkeit in der Membran
    waveSpeed_membrane = 700000

    ## Eigenschaften der Animation
    #  Anzahl der abgebildeten Schwingungsperioden
    periods = 400
    #  Dauer der Animation in Sekunden
    animationDuration = 10
    fps = 24
    #  Anzahl der Frames der gesamten Animation
    num_frames = int(animationDuration * fps)

    ## Werte zur Berechnung der Animation
    #  Maximale Frequenz. Kann auch verändert werden, erhöht jedoch die benötigte Rechenzeit. Grundlage zur Auflösung der Simulation.
    maxFreq = 20000
    #  Zeit in Sekunden der kürzesten Schwingungsperiode, aktuell 0,00005 Sekunden, bzw 5*10^-5
    shortestT_ = 1/maxFreq
    #  100 Simulationsschritte pro 5*10^⁻5
    simStepPerShortestCycle = 100
    #  Dauer eines Simulationsschrittes aktuell 5*10^-7
    simulation_dt = shortestT_ / simStepPerShortestCycle
    simulation_steps = int((T_ * (periods / 10))/simulation_dt)

    selected_frames = selectFrame(simulation_steps, num_frames)
    
    source_position = (0, 0, 100)
    membran_center = (0, 0, 0)
    membrane = MembraneSurface(membran_center, radius, frequency_hz, amplitude, source_position)

    output_steps = max(1, int(simulation_steps / 20))
    
    for step in range(simulation_steps):
        if step % output_steps == 0:
            print(f"Simulation progress: {round(step/(simulation_steps/100))}/{simulation_steps/(simulation_steps/100)}% completed.")
   
        # Berechnung der externen Welle
        profiler.enable()
        membrane.update_Z(current_time, simulation_dt, elasticity, damping, waveSpeed_membrane)
        profiler.disable()
        current_time += simulation_dt
        
        if step in selected_frames:
            X, Y, Z = membrane.get_XYZ()
            frame_data = [(X[i][j], Y[i][j], Z[i][j]) for i in range(len(Z)) for j in range(len(Z[i])) if membrane.mask[i][j]]
            frames.append(frame_data)
    profiler.print_stats()
    export_data_for_blender(frames)
    animatePlay(frames, membrane, fps)

if __name__ == "__main__":
    main()
