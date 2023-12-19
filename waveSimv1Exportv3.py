import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

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


    def update_Z(self, current_time, dt, elasticity, damping):
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

        #print(self.F_elasticity)
        #print(self.F_damping)
        #print(self.F_excitation)

        F_totalPower = self.F_elasticity + self.F_excitation + self.F_damping
        #print(F_totalPower)
        # Aktualisierung der Beschleunigung und der Geschwindigkeit
        c = 900000                                          # in mm/s
        a = c**2 * (Z_xx + Z_yy) + F_totalPower
        self.velocity += a * dt

        #print(self.velocity)

        # Aktualisierung der Auslenkung
        self.Z += self.velocity * dt

        #print(self.Z)

    def update_ZBK(self, current_time, dt, elasticity, damping):
        # Setze zuerst die Auslenkung der Randpunkte auf 0
        self.Z[self.immovable] = 0

        # Berechne die Auslenkung für die restlichen Punkte
        R = np.sqrt((self.X - self.source[0])**2 + (self.Y - self.source[1])**2 + (self.Z - self.source[2])**2)
        self.Z[self.mask & ~self.immovable] = waveFunctionAir(self.amplitude, self.frequency / (2 * np.pi), R[self.mask & ~self.immovable], current_time)

        # Erstelle Slices für die zentrale Region und die Nachbarn
        central = self.Z[1:-1, 1:-1]
        
        top = self.Z[:-2, 1:-1]
        bottom = self.Z[2:, 1:-1]
        left = self.Z[1:-1, :-2]
        right = self.Z[1:-1, 2:]
        
        
        delta_Z = np.zeros_like(central)  # Initialisiere delta_Z mit Nullen
        delta_Z[~self.immovable[1:-1, 1:-1]] = (top + bottom + left + right)[~self.immovable[1:-1, 1:-1]] - 4 * central[~self.immovable[1:-1, 1:-1]]

        # Wende das Hooke'sche Gesetz und die Dämpfung an
        force = elasticity * delta_Z - damping * central

        # Aktualisiere Z, ignoriere immovable Punkte
        self.Z[1:-1, 1:-1][~self.immovable[1:-1, 1:-1]] += force[~self.immovable[1:-1, 1:-1]]

        def plot_masked_points(self):
            # Berechne den tatsächlichen Mittelpunkt der Maske
            mask_center = (self.center[0] - self.radius, self.center[1] - self.radius)

            masked_points = np.argwhere(self.mask)
            masked_points = [(point[0] + mask_center[0], point[1] + mask_center[1]) for point in masked_points]

            x_values_mask = [point[0] for point in masked_points]
            y_values_mask = [point[1] for point in masked_points]

            immov_points = np.argwhere(self.immovable)
            immov_points = [(point[0] + mask_center[0], point[1] + mask_center[1]) for point in immov_points]

            x_values_immov = [point[0] for point in immov_points]
            y_values_immov = [point[1] for point in immov_points]

            # Überprüfe, ob an einer Koordinate sowohl maskierte als auch unbewegliche Punkte liegen
            overlapping_points = set(masked_points) & set(immov_points)

            x_values_overlap = [point[0] for point in overlapping_points]
            y_values_overlap = [point[1] for point in overlapping_points]

            plt.scatter(x_values_mask, y_values_mask, color='blue', label='Maskierte Punkte')
            plt.scatter(x_values_immov, y_values_immov, color='red', label='Unbewegliche Punkte')
            plt.scatter(x_values_overlap, y_values_overlap, color='green', label='Überlappende Punkte')

            plt.xlabel('X-Achse')
            plt.ylabel('Y-Achse')
            plt.title('Punkte im Koordinatensystem')
            plt.grid(True)
            plt.show()

            sys.exit()
        

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
    ax.set_zlim([-membrane.radius/1000000000, membrane.radius/1000000000])
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

    anim = FuncAnimation(fig, animate, frames=len(frames), interval=1000/(fps), repeat=True)
    plt.show()


def selectFrame(simulation_steps, desired_frames):

    steps_per_frame = simulation_steps / desired_frames

    # Auswahl der Schritte für jeden Frame
    selected_steps = [int(round(steps_per_frame * i)) for i in range(desired_frames)]

    return selected_steps



def main():
    frames = []                                     # Container für die Frames
    current_time = 0                                # Initalisierung des Startzeitpunktes
    frequency_hz = 12000                            # Frequenz der Schwingung
    
    amplitude = 300
    T_ = 1/frequency_hz                             # Dauer der Echtzeitabbildung in ms
    elasticity = 1e-6
    damping = 2e-6
    radius = 30                                     # Radius der Membran

    periodDuration = 0.0016                              # Dauer einer Schwingungsperiode
    periods = 12000                                     # Anzahl der gezeigten Schwingungen
    animationDuration = periodDuration * periods    # Dauer der gesamten Animation in Sekunden
    
    maxFreq = 20000
    shortestT_ = 1/maxFreq
    simStepPerCycle = 1000

    simulation_steps = simStepPerCycle * periods


    fps = 24                                        # FPS der Animation
    num_frames = int(animationDuration * fps)            # Anzahl der Frames der gesamten Animation
    simFrameRatio = 100                           # Anzahl der Simulationsschritte pro Frame

    frame_duration = T_ / (periodDuration * fps)    # Realzeit zwischen Frames
    simulation_dt = shortestT_ / simStepPerCycle  # Realzeit zwischen Simulationsschritt

    selected_frames = selectFrame(simulation_steps, num_frames)
    
    source_position = (0, 0, 20)
    membran_center = (0, 0, 0)
    membrane = MembraneSurface(membran_center, radius, frequency_hz, amplitude, source_position)

    output_steps = max(1, int(simulation_steps / 20))
    
    print(simulation_dt)
    for step in range(simulation_steps):
        if step % output_steps == 0:
            print(f"Simulation progress: {step}/{simulation_steps} steps completed.")
   
        # Berechnung der externen Welle
        # distance = np.sqrt((membrane.X - source_position[0])**2 + (membrane.Y - source_position[1])**2 + (membrane.Z - source_position[2])**2)
        membrane.update_Z(current_time, simulation_dt, elasticity, damping)
        current_time += simulation_dt
        
        if step in selected_frames:
            X, Y, Z = membrane.get_XYZ()
            frame_data = [(X[i][j], Y[i][j], Z[i][j]) for i in range(len(Z)) for j in range(len(Z[i])) if membrane.mask[i][j]]
            frames.append(frame_data)

    export_data_for_blender(frames)
    animatePlay(frames, membrane, fps)

if __name__ == "__main__":
    main()
