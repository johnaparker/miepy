import h5py
import matplotlib.pyplot as plt

time_1 = {}
with h5py.File("out_1.h5", "r") as f:
    for name in f:
        time_1[name] = f[name].value

colors = {"build": "C0", "solve": "C1", "flux": "C2", "force": "C3"}


speedup_vals = [4, 8, 16]

for speedup in speedup_vals:
    fig, ax = plt.subplots()
    ax.axhline(speedup, color="black")

    time_x = {}
    with h5py.File(f"out_{speedup}.h5", "r") as f:
        for name in f:
            time_x[name] = f[name].value

        Nparticles = f.attrs["Nparticles"]

    for name in time_x.keys():
        ax.plot(Nparticles, time_1[name] / time_x[name], color=colors[name], label=name)

    ax.grid(axis="y")
    ax.set_ylim(ymin=0)
    ax.set(xlabel="number of particles", ylabel="speedup (relative to 1 cpu)")
    ax.set_title(f"Performance of parallized GMT using OpenMP ({speedup} cores)")
    ax.legend()

plt.show()
