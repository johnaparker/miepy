"""Brownian dynamics of N=91 Ag nanospheres in a Gaussian beam optical trap.

Hexagonal lattice (L=5 layers, 91 particles) of 75nm Ag spheres
in water, illuminated by a wide Gaussian beam. Optical forces from
miepy drive stoked Brownian dynamics at T=300K.

Usage:
    python examples/hex_optical_trapping.py                  # default CPU backend
    python examples/hex_optical_trapping.py --backend gpu    # GPU solver (requires miepy[gpu])
    python examples/hex_optical_trapping.py --steps 5000 --backend gpu
    python examples/hex_optical_trapping.py --animate              # show 2D animation
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import stoked
import miepy

nm = 1e-9
us = 1e-6


def hexagonal_lattice_layers(L, spacing):
    """Return positions [N, 3] for a hexagonal lattice with L layers."""
    k1 = np.array([1, 0, 0], dtype=float)
    k2 = np.array([np.cos(np.pi / 3), np.sin(np.pi / 3), 0], dtype=float)

    lattice = [np.zeros(3)]
    for step in range(1, L + 1):
        lattice.append(step * k1)
        for direc in [k2 - k1, -k1, -k2, -k2 + k1, k1, k2]:
            for _ in range(step):
                lattice.append(lattice[-1] + direc)
        lattice.pop()

    return np.asarray(lattice, dtype=float) * spacing


def main():
    parser = argparse.ArgumentParser(description="Optical trapping of hexagonal Ag nanoparticle cluster")
    parser.add_argument(
        "--backend", choices=["cpu", "gpu"], default="cpu", help="Solver backend: cpu (default) or gpu (requires miepy[gpu])"
    )
    parser.add_argument("--steps", type=int, default=10000, help="Number of BD steps")
    parser.add_argument("--dt", type=float, default=5.0, help="Time step in microseconds")
    parser.add_argument("--no-plot", action="store_true", help="Skip plotting")
    parser.add_argument("--animate", action="store_true", help="Show 2D animation of particle trajectories")
    args = parser.parse_args()

    # --- Physical parameters ---
    radius = 75 * nm
    spacing = 700 * nm
    wavelength = 800 * nm
    L = 1  # hexagonal layers -> N = 3*L*(L+1) + 1 = 91

    # Gaussian beam: waist much wider than cluster so all particles are illuminated
    cluster_radius = L * spacing
    w0 = 3 * cluster_radius  # beam waist ~3x the cluster extent
    n_medium = 1.33
    k = 2 * np.pi * n_medium / wavelength
    zR = k * w0**2 / 2  # Rayleigh range

    source = miepy.sources.gaussian_beam(
        polarization=[1, 1j],  # circular polarization
        width=w0,
        power=0.03,  # 30 mW
        center=[0, 0, zR],
    )

    water = miepy.materials.water()
    Ag = miepy.materials.Ag()

    # --- Initial positions ---
    initial_positions = hexagonal_lattice_layers(L, spacing)
    N = len(initial_positions)

    # Per-particle radii: center particle is larger
    radii = np.full(N, radius)
    radii[0] = 100 * nm

    print(f"Backend: {args.backend}")
    print(f"Particles: N={N}, L={L}, spacing={spacing / nm:.0f} nm")
    print(f"Beam waist: w0={w0 / nm:.0f} nm, Rayleigh range: zR={zR / nm:.0f} nm")

    # --- Build miepy cluster ---
    cluster = miepy.sphere_cluster(
        position=initial_positions,
        radius=radii,
        material=Ag,
        source=source,
        wavelength=wavelength,
        lmax=2,
        medium=water,
    )

    # --- Set solver backend ---
    miepy.backends.set_backend(args.backend)

    # --- Force callback for stoked ---
    def optical_force(time, position, orientation):
        """Compute optical forces by updating miepy cluster positions."""
        pos_3d = np.column_stack([position, np.zeros(len(position))])
        cluster.update_position(pos_3d)
        cluster.solve()
        F = cluster.force()  # [N, 3], always C++
        return F[:, :2]  # [N, 2] for 2D simulation

    # --- Stoked Brownian dynamics (2D) ---
    dt = args.dt * us
    viscosity = 8.9e-4  # water at ~25 C (Pa*s)

    # Screened electrostatic repulsion (DLVO double-layer)
    electrostatic = stoked.double_layer_sphere(
        radius=radii,
        potential=-40e-3,  # -40 mV surface potential (typical for Ag in water)
    )

    bd = stoked.brownian_dynamics(
        temperature=300,
        dt=dt,
        position=initial_positions[:, :2],  # 2D (x, y)
        drag=stoked.drag_sphere(radius=radii, viscosity=viscosity),
        force=optical_force,
        interactions=electrostatic,
    )

    # --- Warmup JIT if using JAX ---
    if args.backend == "gpu":
        print("Warming up JAX JIT...")
        cluster.solve()
        print("JIT warmup done")

    # --- Run simulation ---
    Nsteps = args.steps
    print(f"Running {Nsteps} steps, dt={dt / us:.1f} us, total time={Nsteps * dt * 1e3:.1f} ms")

    result = bd.run(Nsteps, progress=True)
    trajectory = result.position  # [Nsteps, N, 2]
    print(f"Trajectory shape: {trajectory.shape}")

    # --- Save ---
    outfile = f"hex_trapping_{args.backend}.npy"
    np.save(outfile, trajectory)
    print(f"Saved trajectory to {outfile}")

    if args.no_plot:
        return

    # --- Visualize ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Per-particle colors: center red, rest blue
    pcolors = ["red"] + ["blue"] * (N - 1)

    # Left: particle trajectories
    ax = axes[0]
    for i in range(N):
        ax.plot(trajectory[:, i, 0] / nm, trajectory[:, i, 1] / nm, lw=0.3, alpha=0.5, color=pcolors[i])
    ax.plot(initial_positions[:, 0] / nm, initial_positions[:, 1] / nm, "ko", ms=3, zorder=5, label="initial")
    ax.set(aspect="equal", xlabel="x (nm)", ylabel="y (nm)", title="Particle trajectories")
    ax.legend(fontsize=8)

    # Right: first and last frame
    ax = axes[1]
    for i in range(N):
        ax.plot(trajectory[0, i, 0] / nm, trajectory[0, i, 1] / nm, "o", ms=4, color=pcolors[i], alpha=0.5)
        ax.plot(trajectory[-1, i, 0] / nm, trajectory[-1, i, 1] / nm, "s", ms=4, color=pcolors[i], alpha=0.5)
    ax.plot([], [], "o", color="gray", ms=4, label="t=0")
    ax.plot([], [], "s", color="gray", ms=4, label=f"t={Nsteps * dt * 1e3:.0f} ms")
    ax.set(aspect="equal", xlabel="x (nm)", ylabel="y (nm)", title="Initial vs final positions")
    ax.legend(fontsize=8)

    plt.tight_layout()
    plt.savefig(f"hex_trapping_{args.backend}.png", dpi=150)
    plt.show()

    # --- Animation ---
    if args.animate:
        from stoked.vis import trajectory_animation, circle_patches

        traj_nm = stoked.trajectory(trajectory[::100] / nm)
        patches = circle_patches(radii / nm)
        time = np.arange(0, Nsteps, 100) * dt * 1e3  # milliseconds

        fig_anim, ax_anim = plt.subplots(figsize=(8, 8))
        ax_anim.set(xlabel="x (nm)", ylabel="y (nm)", title="Optical trapping dynamics")
        colors = ["red"] + ["blue"] * (N - 1)
        anim = trajectory_animation(
            traj_nm, patches=patches, ax=ax_anim, colors=colors, time=time, time_unit="ms", trail=50
        )
        plt.show()


if __name__ == "__main__":
    main()
