#!/usr/bin/env python
from __future__ import division
from io import StringIO

import numpy as np
import matplotlib.pyplot as plt


def plot_verlet_times():
    """Plot times generated by run_verlet_times."""
    data = np.genfromtxt("verlet_times.txt")
    rv, real, user, sys = data.T
    plt.title("Verlet list radius vs. total CPU time")
    plt.xlabel(r"Verlet Radius ($\sigma$)")
    plt.ylabel("CPU time (s)")
    plt.plot(rv, user + sys)
    plt.grid()
    plt.savefig("verlet_times", dpi=600)


def plot_efficiency_multi():
    """Plot efficiency data generated by run_deltat_efficiency like the right
    plot in figure 15.3, page 429 of the book.
    """
    tmax = 5
    with open("n_efficiency.txt") as f:
        data = f.read()
    timesteps = data.strip().split("\n\n")
    fig = plt.figure(figsize=(8, 8))
    for delt_data in timesteps:
        delt, d = delt_data.split('\n', 1)
        delt = float(delt.split()[-1])
        n, real, user, sys = np.genfromtxt(StringIO(unicode(d))).T
        simulation_time = np.array([tmax] * len(user))
        cpu_time = user + sys
        efficiency = simulation_time / cpu_time
        plt.semilogy(n, efficiency, label="%.1e" % delt)
    plt.title("MTS efficiency vs number of small steps")
    plt.xlabel("n")
    plt.ylabel(r"Efficiency ($\eta$)")
    plt.legend(ncol=3)
    plt.grid()
    plt.savefig("efficiency_multi", dpi=600)


def plot_efficiency():
    """Plot efficiency data for leapfrog and MTS (like figure 15.3 left)."""
    tmax = 5
    fig = plt.figure(figsize=(10, 8))

    # Leapfrog
    delt, real, user, sys = np.genfromtxt("deltat_efficiency_leapfrog.txt").T
    cpu_time = user + sys
    simulation_time = np.array([tmax] * len(user))
    efficiency = simulation_time / cpu_time
    plt.loglog(delt, efficiency, label="Leapfrog")

    # Multiple Time Step
    with open("n_efficiency.txt") as f:
        data = f.read()
    timesteps = data.strip().split("\n\n")
    dt = []
    mts = np.zeros([10, len(timesteps)])
    for (i, delt_data) in enumerate(timesteps):
        delt, d = delt_data.split('\n', 1)
        delt = float(delt.split()[-1])
        dt.append(delt)
        n, real, user, sys = np.genfromtxt(StringIO(unicode(d))).T
        simulation_time = np.array([tmax] * len(user))
        cpu_time = user + sys
        efficiency = simulation_time / cpu_time
        for (n, e) in enumerate(efficiency):
            mts[n, i] = e

    for (n, s) in enumerate(mts):
        plt.loglog(dt, s, label="MTS n = %d" % (n + 1))
    plt.title("Efficiency vs Time Step")
    plt.xlabel(r"$\Delta t$")
    plt.ylabel(r"$\eta$")
    plt.legend(ncol=3, loc="upper left")
    plt.grid()
    plt.savefig("efficiency", dpi=600)


def plot_energy_drift():
    """Plot results generated by run_energy_drift (like figure 15.2)."""
    data = (("energy_drift_leapfrog.txt", "Leapfrog", 1),
            ("energy_drift_multi_n3.txt", "MTS n = 3", 3),
            ("energy_drift_multi_n4.txt", "MTS n = 4", 4),
            ("energy_drift_multi_n10.txt", "MTS n = 10", 10))
    for (s, label, n) in data:
        dt, energy_drft = np.genfromtxt(s).T
        dt = dt * n  # big time step for MTS
        plt.loglog(dt, energy_drft, label=label)
    plt.title("Energy Drift vs Time step")
    plt.xlabel(r"$\Delta t$")
    plt.ylabel("E")
    plt.grid()
    plt.legend(loc="lower right", ncol=2)
    plt.savefig("energy_drift", dpi=400)
