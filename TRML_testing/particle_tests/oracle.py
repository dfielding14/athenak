#!/usr/bin/env python3
"""Analytic oracle helpers for the TRML particle-tracking tests."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class MeshSpec:
    nx1: int
    nx2: int
    nx3: int
    x1min: float = 0.0
    x1max: float = 1.0
    x2min: float = 0.0
    x2max: float = 1.0
    x3min: float = 0.0
    x3max: float = 1.0

    @property
    def dx1(self) -> float:
        return (self.x1max - self.x1min) / self.nx1

    @property
    def dx2(self) -> float:
        return (self.x2max - self.x2min) / self.nx2

    @property
    def dx3(self) -> float:
        return (self.x3max - self.x3min) / self.nx3

    def center(self, i: int, j: int, k: int) -> tuple[float, float, float]:
        return (
            self.x1min + (i + 0.5) * self.dx1,
            self.x2min + (j + 0.5) * self.dx2,
            self.x3min + (k + 0.5) * self.dx3,
        )

    def cell_indices(self, x: float, y: float, z: float) -> tuple[int, int, int]:
        i = min(self.nx1 - 1, max(0, int(math.floor((x - self.x1min) / self.dx1))))
        j = min(self.nx2 - 1, max(0, int(math.floor((y - self.x2min) / self.dx2))))
        k = min(self.nx3 - 1, max(0, int(math.floor((z - self.x3min) / self.dx3))))
        return i, j, k


def splitmix_unit(seed: int) -> float:
    z = seed & 0xFFFFFFFFFFFFFFFF
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
    z = (z ^ (z >> 31)) & 0xFFFFFFFFFFFFFFFF
    return (z & 0x7FFFFFFF) / float(0x80000000)


def initial_cell_center_particles(mesh: MeshSpec) -> list[tuple[float, float, float]]:
    return [
        mesh.center(i, j, k)
        for k in range(mesh.nx3)
        for j in range(mesh.nx2)
        for i in range(mesh.nx1)
    ]


def outer_x1_particles(mesh: MeshSpec) -> list[tuple[float, float, float]]:
    i = mesh.nx1 - 1
    return [
        mesh.center(i, j, k)
        for k in range(mesh.nx3)
        for j in range(mesh.nx2)
    ]


def analytic_state(
    mesh: MeshSpec, x: float, y: float, z: float, level: int = 0, gamma: float = 5.0 / 3.0
) -> dict[str, float]:
    i, j, k = mesh.cell_indices(x, y, z)
    rho = 10.0 + i + 10.0 * j + 100.0 * k + 1000.0 * level
    temp = 1.0 + 0.01 * i + 0.1 * j + k
    press = rho * temp
    eint = press / (gamma - 1.0)
    scalar0 = x + 2.0 * y + 3.0 * z
    return {
        "rho": rho,
        "press": press,
        "temp": temp,
        "eint": eint,
        "scalar0": scalar0,
    }


def uniform_state(
    rho: float = 1.0, temp: float = 2.0, scalar0: float = 0.25,
    gamma: float = 5.0 / 3.0
) -> dict[str, float]:
    press = rho * temp
    return {
        "rho": rho,
        "press": press,
        "temp": temp,
        "eint": press / (gamma - 1.0),
        "scalar0": scalar0,
    }


def predicted_x1_mc_push(
    tag: int, x: float, mesh: MeshSpec, probability: float, random_seed: int,
    cycle_used_by_pusher: int = 0, periodic: bool = True
) -> float:
    draw = splitmix_unit(tag * 7919 + cycle_used_by_pusher * 104729 + random_seed)
    if draw >= probability:
        return x
    x_new = x + mesh.dx1
    if periodic and x_new >= mesh.x1max:
        x_new -= mesh.x1max - mesh.x1min
    return x_new
