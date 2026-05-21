#!/usr/bin/env python3
"""Generate the streaming-instability eigenmode used by the particle tests."""

import argparse

import numpy as np


def nsh_state(epsilon, tau_s, eta_vk):
    """Return the gas and particle NSH drift velocities used by the pgen."""
    denom = (1.0 + epsilon) ** 2 + tau_s ** 2
    gas = np.array(
        [
            2.0 * epsilon * tau_s * eta_vk / denom,
            -(1.0 + epsilon + tau_s ** 2) * eta_vk / denom,
            0.0,
        ]
    )
    dust = np.array(
        [
            -2.0 * tau_s * eta_vk / denom,
            -(1.0 + epsilon) * eta_vk / denom,
            0.0,
        ]
    )
    return gas, dust


def streaming_matrix(args):
    """Build gamma q = A q for the compressible, axisymmetric dust-gas system."""
    rho_g = args.rho0
    rho_d = args.epsilon * args.rho0
    omega = args.omega0
    t_stop = args.tau_s / omega
    kx = 2.0 * np.pi * args.nwx / args.lx
    kz = 2.0 * np.pi * args.nwz / args.lz
    gas, dust = nsh_state(args.epsilon, args.tau_s, args.eta_vk)
    delta_v = dust - gas

    mat = np.zeros((8, 8), dtype=complex)
    irhog, iugx, iugy, iugz, irhod, ivdx, ivdy, ivdz = range(8)

    mat[irhog, irhog] += -1j * kx * gas[0]
    mat[irhog, iugx] += -1j * rho_g * kx
    mat[irhog, iugz] += -1j * rho_g * kz

    mat[irhod, irhod] += -1j * kx * dust[0]
    mat[irhod, ivdx] += -1j * rho_d * kx
    mat[irhod, ivdz] += -1j * rho_d * kz

    mat[iugx, iugx] += -1j * kx * gas[0] - args.epsilon / t_stop
    mat[iugx, iugy] += 2.0 * omega
    mat[iugx, ivdx] += args.epsilon / t_stop
    mat[iugx, irhog] += -1j * kx * args.cs ** 2 / rho_g
    mat[iugx, irhog] += -args.epsilon * delta_v[0] / (t_stop * rho_g)
    mat[iugx, irhod] += args.epsilon * delta_v[0] / (t_stop * rho_d)

    mat[iugy, iugx] += -(2.0 - args.qshear) * omega
    mat[iugy, iugy] += -1j * kx * gas[0] - args.epsilon / t_stop
    mat[iugy, ivdy] += args.epsilon / t_stop
    mat[iugy, irhog] += -args.epsilon * delta_v[1] / (t_stop * rho_g)
    mat[iugy, irhod] += args.epsilon * delta_v[1] / (t_stop * rho_d)

    mat[iugz, iugz] += -1j * kx * gas[0] - args.epsilon / t_stop
    mat[iugz, ivdz] += args.epsilon / t_stop
    mat[iugz, irhog] += -1j * kz * args.cs ** 2 / rho_g

    mat[ivdx, iugx] += 1.0 / t_stop
    mat[ivdx, ivdx] += -1j * kx * dust[0] - 1.0 / t_stop
    mat[ivdx, ivdy] += 2.0 * omega

    mat[ivdy, iugy] += 1.0 / t_stop
    mat[ivdy, ivdx] += -(2.0 - args.qshear) * omega
    mat[ivdy, ivdy] += -1j * kx * dust[0] - 1.0 / t_stop

    mat[ivdz, iugz] += 1.0 / t_stop
    mat[ivdz, ivdz] += -1j * kx * dust[0] - 1.0 / t_stop
    return mat


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho0", type=float, default=1.0)
    parser.add_argument("--epsilon", type=float, default=1.0)
    parser.add_argument("--tau-s", dest="tau_s", type=float, default=0.3)
    parser.add_argument("--eta-vk", dest="eta_vk", type=float, default=0.5)
    parser.add_argument("--omega0", type=float, default=1.0)
    parser.add_argument("--qshear", type=float, default=1.5)
    parser.add_argument("--cs", type=float, default=1.0)
    parser.add_argument("--lx", type=float, default=1.0)
    parser.add_argument("--lz", type=float, default=1.0)
    parser.add_argument("--nwx", type=int, default=2)
    parser.add_argument("--nwz", type=int, default=4)
    return parser.parse_args()


def main():
    args = parse_args()
    eigvals, eigvecs = np.linalg.eig(streaming_matrix(args))
    index = int(np.argmax(eigvals.real))
    eigval = eigvals[index]
    eigvec = eigvecs[:, index] / eigvecs[4, index]

    print(f"problem/expected_growth = {eigval.real:.17g}")
    print(f"problem/expected_frequency = {eigval.imag:.17g}")
    names = [
        "rhog",
        "vgx",
        "vgy",
        "vgz",
        "rhod",
        "vdx",
        "vdy",
        "vdz",
    ]
    for name, value in zip(names, eigvec):
        print(f"problem/eigen_{name}_re = {value.real:.17e}")
        print(f"problem/eigen_{name}_im = {value.imag:.17e}")


if __name__ == "__main__":
    main()
