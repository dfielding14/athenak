#!/usr/bin/env python3
## Alireza Rashti - Oct 2024 (C)
## usage:
## $ ./me -h
##

import sys
import os
import numpy as np
import math
import argparse
import re
import h5py

# import matplotlib.pyplot as plt
# import glob
# import sympy
## ---------------------------------------------------------------------- ##

## field names
g_field_names = [
    "gxx",
    "gxy",
    "gxz",
    "gyy",
    "gyz",
    "gzz",
    "betax",
    "betay",
    "betaz",
    "alp",
]


def parse_cli():
  """
    arg parser
    """
  p = argparse.ArgumentParser(
      description="convert Athenak CCE dumps to Specter CCE")
  p.add_argument("-f_h5", type=str, required=True, help="/path/to/cce/h5/dumps")
  p.add_argument("-d_out", type=str, required=True, help="/path/to/output/dir")
  p.add_argument(
      "-t_deriv",
      type=str,
      default="Fourier",
      help="method to take the time derivative of fields:{Fourier}",
  )

  args = p.parse_args()
  return args


def load(fpath: str, field_name: str, attr: dict) -> list:
  """
    read the field accroding to attr.
    return convention:
      ret[i] = [dump_time_value,
                2d_array_coeffs_for_given_time_real,
                2d_array_coeffs_for_given_time_imag],

      where i indicates the dump number.
    """

  ret = []
  if attr["file_type"] == "h5":
    lev_t = attr["lev_t"]
    max_n = attr["max_n"]
    max_l = attr["max_l"]
    with h5py.File(fpath, "r") as h5f:
      # get shape and dtype
      key = f"{0}"
      h5_re = h5f[f"{key}/{field_name}/re"]
      h5_im = h5f[f"{key}/{field_name}/im"]
      re = np.empty_like(h5_re)
      im = np.empty_like(h5_im)

      # read & save
      for i in range(0, lev_t):
        key = f"{i}"
        t = h5f[key].attrs["Time"][0]
        h5_re = h5f[f"{key}/{field_name}/re"]
        h5_im = h5f[f"{key}/{field_name}/im"]
        re[:] = h5f[f"{key}/{field_name}/re"]
        im[:] = h5f[f"{key}/{field_name}/im"]
        # save
        ret.append([t, re, im])
  else:
    raise ValueError("no such option")

  return ret


def lm_mode(l, m):
  """
    l and m mode convention
    """
  return l * l + l + m


def get_attribute(fpath: str,
                  field_name: str = "gxx",
                  type: str = "h5",
                  args=None) -> dict:
  """
    find attributes such as num. of time level, and n,l,m in C_nlm
    """
  attr = {}
  if type == "h5":
    attr["file_type"] = "h5"
    with h5py.File(fpath, "r") as h5f:
      # find attribute about num. of time level, and n,l,m in C_nlm
      attr["lev_t"] = len(h5f.keys()) - 1
      attr["max_n"] = h5f["1/gxx/re"].shape[0]
      attr["max_l"] = int(math.sqrt(h5f["1/gxx/re"].shape[1]))

  else:
    raise ValueError("no such option")

  return attr


def time_derivative_fourier(field_name: str, db: dict, args):
  """
    return the time derivative of the given field using Fourier method
    """

  # populate
  # note field[i] = [time, re, im]
  field = db[field_name]["value"]
  len_all_modes = field[1].shape[0] * field[1].shape[1]
  len_all_time = len(field[0])
  coeff = np.empty(shape=(len_all_time))
  f_coeff = np.empty(shape=(len_all_modes, len_all_time))

  ## it's 2d and not 1d ???
  for nlm in range(0, len_all_modes):
    # fetch data along time axis for each nlm coeffs
    for i in range(0, len_all_time):
      coeff[i] = complex(field[i][1][nlm], field[i][2][nlm])

    # F. transform
    f_coeff[nlm] = np.fft.fft(coeff[i])

  # F. derivative
  wm = math.pi * 2.0 / (nmax * dt)
  f_coeff[0] = 0
  for nlm in range(0, len_all_modes):
    for i in range(1, len_all_time // 2 + 1):
      omega = i * wm

      re = np.real(f_coeff[i])
      im = np.imag(f_coeff[i])

      re2 = np.real(f_coeff[len_all_time - i])
      im2 = np.imag(f_coeff[len_all_time - i])

      f_coeff[i] = np.complex(-im * omega, re * omega)
      f_coeff[len_all_time - i] = np.complex(-im2 * omega, re2 * omega)

  # F. inverse
  for nlm in range(0, len_all_modes):
    d_coeff[nlm] = np.fft.ifft(f_coeff[nlm])


def time_derivative(field_name: str, db: dict, args):
  """
    return the time derivative of the given field
    """

  if args["t_deriv"] == "Fourier":
    return time_derivative_fourier(field_name, db, args)
  else:
    raise ValueError("no such option")


def main(args):
  """
    create output required by Specter code
    ref: https://spectre-code.org/tutorial_cce.html
    """

  # find attribute for an arbitrary field
  attr = get_attribute(args["f_h5"], "gxx", "h5")

  # for each field
  field_name = "gxx"

  # load data
  db[field_name] = {}
  db[field_name]["value"] = load(args["f_h5"], field_name, attr)

  # print(dat)

  # time derivative
  db[field_name]["Dt"] = time_derivative(field_name, db, args)

  # radial derivative

  # save


if __name__ == "__main__":
  args = parse_cli()
  main(args.__dict__)
