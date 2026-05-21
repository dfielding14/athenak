#!/usr/bin/env python3
"""Validate and optionally generate fixed-format AthenaK cooling tables."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import isfinite, pow
from pathlib import Path
import sys


VALID_UNITS = {"code", "cgs"}
VALID_AXIS_KINDS = {"temperature", "density", "scalar", "metallicity"}
VALID_SCALES = {"linear", "log10"}
VALID_BOUNDS = {"zero", "clamp", "fatal"}
VALID_VALUE_KINDS = {"lambda", "gamma", "modifier"}
DATA_ORDER = "c_row_major_last_axis_fastest"
MAX_AXES = 3


class TableError(RuntimeError):
  pass


@dataclass
class Axis:
  index: int
  kind: str
  units: str
  scale: str
  xmin: float
  xmax: float
  n: int
  options: dict[str, str]


@dataclass
class CoolingTable:
  path: Path
  ndim: int
  value_kind: str
  value_units: str
  value_scale: str
  bounds: str
  axes: list[Axis]
  values: list[float]


def parse_float(token: str, context: str) -> float:
  try:
    value = float(token)
  except ValueError as exc:
    raise TableError(f"{context}: could not parse real value {token!r}") from exc
  if not isfinite(value):
    raise TableError(f"{context}: value must be finite")
  return value


def parse_int(token: str, context: str) -> int:
  try:
    value = int(token)
  except ValueError as exc:
    raise TableError(f"{context}: could not parse integer value {token!r}") from exc
  return value


def noncomment_tokens(path: Path) -> list[tuple[int, list[str]]]:
  lines: list[tuple[int, list[str]]] = []
  for lineno, raw in enumerate(path.read_text().splitlines(), start=1):
    body = raw.split("#", 1)[0].strip()
    if body:
      lines.append((lineno, body.split()))
  return lines


def parse_axis(path: Path, lineno: int, key: str, fields: list[str]) -> Axis:
  try:
    index = int(key[4:])
  except ValueError as exc:
    raise TableError(f"{path}:{lineno}: malformed axis key {key!r}") from exc
  if index < 0 or index >= MAX_AXES:
    raise TableError(f"{path}:{lineno}: axis index must be 0, 1, or 2")
  if len(fields) < 5:
    raise TableError(
        f"{path}:{lineno}: expected axisN kind units [linear|log10] xmin xmax n"
    )
  kind = fields[0].lower()
  units = fields[1].lower()
  if kind not in VALID_AXIS_KINDS:
    raise TableError(f"{path}:{lineno}: bad axis kind {fields[0]!r}")
  if units not in VALID_UNITS:
    raise TableError(f"{path}:{lineno}: bad axis units {fields[1]!r}")

  offset = 2
  scale = "log10"
  if fields[offset].lower() in VALID_SCALES:
    scale = fields[offset].lower()
    offset += 1
  else:
    try:
      float(fields[offset])
    except ValueError as exc:
      raise TableError(
          f"{path}:{lineno}: optional axis scale must be linear or log10; "
          f"got {fields[offset]!r}"
      ) from exc
  if len(fields) < offset + 3:
    raise TableError(f"{path}:{lineno}: incomplete axis extent")
  xmin = parse_float(fields[offset], f"{path}:{lineno} axis xmin")
  xmax = parse_float(fields[offset + 1], f"{path}:{lineno} axis xmax")
  n = parse_int(fields[offset + 2], f"{path}:{lineno} axis n")
  if xmax <= xmin:
    raise TableError(f"{path}:{lineno}: axis xmax must exceed xmin")
  if n < 2:
    raise TableError(f"{path}:{lineno}: axis n must be at least 2")
  if scale == "log10":
    try:
      raw_min = pow(10.0, xmin)
      raw_max = pow(10.0, xmax)
    except OverflowError as exc:
      raise TableError(
          f"{path}:{lineno}: log10 axis extent overflows raw coordinates"
      ) from exc
    if raw_min <= 0.0 or raw_max <= 0.0 or not isfinite(raw_min + raw_max):
      raise TableError(
          f"{path}:{lineno}: log10 axis extent must map to finite positive "
          "raw coordinates"
      )
  options: dict[str, str] = {}
  for token in fields[offset + 3:]:
    if "=" not in token:
      if kind == "density":
        options["density_kind"] = token
        continue
      raise TableError(f"{path}:{lineno}: malformed axis option {token!r}")
    opt_key, opt_value = token.split("=", 1)
    if opt_key not in {"density_kind", "scalar_index"}:
      raise TableError(f"{path}:{lineno}: unknown axis option {opt_key!r}")
    options[opt_key] = opt_value
  if kind == "density" and options.get("density_kind", "mass_density") not in {
      "mass_density", "number_density", "hydrogen_number_density"}:
    raise TableError(f"{path}:{lineno}: bad density_kind option")
  if kind in {"scalar", "metallicity"} and "scalar_index" in options:
    if parse_int(options["scalar_index"], f"{path}:{lineno} scalar_index") < 0:
      raise TableError(f"{path}:{lineno}: scalar_index must be non-negative")
  return Axis(index, kind, units, scale, xmin, xmax, n, options)


def parse_table(path: Path, expect_value_kind: str | None = None) -> CoolingTable:
  rows = noncomment_tokens(path)
  if not rows:
    raise TableError(f"{path}: empty table")
  first_line, first = rows[0]
  if first != ["ATHENAK_COOLING_TABLE", "1"]:
    raise TableError(
        f"{path}:{first_line}: first non-comment line must be "
        "ATHENAK_COOLING_TABLE 1"
    )

  ndim: int | None = None
  value_kind: str | None = None
  value_units = "code"
  value_scale = "linear"
  bounds = "zero"
  axes: dict[int, Axis] = {}
  saw_data = False
  value_tokens: list[tuple[int, str]] = []

  for lineno, fields in rows[1:]:
    if saw_data:
      value_tokens.extend((lineno, token) for token in fields)
      continue
    key = fields[0]
    values = fields[1:]
    if key == "data":
      if values:
        raise TableError(f"{path}:{lineno}: data line must not contain values")
      saw_data = True
    elif key == "ndim":
      if len(values) != 1:
        raise TableError(f"{path}:{lineno}: ndim expects one value")
      ndim = parse_int(values[0], f"{path}:{lineno} ndim")
      if ndim < 1 or ndim > MAX_AXES:
        raise TableError(f"{path}:{lineno}: ndim must be 1, 2, or 3")
    elif key == "value_kind":
      if len(values) != 1:
        raise TableError(f"{path}:{lineno}: value_kind expects one value")
      value_kind = values[0].lower()
      if value_kind not in VALID_VALUE_KINDS:
        raise TableError(f"{path}:{lineno}: bad value_kind {values[0]!r}")
    elif key == "value_units":
      if len(values) != 1 or values[0].lower() not in VALID_UNITS:
        raise TableError(f"{path}:{lineno}: value_units must be code or cgs")
      value_units = values[0].lower()
    elif key == "value_scale":
      if len(values) != 1 or values[0].lower() not in VALID_SCALES:
        raise TableError(f"{path}:{lineno}: value_scale must be linear or log10")
      value_scale = values[0].lower()
    elif key == "bounds":
      if len(values) != 1 or values[0].lower() not in VALID_BOUNDS:
        raise TableError(f"{path}:{lineno}: bounds must be zero, clamp, or fatal")
      bounds = values[0].lower()
    elif key == "data_order":
      if len(values) != 1 or values[0] != DATA_ORDER:
        raise TableError(f"{path}:{lineno}: data_order must be {DATA_ORDER}")
    elif key.startswith("axis"):
      axis = parse_axis(path, lineno, key, values)
      if axis.index in axes:
        raise TableError(f"{path}:{lineno}: duplicate axis{axis.index}")
      axes[axis.index] = axis
    else:
      raise TableError(f"{path}:{lineno}: unknown header token {key!r}")

  if ndim is None:
    raise TableError(f"{path}: missing ndim")
  if value_kind is None:
    raise TableError(f"{path}: missing value_kind")
  if expect_value_kind is not None and value_kind != expect_value_kind.lower():
    raise TableError(
        f"{path}: value_kind is {value_kind}, expected {expect_value_kind.lower()}"
    )
  if not saw_data:
    raise TableError(f"{path}: missing data marker")
  for axis_index in range(ndim):
    if axis_index not in axes:
      raise TableError(f"{path}: missing axis{axis_index}")
  for axis_index in axes:
    if axis_index >= ndim:
      raise TableError(f"{path}: axis{axis_index} is defined beyond ndim={ndim}")

  ordered_axes = [axes[index] for index in range(ndim)]
  expected_values = 1
  for axis in ordered_axes:
    expected_values *= axis.n
  values = [parse_float(token, f"{path}:{lineno} data") for lineno, token in value_tokens]
  if len(values) != expected_values:
    raise TableError(
        f"{path}: has {len(values)} values but expected {expected_values}"
    )
  return CoolingTable(path, ndim, value_kind, value_units, value_scale, bounds,
                      ordered_axes, values)


def write_example(path: Path, ndim: int, value_kind: str) -> None:
  if ndim < 1 or ndim > MAX_AXES:
    raise TableError("--ndim must be 1, 2, or 3")
  value_kind_token = {
      "lambda": "Lambda",
      "gamma": "Gamma",
      "modifier": "modifier",
  }[value_kind.lower()]
  axes = [
      "axis0 temperature code linear 0.0 2.0 3",
      "axis1 density code linear 0.0 2.0 3 density_kind=mass_density",
      "axis2 scalar code linear 0.0 1.0 3 scalar_index=0",
  ]
  values = []
  for i in range(3):
    for j in range(3 if ndim > 1 else 1):
      for k in range(3 if ndim > 2 else 1):
        if value_kind.lower() == "modifier":
          values.append(1.0 + 0.1*i + 0.01*j + 0.001*k)
        elif value_kind.lower() == "gamma":
          values.append(1.0e-3*(1.0 + i + j + k))
        else:
          values.append(1.0e-2*(1.0 + i + j + k))
  lines = [
      "ATHENAK_COOLING_TABLE 1",
      f"ndim {ndim}",
      f"value_kind {value_kind_token}",
      "value_units code",
      "value_scale linear",
      "bounds zero",
      f"data_order {DATA_ORDER}",
      *axes[:ndim],
      "data",
      *(f"{value:.16e}" for value in values),
  ]
  path.parent.mkdir(parents=True, exist_ok=True)
  path.write_text("\n".join(lines) + "\n")


def plot_table(table: CoolingTable, output: Path) -> None:
  import matplotlib
  matplotlib.use("Agg")
  import matplotlib.pyplot as plt
  import numpy as np

  shape = tuple(axis.n for axis in table.axes)
  data = np.array(table.values).reshape(shape)
  if table.value_scale == "log10":
    data = np.power(10.0, data)
  coords = [np.linspace(axis.xmin, axis.xmax, axis.n) for axis in table.axes]
  fig, ax = plt.subplots(figsize=(6.2, 4.3), constrained_layout=True)
  if table.ndim == 1:
    ax.plot(coords[0], data, marker="o")
    ax.set_xlabel(f"{table.axes[0].kind} ({table.axes[0].scale})")
    ax.set_ylabel(table.value_kind)
  else:
    image = data if table.ndim == 2 else data[:, :, table.axes[2].n//2]
    mesh = ax.pcolormesh(coords[1], coords[0], image, shading="auto")
    ax.set_xlabel(f"{table.axes[1].kind} ({table.axes[1].scale})")
    ax.set_ylabel(f"{table.axes[0].kind} ({table.axes[0].scale})")
    if table.ndim == 3:
      ax.set_title(f"central {table.axes[2].kind} slice")
    fig.colorbar(mesh, ax=ax, label=table.value_kind)
  output.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(output, dpi=160)
  plt.close(fig)


def main(argv: list[str] | None = None) -> int:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("table", nargs="?", type=Path, help="table to validate")
  parser.add_argument("--expect-value-kind", choices=sorted(VALID_VALUE_KINDS))
  parser.add_argument("--plot", type=Path, help="write a quick-look plot")
  parser.add_argument("--write-example", type=Path, help="write an example table")
  parser.add_argument("--ndim", type=int, default=1, help="example table dimension")
  parser.add_argument("--value-kind", default="lambda",
                      choices=sorted(VALID_VALUE_KINDS),
                      help="example table value kind")
  args = parser.parse_args(argv)

  try:
    if args.write_example is not None:
      write_example(args.write_example, args.ndim, args.value_kind)
      table_path = args.write_example
    else:
      if args.table is None:
        parser.error("table is required unless --write-example is used")
      table_path = args.table
    table = parse_table(table_path, args.expect_value_kind)
    if args.plot is not None:
      plot_table(table, args.plot)
  except TableError as exc:
    print(f"ERROR: {exc}", file=sys.stderr)
    return 1
  print(
      f"OK: {table.path} ndim={table.ndim} value_kind={table.value_kind} "
      f"values={len(table.values)}"
  )
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
