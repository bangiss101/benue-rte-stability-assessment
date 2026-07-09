"""Compute Fourier-domain RTE transfer-function gain statistics.

The script estimates strict and robust condition numbers for a
Reduction-to-Equator operator using the supplied magnetic inclination,
declination, and grid cell size.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute RTE gain statistics and condition numbers.")
    parser.add_argument(
        "--input",
        default=r"C:\data\TMI_50m_minimum_curvature.csv",
        help="Input regular-grid XYZ CSV file with X,Y,Z columns.",
    )
    parser.add_argument(
        "--output",
        default=r"C:\data\RTE_condition_number_report.csv",
        help="Output CSV report file.",
    )
    parser.add_argument("--cell-size", type=float, default=50.0, help="Grid cell size in map units.")
    parser.add_argument("--inclination", type=float, default=-8.94, help="Observed magnetic inclination in degrees.")
    parser.add_argument("--declination", type=float, default=0.21, help="Observed magnetic declination in degrees.")
    parser.add_argument("--target-inclination", type=float, default=0.0, help="Target RTE inclination in degrees.")
    parser.add_argument("--target-declination", type=float, default=None, help="Target RTE declination in degrees. Default: same as observed declination.")
    parser.add_argument("--gain-floor-percentile", type=float, default=1.0, help="Percentile used for robust gain floor. Default: 1.")
    return parser.parse_args()


def read_grid_xyz(path: Path) -> tuple[pd.Index, pd.Index, np.ndarray]:
    df = pd.read_csv(path)
    df.columns = [c.strip().upper() for c in df.columns]
    required = {"X", "Y", "Z"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV must contain X, Y, Z columns. Found: {df.columns.tolist()}")

    df = df[["X", "Y", "Z"]].copy()
    for col in ["X", "Y", "Z"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["X", "Y"])
    df["X"] = df["X"].round(3)
    df["Y"] = df["Y"].round(3)

    x_unique = np.sort(df["X"].unique())
    y_unique = np.sort(df["Y"].unique())
    if len(x_unique) < 2 or len(y_unique) < 2:
        raise ValueError("Input must contain at least two unique X and Y coordinates.")

    grid = df.pivot_table(index="Y", columns="X", values="Z", aggfunc="mean")
    grid = grid.reindex(index=y_unique, columns=x_unique)
    return x_unique, y_unique, grid.to_numpy(dtype=float)


def direction_cosines(inclination_deg: float, declination_deg: float) -> tuple[float, float, float]:
    inclination = np.deg2rad(inclination_deg)
    declination = np.deg2rad(declination_deg)
    mx = np.cos(inclination) * np.sin(declination)
    my = np.cos(inclination) * np.cos(declination)
    mz = np.sin(inclination)
    return mx, my, mz


def main() -> None:
    args = parse_args()
    target_declination = args.declination if args.target_declination is None else args.target_declination

    x_unique, y_unique, _ = read_grid_xyz(Path(args.input))
    nx = len(x_unique)
    ny = len(y_unique)
    print(f"Grid dimensions: nx={nx}, ny={ny}")

    kx = 2 * np.pi * np.fft.fftfreq(nx, d=args.cell_size)
    ky = 2 * np.pi * np.fft.fftfreq(ny, d=args.cell_size)
    kx_grid, ky_grid = np.meshgrid(kx, ky)
    wavenumber = np.sqrt(kx_grid**2 + ky_grid**2)
    valid = wavenumber > 0

    mx, my, mz = direction_cosines(args.inclination, args.declination)
    mtx, mty, mtz = direction_cosines(args.target_inclination, target_declination)

    with np.errstate(divide="ignore", invalid="ignore"):
        observed = mz + 1j * (mx * kx_grid + my * ky_grid) / wavenumber
        target = mtz + 1j * (mtx * kx_grid + mty * ky_grid) / wavenumber
        transfer = target / observed

    transfer[~valid] = np.nan
    gain = np.abs(transfer)
    gain_valid = gain[np.isfinite(gain) & valid]
    if gain_valid.size == 0:
        raise ValueError("No finite gain values were computed.")

    gain_min = np.nanmin(gain_valid)
    gain_p1 = np.nanpercentile(gain_valid, args.gain_floor_percentile)
    gain_p5 = np.nanpercentile(gain_valid, 5.0)
    gain_p50 = np.nanpercentile(gain_valid, 50.0)
    gain_p95 = np.nanpercentile(gain_valid, 95.0)
    gain_p99 = np.nanpercentile(gain_valid, 99.0)
    gain_max = np.nanmax(gain_valid)

    kappa_strict = gain_max / gain_min if gain_min > 0 else np.inf
    kappa_robust_p1 = gain_max / gain_p1 if gain_p1 > 0 else np.inf
    kappa_robust_p5 = gain_max / gain_p5 if gain_p5 > 0 else np.inf

    report = pd.DataFrame(
        {
            "Parameter": [
                "Observed inclination I",
                "Observed declination D",
                "Target inclination It",
                "Target declination Dt",
                "Cell size",
                "Grid nx",
                "Grid ny",
                "Gain minimum",
                "Gain P1",
                "Gain P5",
                "Gain median P50",
                "Gain P95",
                "Gain P99",
                "Gain maximum",
                "Strict condition number max/min",
                "Robust condition number max/P1",
                "Robust condition number max/P5",
            ],
            "Value": [
                args.inclination,
                args.declination,
                args.target_inclination,
                target_declination,
                args.cell_size,
                nx,
                ny,
                gain_min,
                gain_p1,
                gain_p5,
                gain_p50,
                gain_p95,
                gain_p99,
                gain_max,
                kappa_strict,
                kappa_robust_p1,
                kappa_robust_p5,
            ],
        }
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(output_path, index=False)

    print("\nRTE condition-number assessment completed.")
    print("----------------------------------------")
    print(f"Observed inclination: {args.inclination}°")
    print(f"Observed declination: {args.declination}°")
    print(f"Target inclination: {args.target_inclination}°")
    print(f"Target declination: {target_declination}°")
    print(f"Gain max: {gain_max:.6f}")
    print(f"Gain P1: {gain_p1:.6f}")
    print(f"Gain P5: {gain_p5:.6f}")
    print(f"Gain median P50: {gain_p50:.6f}")
    print(f"Robust condition number max/P1: {kappa_robust_p1:.6f}")
    print(f"Robust condition number max/P5: {kappa_robust_p5:.6f}")
    print(f"\nSaved report: {output_path}")


if __name__ == "__main__":
    main()
