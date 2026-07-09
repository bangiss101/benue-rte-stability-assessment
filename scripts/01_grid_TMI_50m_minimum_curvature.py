"""Regrid TMI XYZ data to a regular grid using GMT/pyGMT.

The workflow applies blockmean followed by surface, which is a
continuous-curvature/minimum-curvature-type interpolation commonly used for
potential-field gridding.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pygmt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Regrid TMI XYZ data to a regular grid using pyGMT blockmean + surface."
    )
    parser.add_argument(
        "--input",
        default=r"C:\data\TMI.csv",
        help="Input XYZ CSV file with X,Y,Z columns. Default: C:\\data\\TMI.csv",
    )
    parser.add_argument(
        "--out-csv",
        default=r"C:\data\TMI_50m_minimum_curvature.csv",
        help="Output regular-grid XYZ CSV file.",
    )
    parser.add_argument(
        "--out-nc",
        default=r"C:\data\TMI_50m_minimum_curvature.nc",
        help="Output NetCDF grid file.",
    )
    parser.add_argument(
        "--cell-size",
        type=float,
        default=50.0,
        help="Output grid cell size in map units, normally metres. Default: 50.",
    )
    parser.add_argument(
        "--tension",
        type=float,
        default=0.25,
        help="pyGMT surface tension factor. Default: 0.25.",
    )
    return parser.parse_args()


def read_xyz(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip().upper() for c in df.columns]
    required = {"X", "Y", "Z"}
    if not required.issubset(df.columns):
        raise ValueError(f"CSV must contain X, Y, Z columns. Found: {df.columns.tolist()}")

    df = df[["X", "Y", "Z"]].copy()
    for col in ["X", "Y", "Z"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["X", "Y", "Z"])
    if df.empty:
        raise ValueError("No valid numeric X,Y,Z rows were found in the input file.")
    return df


def main() -> None:
    args = parse_args()
    input_csv = Path(args.input)
    output_grid_nc = Path(args.out_nc)
    output_xyz_csv = Path(args.out_csv)
    cell_size = float(args.cell_size)

    output_grid_nc.parent.mkdir(parents=True, exist_ok=True)
    output_xyz_csv.parent.mkdir(parents=True, exist_ok=True)

    df = read_xyz(input_csv)

    xmin = np.floor(df["X"].min() / cell_size) * cell_size
    xmax = np.ceil(df["X"].max() / cell_size) * cell_size
    ymin = np.floor(df["Y"].min() / cell_size) * cell_size
    ymax = np.ceil(df["Y"].max() / cell_size) * cell_size

    region = [xmin, xmax, ymin, ymax]
    spacing = str(cell_size)

    print("Input rows:", len(df))
    print("Region:", region)
    print("Spacing:", spacing, "map units")
    print("Running blockmean...")
    df_block = pygmt.blockmean(data=df, region=region, spacing=spacing)

    print("Blockmean rows:", len(df_block))
    print("Running surface gridding...")
    grid = pygmt.surface(data=df_block, region=region, spacing=spacing, tension=args.tension)

    grid.to_netcdf(output_grid_nc)

    print("Exporting grid to XYZ CSV...")
    xyz = grid.to_dataframe(name="Z").reset_index()
    xyz = xyz.rename(columns={"x": "X", "y": "Y"})
    xyz = xyz[["X", "Y", "Z"]]
    xyz.to_csv(output_xyz_csv, index=False)

    print("Done.")
    print(f"Output grid: {output_grid_nc}")
    print(f"Output XYZ CSV: {output_xyz_csv}")
    print(f"Output cell size: {cell_size}")
    print(f"Output rows: {len(xyz):,}")


if __name__ == "__main__":
    main()
