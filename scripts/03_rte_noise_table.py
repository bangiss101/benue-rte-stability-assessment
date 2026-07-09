"""Compare high-frequency variability before and after RTE.

The high-frequency component is estimated as the residual obtained after
subtracting a moving-window local mean from each grid. The diagnostic reports
TMI and RTE high-frequency standard deviations and their ratio.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute TMI/RTE high-frequency noise-assessment table.")
    parser.add_argument(
        "--tmi",
        default=r"C:\data\TMI_50m_minimum_curvature.csv",
        help="Input TMI XYZ CSV file with X,Y,Z columns.",
    )
    parser.add_argument(
        "--rte",
        default=r"C:\data\RTE_50m.csv",
        help="Input RTE XYZ CSV file with X,Y,Z columns.",
    )
    parser.add_argument(
        "--output",
        default=r"C:\data\RTE_noise_assessment_table.csv",
        help="Output CSV report file.",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=15,
        help="Moving-window size in grid cells. Default: 15.",
    )
    return parser.parse_args()


def read_xyz_grid(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip().upper() for c in df.columns]
    required = {"X", "Y", "Z"}
    if not required.issubset(df.columns):
        raise ValueError(f"{path} must contain X, Y, Z columns. Found: {df.columns.tolist()}")

    df = df[["X", "Y", "Z"]].copy()
    for col in ["X", "Y", "Z"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=["X", "Y"])
    df["X"] = df["X"].round(3)
    df["Y"] = df["Y"].round(3)

    grid = df.pivot_table(index="Y", columns="X", values="Z", aggfunc="mean")
    if grid.empty:
        raise ValueError(f"No valid grid values were found in {path}.")
    return grid


def basic_stats(array: np.ndarray) -> dict[str, float]:
    return {
        "Minimum": np.nanmin(array),
        "Maximum": np.nanmax(array),
        "Mean": np.nanmean(array),
        "Standard deviation": np.nanstd(array),
        "P1": np.nanpercentile(array, 1),
        "P99": np.nanpercentile(array, 99),
        "Robust amplitude range P99-P1": np.nanpercentile(array, 99) - np.nanpercentile(array, 1),
    }


def high_frequency_std(array: np.ndarray, size: int) -> float:
    if size < 3:
        raise ValueError("window-size must be at least 3 cells.")

    valid = np.isfinite(array)
    array_zero = np.where(valid, array, 0.0)
    weight = valid.astype(float)

    local_sum = uniform_filter(array_zero, size=size, mode="nearest")
    local_weight = uniform_filter(weight, size=size, mode="nearest")
    local_mean = local_sum / np.maximum(local_weight, 1e-12)
    residual = array - local_mean
    return float(np.nanstd(residual))


def main() -> None:
    args = parse_args()

    print("Reading TMI grid...")
    tmi_grid = read_xyz_grid(Path(args.tmi))
    print("Reading RTE grid...")
    rte_grid = read_xyz_grid(Path(args.rte))

    print("Original TMI shape:", tmi_grid.shape)
    print("Original RTE shape:", rte_grid.shape)

    common_y = tmi_grid.index.intersection(rte_grid.index)
    common_x = tmi_grid.columns.intersection(rte_grid.columns)
    if len(common_x) == 0 or len(common_y) == 0:
        raise ValueError("No common X-Y overlap between TMI and RTE grids.")

    tmi_common = tmi_grid.loc[common_y, common_x].to_numpy(dtype=float)
    rte_common = rte_grid.loc[common_y, common_x].to_numpy(dtype=float)
    print("Common overlap shape:", tmi_common.shape)

    valid_overlap = np.isfinite(tmi_common) & np.isfinite(rte_common)
    n_total = tmi_common.size
    n_valid = int(np.sum(valid_overlap))
    valid_percentage = (n_valid / n_total) * 100

    print("Total common-overlap cells:", n_total)
    print("Valid overlapping cells:", n_valid)
    print("Valid overlap percentage:", valid_percentage)
    if n_valid == 0:
        raise ValueError("No overlapping valid cells between TMI and RTE grids.")

    tmi_valid = np.where(valid_overlap, tmi_common, np.nan)
    rte_valid = np.where(valid_overlap, rte_common, np.nan)

    tmi_stats = basic_stats(tmi_valid)
    rte_stats = basic_stats(rte_valid)
    tmi_hf_std = high_frequency_std(tmi_valid, size=args.window_size)
    rte_hf_std = high_frequency_std(rte_valid, size=args.window_size)
    noise_ratio = rte_hf_std / tmi_hf_std

    rows = [
        {
            "Grid": "TMI before RTE",
            **tmi_stats,
            "High-frequency standard deviation": tmi_hf_std,
            "Valid cells used": n_valid,
            "Valid overlap percentage": valid_percentage,
        },
        {
            "Grid": "RTE after transformation",
            **rte_stats,
            "High-frequency standard deviation": rte_hf_std,
            "Valid cells used": n_valid,
            "Valid overlap percentage": valid_percentage,
        },
        {
            "Grid": "RTE/TMI high-frequency noise ratio",
            "Minimum": np.nan,
            "Maximum": np.nan,
            "Mean": np.nan,
            "Standard deviation": np.nan,
            "P1": np.nan,
            "P99": np.nan,
            "Robust amplitude range P99-P1": np.nan,
            "High-frequency standard deviation": noise_ratio,
            "Valid cells used": n_valid,
            "Valid overlap percentage": valid_percentage,
        },
    ]

    table = pd.DataFrame(rows)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, index=False)

    print("\nRTE noise assessment completed.")
    print("--------------------------------")
    print(f"TMI high-frequency std: {tmi_hf_std:.6f}")
    print(f"RTE high-frequency std: {rte_hf_std:.6f}")
    print(f"RTE/TMI noise ratio: {noise_ratio:.6f}")
    print(f"Valid overlapping cells: {n_valid:,}")
    print(f"Saved table: {output_path}")


if __name__ == "__main__":
    main()
