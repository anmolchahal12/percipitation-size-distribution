"""Plot measured CMSX-4 gamma-prime equivalent diameters (Figures 7 and 8)."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import norm

CONDITIONS = (15, 60, 165)
SHEETS = {15: "15_min_0p25_h", 60: "60_min_1_h", 165: "165_min_2p75_h"}
COLORS = ("#007F86", "#CF7531", "#7563A6")


def load_data(path: Path) -> pd.DataFrame:
    """Read tidy CSV or the supplied workbook's particle sheets, never its summary."""
    if path.suffix.lower() == ".csv":
        frame = pd.read_csv(path)
    elif path.suffix.lower() == ".xlsx":
        parts = []
        for minutes, sheet in SHEETS.items():
            raw = pd.read_excel(path, sheet_name=sheet, header=1)
            if "D_eq (nm)" not in raw:
                raise ValueError(f"{sheet}: missing D_eq (nm) column")
            raw = raw.dropna(how="all")
            parts.append(pd.DataFrame({"aging_time_min": minutes,
                                       "diameter_nm": raw["D_eq (nm)"]}))
        frame = pd.concat(parts, ignore_index=True)
    else:
        raise ValueError("Input must be a .csv or .xlsx file")
    required = ["aging_time_min", "diameter_nm"]
    if not set(required).issubset(frame):
        raise ValueError("CSV requires aging_time_min and diameter_nm columns")
    frame = frame[required].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(frame.to_numpy()).all() or (frame.diameter_nm <= 0).any():
        raise ValueError("Every row must have a finite, positive diameter and finite aging time")
    if set(frame.aging_time_min) != set(CONDITIONS):
        raise ValueError("Provide all three aging conditions: 15, 60, 165 minutes")
    for minutes, group in frame.groupby("aging_time_min"):
        if len(group) < 2 or group.diameter_nm.std(ddof=1) == 0:
            raise ValueError(f"{minutes:g} min needs at least two distinct diameters")
    return frame


def summarize(frame: pd.DataFrame) -> pd.DataFrame:
    """Use sample standard deviation (ddof=1) for descriptive Gaussian curves."""
    rows = []
    for minutes in CONDITIONS:
        d = frame.loc[frame.aging_time_min == minutes, "diameter_nm"].to_numpy()
        rows.append(dict(aging_time_min=minutes, n=len(d), mean_nm=np.mean(d),
                         sample_std_nm=np.std(d, ddof=1), median_nm=np.median(d)))
    return pd.DataFrame(rows)


def empirical_cdf(d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values, counts = np.unique(d, return_counts=True)
    return values, np.cumsum(counts) / len(d)


def plot_figures(frame: pd.DataFrame, output: Path, bins: int = 15,
                 label: str = "") -> pd.DataFrame:
    output.mkdir(parents=True, exist_ok=True)
    stats = summarize(frame)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "savefig.facecolor": "white"})
    xmax = max(frame.diameter_nm.max() * 1.08,
               (stats.mean_nm + 4 * stats.sample_std_nm).max())
    x = np.linspace(0, xmax, 1000)
    edges = np.linspace(0, xmax, bins + 1)
    suffix = f" · {label}" if label else ""

    fig7, axes = plt.subplots(1, 3, figsize=(15, 5.3), sharex=True, sharey=True,
                             layout="constrained")
    fig7.suptitle("Figure 7 | γ′ precipitate size distributions" + suffix, fontsize=15)
    fig8, (pdf_ax, cdf_ax) = plt.subplots(1, 2, figsize=(13, 5.5), layout="constrained")
    fig8.suptitle("Figure 8 | Comparison across aging conditions" + suffix, fontsize=15)
    for ax, color, row in zip(axes, COLORS, stats.itertuples()):
        d = frame.loc[frame.aging_time_min == row.aging_time_min, "diameter_nm"].to_numpy()
        density = norm.pdf(x, loc=row.mean_nm, scale=row.sample_std_nm)
        ax.hist(d, bins=edges, density=True, color=color, alpha=0.28,
                edgecolor="white", label="Measured diameters")
        ax.plot(x, density, color=color, lw=2.4, label="Gaussian approximation")
        ax.axvline(row.mean_nm, color="#283747", ls="--", lw=1.6,
                   label=f"Mean: {row.mean_nm:.1f} nm")
        ax.axvline(row.median_nm, color=color, ls=":", lw=2,
                   label=f"Median: {row.median_nm:.1f} nm")
        ax.set_title(f"{row.aging_time_min} min · N = {row.n}\n"
                     f"Sample SD = {row.sample_std_nm:.1f} nm", fontsize=12)
        ax.legend(fontsize=8, loc="upper right")
        pdf_ax.plot(x, density, color=color, lw=2.4,
                    label=f"{row.aging_time_min} min · μ = {row.mean_nm:.1f}, s = {row.sample_std_nm:.1f} nm")
        values, cumulative = empirical_cdf(d)
        cdf_ax.step(np.r_[0, values, xmax], np.r_[0, cumulative, 1], where="post",
                    color=color, lw=2, label=f"{row.aging_time_min} min (N = {row.n})")
        cdf_ax.plot(x, norm.cdf(x, row.mean_nm, row.sample_std_nm),
                    color=color, ls="--", alpha=0.7)
        cdf_ax.plot(row.median_nm, .5, "D", color=color, ms=5)
    axes[0].set_ylabel("Probability density (nm⁻¹)")
    pdf_ax.set(title="(a) Gaussian probability density", ylabel="Probability density (nm⁻¹)")
    cdf_ax.set(title="(b) Cumulative size distribution", ylabel="Cumulative fraction", ylim=(0, 1.04))
    pdf_ax.legend(fontsize=9)
    handles, labels = cdf_ax.get_legend_handles_labels()
    handles += [Line2D([], [], color="#45515B", ls="--"),
                Line2D([], [], color="#45515B", marker="D", ls="")]
    cdf_ax.legend(handles, labels + ["Gaussian CDF", "Median D₅₀"], fontsize=9, loc="lower right")
    for ax in [*axes, pdf_ax, cdf_ax]:
        ax.set_xlabel("Equivalent circle diameter (nm)")
        ax.set_xlim(0, xmax)
        ax.grid(axis="y", alpha=0.15)
        ax.set_axisbelow(True)
    for fig, stem in [(fig7, "figure_7_probability_density"), (fig8, "figure_8_gaussian_cdf")]:
        for extension in ("png", "svg"):
            fig.savefig(output / f"{stem}.{extension}", dpi=220)
        plt.close(fig)
    stats.to_csv(output / "summary_statistics.csv", index=False)
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Tidy CSV or supplied-format XLSX")
    parser.add_argument("--output", type=Path, default=Path("figures"))
    parser.add_argument("--bins", type=int, default=15)
    parser.add_argument("--label", default="", help="Optional dataset label printed on both figures")
    args = parser.parse_args()
    if args.bins < 1:
        parser.error("--bins must be a positive integer")
    try:
        frame = load_data(args.input)
        print(plot_figures(frame, args.output, args.bins, args.label).to_string(index=False))
    except (ValueError, OSError) as exc:
        parser.exit(2, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
