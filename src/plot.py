import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import tikzplotlib

sns.set_theme()


def tweak_tikz_plots(filename):
    """Tweaks the tikz plots to make them look better

    Parameters
    ----------
        filename : str
            The filename of the tikz plot to be tweaked
    """
    with open(filename, "r") as f:
        lines = f.readlines()

    with open(filename, "w") as f:
        for line in lines:
            if "majorticks" in line:
                f.write(line.replace("false", "true"))
            elif "addplot" in line:
                f.write(line.replace("semithick", "thick"))
            else:
                f.write(line)


def save_tikz(filename):
    """Saves the plot as a tikz-tex file

    Parameters
    ----------
        filename : str
            The filename of the tikz plot to be saved
    """
    plt.grid(True)
    tikzplotlib.clean_figure()
    tikzplotlib.save(filename)
    tweak_tikz_plots(filename)
    plt.clf()


def plot_burn_in_times(l_values, plot_lengths, temperatures=[1, 2.4]):
    """Plots the burn in times for different temperatures and lattice sizes

    Parameters
    ----------
        l_values : list
            The lattice sizes to be used
        plot_lengths : list
            The lengths of the burn in times to be plotted
        temperatures : list
            The temperatures to be used
    """
    assert (
        len(l_values) == len(plot_lengths) / 2
    ), "l_values and plot_lenghts must have the same length"

    for i, L in enumerate(l_values):
        for j, value_type in enumerate(["Energy", "Magnetization"]):
            for randomness, order_type in zip(
                ["random", "nonrandom"], ["Unordered", "Ordered"]
            ):
                for T in temperatures:
                    filename = f"burn_in_L_{L}_T_{T:.6f}_{randomness}.csv"

                    df = pd.read_csv("output/" + filename)

                    if value_type == "Energy":
                        plt.plot(
                            df.N[: plot_lengths[2 * i + j]],
                            df.expected_E[: plot_lengths[2 * i + j]],
                            label=f"{order_type}, $T={T}$",
                        )
                    else:
                        plt.plot(
                            df.N[: plot_lengths[2 * i + j]],
                            df.expected_M[: plot_lengths[2 * i + j]],
                            label=f"{order_type}, $T={T}$",
                        )

            plt.legend()
            estimated_what = (
                r"$\langle\epsilon\rangle$"
                if value_type == "Energy"
                else r"$\langle|m|\rangle$"
            )
            plt.title(rf"Estimated {estimated_what} for L={L}")
            if value_type == "Energy":
                plt.ylabel(rf"{estimated_what} $[J/k_B]$")
            else:
                plt.ylabel(rf"{estimated_what} $[1]$")
            plt.xlabel(r"$N$")

            filename = f"plots/burn_in/burn_in_L_{L}_{value_type}.tex"
            save_tikz(filename)


def plot_probability_distribution():
    """Plots probability distribution for different temperatures"""
    L = "20"
    for T in ["1.0", "2.1", "2.4"]:
        df = pd.read_csv(f"output/samples_L={L}_T={T}.csv")
        plt.title(fr"Estimated probability distribution of $\epsilon$ at $T={T}J/k_B$")
        plt.xlabel(r"$\epsilon$ $[J]$")
        plt.ylabel(fr"$p(\epsilon; {T} J / k_B)$")
        plt.hist(df.epsilon, bins="auto", density=True)
        filename = f"plots/distributions/epsilon_L={L}_T={T}.tex"
        save_tikz(filename)
        print(f"Variance at T={T}: {df.epsilon.var()}")
        print(f"Expected value at T={T}: {df.epsilon.mean()}")


def plot_values():
    """Plots expected values for different lattice sizes"""
    all_L = range(40, 160, 20)

    dfs = {L: pd.read_csv(f"output/values_L={L}.csv") for L in all_L}
    for i, (value, ylabel, unit) in enumerate(
        zip(
            ["<epsilon>", "<|m|>", "C_v", "chi"],
            [r"\langle\epsilon\rangle", r"\langle|m|\rangle", "C_v", r"\chi"],
            ["J/k_B", "1", "k_B", "1 / J"],
        )
    ):
        for L in all_L:
            df = dfs[L]
            df.sort_values("T", inplace=True, ignore_index=True)
            plt.plot(df["T"], df[value], label=f"L={L}")

        plt.title(f"Estimated values of ${ylabel}$ for different T")
        plt.ylabel(f"${ylabel}$ [${unit}$]")
        plt.xlabel("$T$ $[J / k_B]$")
        plt.legend()

        filename = f"plots/values/values_{value}.tex"
        save_tikz(filename)


def estimate_T_inf(value: str):
    """Estimates T_inf based observed values

    Parameters
    ----------
        value : str
            The value to base the estimate of T_inf
    """
    y = []
    x = []
    label = value if value == "C_v" else r"\chi"
    for L in range(40, 160, 20):
        df = pd.read_csv(f"output/values_zoom_L={L}.csv")
        argmax = df[value].idxmax()
        y.append(df.iloc[argmax]["T"])
        x.append(1 / L)
    linear_fit = stats.linregress(x, y)
    plt.title(
        r"Observations of $T_c(L)$ against $L^{-1}$ and linear fit to find $T_c(\infty)$"
    )
    plt.scatter(x, y, label=r"Observed $T_c$")
    estimate = linear_fit.intercept
    plt.scatter(
        [0] * len(x),
        [estimate] * len(x),
        s=40,
        label=fr"$T_c(\infty) = {estimate: .5f}$",
    )
    plt.plot([0, max(x)], estimate + linear_fit.slope * np.asarray([0, max(x)]))
    plt.legend()
    plt.title(fr"Estimating $T_c(\infty)$ using ${label}$")
    plt.ylabel("$T_c$ [$J / k_B$]")
    plt.xlabel("$L^{-1}$ [1]")

    filename = f"plots/T_inf/T_inf_{value}.tex"
    save_tikz(filename)


def main():
    """Main function"""
    plot_burn_in_times(l_values=[20, 100], plot_lengths=[6_000, 20_000, 1_500, 3_000])
    plot_probability_distribution()
    plot_values()
    estimate_T_inf("C_v")
    estimate_T_inf("chi")


if __name__ == "__main__":
    main()
