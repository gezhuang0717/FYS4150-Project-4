from pandas.io.parsers import read_csv
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import subprocess

sns.set_theme()


def plot_abs_m_unordered(L=20):
    """Plots expected value for epsilon and abs(m) for T=1.0
    and unordered initial spins

    Parameters
    ----------
        L : int
            Lattice size
    """
    filename = "burn_in_L_" + str(L) + "_T_1.000000_random.csv"
    df = pd.read_csv("output/" + filename)
    plt.title("$<|m|>$ for T=1.0 and unordered initial spins")
    plt.xlabel("N")
    plt.ylabel("$<|m|>$")
    plt.plot(df.N, df.expected_M)
    plt.savefig("plots/burn_in/magnetization_for_L_equals_" + str(L) + ".pdf")


def plot_burn_in_times(L=20):
    """Plots burn in time for different temperatures

    Parameters
    ----------
        L : int
            Lattice size
    """
    filenames_ordered = [
        "burn_in_L_" + str(L) + "_T_1.000000_nonrandom.csv",
        "burn_in_L_" + str(L) + "_T_2.400000_nonrandom.csv",
    ]
    filenames_unordered = [
        "burn_in_L_" + str(L) + "_T_1.000000_random.csv",
        "burn_in_L_" + str(L) + "_T_2.400000_random.csv",
    ]

    temps = [1, 2.4]
    fig, axs = plt.subplots(2, 2)
    for i, (T, filename) in enumerate(zip(temps, filenames_unordered)):
        df = pd.read_csv("output/" + filename)

        axs[i][0].plot(df.N, df.expected_E)
        axs[i][0].set_title("$<\epsilon>$ for T=" + str(T) + "J / $k_B$")
        axs[i][0].set(xlabel=("N"), ylabel=("$<\epsilon>$ [J]"))

        axs[i][1].plot(df.N, df.expected_M)
        axs[i][1].set_title("$<|m|>$ for T=" + str(T) + "J / $k_B$")
        axs[i][1].set(xlabel=("N"), ylabel=("$<|m|>$ [1]"))
    fig.suptitle("Expected values for unordered initial spins")
    plt.tight_layout()
    plt.savefig("plots/burn_in/burn_in_time_unordered_L_equals_" + str(L) + ".pdf")
    plt.cla()

    fig, axs = plt.subplots(2, 2)
    for i, (T, filename) in enumerate(zip(temps, filenames_ordered)):
        df = pd.read_csv("output/" + filename)

        axs[i][0].plot(df.N, df.expected_E)
        axs[i][0].set_title("$<\epsilon>$ for T=" + str(T) + "J / $k_B$")
        axs[i][0].set(xlabel=("N"), ylabel=("$<\epsilon>$ [J]"))

        axs[i][1].plot(df.N, df.expected_M)
        axs[i][1].set_title("$<|m|>$ for T=" + str(T) + "J / $k_B$")
        axs[i][1].set(xlabel=("N"), ylabel=("$<|m|>$ [1]"))
    fig.suptitle("Expected values for ordered initial spins")
    plt.tight_layout()
    plt.savefig("plots/burn_in/burn_in_time_ordered_L_equals_" + str(L) + ".pdf")
    plt.cla()
    plt.clf()



def plot_probability_distribution():
    """Plots probability distribution for different temperatures"""
    L = "20"
    for T in ["1.0", "2.1", "2.4"]:
        df = pd.read_csv(f"output/samples_L={L}_T={T}.csv")
        plt.title(fr"Estimated probability distribution of $\epsilon$ at $T={T}J/k_B$")
        plt.xlabel(r"$\epsilon [J]$")
        plt.ylabel(fr"$p(\epsilon; {T} J / k_B)$")
        plt.hist(df.epsilon, bins="auto", density=True)
        plt.savefig(f"plots/distributions/epsilon_L={L}_T={T}.pdf")
        plt.cla()
        print(f"Variance at T={T}: {df.epsilon.var()}")
        print(f"Expected value at T={T}: {df.epsilon.mean()}")




def plot_values():
    """Plots expected values for different lattice sizes"""
    fig, axs = plt.subplots(2, 2, sharex=True)
    fig.suptitle("Plotting estimated values for different sizes of the Ising model")
    dfs = {L: pd.read_csv(f"output/values_L={L}.csv") for L in range(20, 160, 20)}
    for i, (value, ylabel, unit) in enumerate(zip(["<epsilon>", "<|m|>", "C_v", "chi"], [r"<\epsilon>", "<|m|>", "C_v", r"\chi"], ["J", "1", "k_B", "1 / J"])):
        plt.sca(axs[i // 2] [i % 2])
        for L in range(40, 160, 20):
            df = dfs[L]
            df.sort_values("T", inplace=True, ignore_index=True)
            plt.plot(df["T"], df[value], label=f"L={L}")
        plt.ylabel(f"${ylabel}$ [${unit}$]")
        plt.xlabel("T")
        plt.legend(prop={'size': 6})
    plt.tight_layout()
    plt.savefig(f"plots/values/values.pdf")
    plt.clf()


def estimate_T_inf(value: str):
    """Estimates T_inf for different lattice sizes"""
    y = []
    x = []
    label = value if value == "C_v" else r"\chi"
    for L in range(40, 160, 20):
        df = pd.read_csv(f"output/values_zoom_L={L}.csv")
        argmax = df[value].idxmax()
        y.append(df.iloc[argmax]['T'])
        x.append(1 / L)
    linear_fit = stats.linregress(x, y)
    plt.cla()
    plt.title(r"Observations of $T_c(L)$ against $L^{-1}$ and linear fit to find $T_c(\infty)$")
    plt.scatter(x, y, label=r"Observed $T_c$")
    estimate = linear_fit.intercept
    plt.plot([0] + x, estimate + linear_fit.slope * np.asarray([0] + x))
    plt.scatter([0], [estimate], s=40, label=fr"$T_c(\infty) = {estimate: .3f}$")
    plt.legend()
    plt.title(fr"Estimating $T_c(\infty)$ using ${label}$")

    plt.yticks([estimate], [f"{estimate: .3f}"])
    plt.ylabel("$T_c$ [$J / k_B$]")
    plt.xlabel("$L^{-1}$ [1]")
    plt.grid(True, color="k")
    print(f"Standard error {value}: {linear_fit.stderr}")
    plt.savefig(f"plots/T_inf/estimating_T_inf_{value}.pdf")


def main():
    plot_burn_in_times()
    plot_probability_distribution()
    plot_values()
    estimate_T_inf("C_v")
    estimate_T_inf("chi")


if __name__ == "__main__":
    main()
