from pandas.io.parsers import read_csv
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as sts
import subprocess

sns.set_theme()

def plot_abs_m_unordered(L=20):
    """Plots expected value for epsilon and abs(m) for T=1.0
    and unordered initial spins"""
    filename = "burn_in_L_" + str(L) + "_T_1.000000_random.csv"
    df = pd.read_csv("output/" + filename)
    plt.title("$<|m|>$ for T=1.0 and unordered initial spins")
    plt.xlabel("N")
    plt.ylabel("$<|m|>$")
    plt.plot(df.N, df.expected_M)
    plt.savefig("plots/burn_in/magnetization_for_L_equals_" + str(L) +".pdf")

def plot_burn_in_times(L=20):
    filenames_ordered = [
        "burn_in_L_" + str(L) + "_T_1.000000_nonrandom.csv",

        "burn_in_L_" + str(L) + "_T_2.400000_nonrandom.csv",
        
    ]
    filenames_unordered = [
        "burn_in_L_" + str(L) + "_T_1.000000_random.csv",
        "burn_in_L_" + str(L) + "_T_2.400000_random.csv"
    ]

    temps = [1, 2.4]
    fig, axs = plt.subplots(2, 2)
    for i, (T, filename) in enumerate(zip(temps, filenames_unordered)):
        df = pd.read_csv("output/" + filename)

        axs[i][0].plot(df.N, df.expected_E)
        axs[i][0].set_title("$<\epsilon>$ for T=" + str(T))
        axs[i][0].set(xlabel=("N"), ylabel=("$<\epsilon>$ [J/$k_B$]"))

        axs[i][1].plot(df.N, df.expected_M)
        axs[i][1].set_title("$<|m|>$ for T=" + str(T))
        axs[i][1].set(xlabel=("N"), ylabel=("$<|m|>$ [1]"))
    fig.suptitle("Expected values for unordered initial spins")
    plt.tight_layout()
    plt.savefig("plots/burn_in/burn_in_time_unordered_L_equals_" + str(L) + ".pdf")
    plt.cla()

    fig, axs = plt.subplots(2, 2)
    for i, (T, filename) in enumerate(zip(temps, filenames_ordered)):
        df = pd.read_csv("output/" + filename)

        axs[i][0].plot(df.N, df.expected_E)
        axs[i][0].set_title("$<\epsilon>$ for T=" + str(T))
        axs[i][0].set(xlabel=("N"), ylabel=("$<\epsilon>$ [J/$k_B$]"))

        axs[i][1].plot(df.N, df.expected_M)
        axs[i][1].set_title("$<|m|>$ for T=" + str(T))
        axs[i][1].set(xlabel=("N"), ylabel=("$<|m|>$ [1]"))
    fig.suptitle("Expected values for ordered initial spins")
    plt.tight_layout()
    plt.savefig("plots/burn_in/burn_in_time_ordered_L_equals_" + str(L) + ".pdf")
    plt.cla()


def plot_probability_distribution():
    L = "20"
    for T in ['1.0', '2.1', '2.4']:
        df = pd.read_csv(f"output/samples_L={L}_T={T}.csv")
        plt.title(fr"Estimated probability distribution of $\epsilon$ at T={T}")
        plt.xlabel(r"$\epsilon$")
        plt.ylabel(fr"$p(\epsilon; {T})$")
        plt.hist(df.epsilon, bins="auto", density=True)
        plt.savefig(f"plots/distributions/epsilon_L={L}_T={T}.pdf")
        plt.cla()
        print(f"Variance at T={T}: {df.epsilon.var()}")



def plot_values():
    for L in range(20, 180, 20):
        df = pd.read_csv(f"output/values_T=[2.1,2.4]_L={L}.csv")
        df.sort_values("T", inplace=True, ignore_index=True)
        df.plot(x="T", y="C_v", title=f"L={L}")
        plt.show()


def estimate_T_inf():
    y = []
    x = []
    for L in range(20, 160, 20):
        df = pd.read_csv(f"output/values_zoom_L={L}.csv")
        quad_fit_C_v = np.poly1d(np.polyfit(df["T"].to_numpy(), df["C_v"].to_numpy(), 2))
        quad_fit_chi = np.poly1d(np.polyfit(df["T"].to_numpy(), df["chi"].to_numpy(), 2))
        T_c_C_v = quad_fit_C_v.deriv().roots[0]
        T_c_chi = quad_fit_chi.deriv().roots[0]
        y.append((T_c_C_v + T_c_chi) / 2)
        x.append(1 / L)
    linear_fit = sts.linregress(x, y)
    plt.title(r"Observations of $T_c(L)$ against $L^{-1}$ and linear fit to find $T_c(\infty)$")
    plt.scatter(x, y, label="Observations")
    estimate = linear_fit.intercept
    plt.plot([0] + x, estimate + linear_fit.slope * np.asarray([0] + x))
    plt.plot([0, 0], [estimate, max(y)], color="black")
    plt.yticks([estimate], [f"{estimate: .3f}"])
    plt.xlabel([0])
    plt.scatter([0], [estimate], s=40, label=fr"$T_c(\infty) = {estimate: .3f}$")
    plt.legend()
    plt.ylabel(r"$T_c$")
    plt.xlabel(r"$L^{-1}$")
    print(estimate)
    plt.savefig("plots/T_inf/estimating_T_inf.pdf")

def main():
    plot_abs_m_unordered()
    # plot_burn_in_times()
    # plot_burn_in_times(100)
    # plot_probability_distribution()
    # plot_values_and_print_max()
    # estimate_T_inf()


if __name__ == "__main__":
    main()
