import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as sts
import subprocess

sns.set_theme()


def plot_burn_in_time(L=20):
    filenames = [
        "burn_in_L_" + str(L) + "_T_1.000000_nonrandom.csv",
        "burn_in_L_" + str(L) + "_T_1.000000_random.csv",
        "burn_in_L_" + str(L) + "_T_2.400000_nonrandom.csv",
        "burn_in_L_" + str(L) + "_T_2.400000_random.csv",
    ]
    temps = [1, 1, 2.4, 2.4]
    spin_orientation = ["ordered", "unordered", "ordered", "unordered"]

    for filename, T, orientation in zip(filenames, temps, spin_orientation):
        df = pd.read_csv("output/" + filename)
        plt.plot(df.N, df.expected_E)
        plt.title(
            "Calculated $<\epsilon>$ for T="
            + str(T)
            + ", L="
            + str(L)
            + ", and "
            + orientation
            + " initial spins"
        )
        plt.xlabel("N")
        plt.ylabel("$<\epsilon>$")
        plt.savefig(
            "plots/burn_in/expected_E" 
            + str(L) 
            + "_T_ " 
            + str(T) 
            + "_" 
            + orientation 
            + ".pdf")
        # plt.show()
        plt.cla()
        plt.plot(df.N, df.expected_M)
        plt.title(
            "Calculated <|m|> for T="
            + str(T)
            + ", L="
            + str(L)
            + ", and "
            + orientation
            + " initial spins"
        )
        plt.xlabel("N")
        plt.ylabel("<|m|>")
        plt.savefig(
            "plots/burn_in/expected_m_abs" 
            + str(L) 
            + "_T_ " 
            + str(T) 
            + "_" 
            + orientation 
            + ".pdf")
        # plt.show()
        plt.cla()


def plot_probability_distribution():
    df = pd.read_csv("output/distribution_epsilon_L=20_T=1.csv")
    print(df.epsilon)
    print(df.p)
    plt.bar(df.epsilon, df.p)
    plt.show()


def plot_values():
    for L in range(20, 180, 20):
        print(L)
        df = pd.read_csv(f"output/values_T=[2.1,2.4]_L={L}.csv")
        df.sort_values("T", inplace=True, ignore_index=True)
        df.plot(x="T", y="C_v", title=f"L={L}")
        plt.show()


def estimate_T_inf():
    y = []
    x = []
    for L in range(20, 180, 20):
        df = pd.read_csv(f"output/values_zoom_L={L}.csv")
        argmax_C_v = df.C_v.idxmax()
        argmax_chi = df.chi.idxmax()
        y.append((df.loc[argmax_C_v]['T'] + df.loc[argmax_chi]['T']) / 2)
        x.append(1 / L)
    linear_fit = sts.linregress(x, y)
    plt.scatter(x, y)
    estimate = linear_fit.intercept
    plt.plot([0] + x, estimate + linear_fit.slope * np.asarray([0] + x))
    plt.plot([0, 0], [estimate, max(y)], color="black")
    plt.yticks([estimate])
    plt.scatter([0], [estimate], s=40, label=f"{estimate: .3f}")
    plt.legend()
    plt.ylabel(r"$T_c$")
    plt.xlabel(r"$L^{-1}$")
    print(estimate)
    plt.savefig("plots/T_inf/estimating_T_inf.pdf")

def main():
    #plot_burn_in_time()
    #plot_probability_distribution()
    #plot_values_and_print_max()
    estimate_T_inf()


if __name__ == "__main__":
    main()
