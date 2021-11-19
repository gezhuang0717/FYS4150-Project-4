import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_burn_in_time():
    filenames = ["burn_in_T_1.000000_nonrandom.csv",
                 "burn_in_T_1.000000_random.csv",
                 "burn_in_T_2.400000_nonrandom.csv",
                 "burn_in_T_2.400000_random.csv"]
    temps = [1, 1, 2.4, 2.4]
    spin_orientation = ["ordered", "unordered", "ordered", "unordered"]
    
    for filename, T, orientation in zip(filenames, temps, spin_orientation):
        df = pd.read_csv("output/" + filename)
        plt.plot(df.N, df.expected_E)
        plt.title("Calculated $<\epsilon>$ for T=" + str(T) + " and " + orientation + " initial spins")
        plt.xlabel("N")
        plt.ylabel("$<\epsilon>$")
        plt.show()
        plt.plot(df.N, df.expected_M)
        plt.title("Calculated <|m|> for T=" + str(T) + " and " + orientation + " initial spins")
        plt.xlabel("N")
        plt.ylabel("<|m|>")
        plt.show()

def plot_probability_distribution():
    df = pd.read_csv("output/distribution_epsilon_L=20_T=1.csv")
    print(df.epsilon)
    print(df.p)
    plt.bar(df.epsilon, df.p)
    plt.show()


def plot_values():
    for L in range(20, 120, 20):
        print(L)
        df = pd.read_csv(f"output/values_L={L}.csv")
        df.sort_values("T", inplace=True)
        df.plot(x="T", y="chi", title=f"L={L}")
        plt.show()

def main():
    plot_burn_in_time()
    #plot_probability_distribution()
    #plot_values()

if __name__ == "__main__":
    main()