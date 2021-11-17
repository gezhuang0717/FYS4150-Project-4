import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_burn_in_time():
    filenames = ["burn_in_T_1.000000_nonrandom.csv",
                 "burn_in_T_1.000000_random.csv",
                 "burn_in_T_2.400000_nonrandom.csv",
                 "burn_in_T_2.400000_random.csv"]
    temps = [1, 1, 2.4, 2.4]
    spin_orientation = ["uniform", "random", "uniform", "random"]
    
    for filename, T, orientation in zip(filenames, temps, spin_orientation):
        df = pd.read_csv("output/" + filename)
        plt.plot(df.N, df.expected_E)
        plt.title("Calculated <E> for T=" + str(T) + " and " + orientation + " initial spins")
        plt.xlabel("N")
        plt.ylabel("<E>")
        plt.show()
        plt.plot(df.N, df.expected_M)
        plt.title("Calculated <|M|> for T=" + str(T) + " and " + orientation + " initial spins")
        plt.xlabel("N")
        plt.ylabel("<M>")
        plt.show()

def plot_probability_distribution():
    df = pd.read_csv("output/distribution_epsilon_L=20_T=1.csv")
    print(df.epsilon)
    print(df.p)
    df.plot.hist(bins=5)
    #plt.hist(df.p, df.epsilon)
    plt.show()


def main():
    plot_burn_in_time()
    # plot_probability_distribution()

if __name__ == "__main__":
    main()