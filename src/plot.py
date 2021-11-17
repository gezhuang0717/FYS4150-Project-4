import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_burn_in_time():
    df = pd.read_csv("output/burn_in.csv")
    # Avoid plotting tail because of large error (?
    plt.plot(df.N, df.expected_E)
    plt.show()
    plt.plot(df.N, df.expected_M)
    plt.show()

def plot_probability_distribution():
    df = pd.read_csv("output/distribution_epsilon_L=20_T=1.csv")
    print(df.epsilon)
    print(df.p)
    df.plot.hist(bins=5)
    #plt.hist(df.p, df.epsilon)
    plt.show()


def main():
    #plot_burn_in_time()
    plot_probability_distribution()

if __name__ == "__main__":
    main()