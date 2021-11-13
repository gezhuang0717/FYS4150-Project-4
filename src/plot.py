import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_burn_in_time():
    df = pd.read_csv("output/burn_in.csv")
    N = len(df.burn_in)
    # Avoid plotting tail because of large error (?)
    N_plot = int(N*4/5)
    plt.plot(df.burn_in[:N_plot], df.expected_E[:N_plot])
    plt.show()
    plt.plot(df.burn_in[:N_plot], df.expected_M[:N_plot])
    plt.show()

def main():
    plot_burn_in_time()

if __name__ == "__main__":
    main()