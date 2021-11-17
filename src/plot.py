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

def main():
    plot_burn_in_time()

if __name__ == "__main__":
    main()