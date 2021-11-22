"""
Take the results from initial scanning at a large interval of temperatures,
and zoom in at a neighbourhood containting the values for T when C_v and
chi attains their maximum.
"""

import pandas as pd
import subprocess


def zoom():
    """
    Take the results from initial scanning at a large interval of temperatures,
    and zoom in at a neighbourhood containting the values for T when C_v and
    chi attains their maximum.
    """
    seed = 9642
    for L in range(20, 160, 20):
        print(L)
        df = pd.read_csv(f"output/values_L={L}.csv")
        df.sort_values("T", inplace=True, ignore_index=True)
        margin = 2
        argmax_C_v = df.C_v.idxmax()
        argmax_chi = df.chi.idxmax()
        T_min = df.loc[max(0, min(argmax_C_v, argmax_chi) - margin)]["T"]
        T_max = df.loc[min(max(argmax_C_v, argmax_chi) + margin, len(df["T"]) - 1)]["T"]
        print(
            f"L = {L} - argmax C_v: {df.loc[argmax_C_v]['T']}, argmax chi {df.loc[argmax_chi]['T']}"
        )
        print(f"Look between temperatures {T_min} and {T_max}")
        subprocess.run(["./runner", "-z", str(L), str(T_min), str(T_max), str(seed)])
        seed += 1


def main():
    zoom()


if __name__ == "__main__":
    main()
