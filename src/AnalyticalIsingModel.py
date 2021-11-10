import pandas as pd
import numpy as np

class AnalyticalIsingModel:
    _N = 4
    _L = 2
    def __init__(self, filename: str, temprature: float):
        self._T = temprature
        self._beta = 1 / self._T # Boltzman-constant?
        self._states = pd.read_csv(filename)
        self._Z = np.sum(np.exp(-self._beta * self._states["E(s)"]))
        self.expected_epsilon = np.sum(self._states["E(s)"] * self.p(self._states["E(s)"])) / self._N

    def p(self, E_s):
        return 1/self._Z * np.exp(-self._beta * E_s)

    

if __name__ == "__main__":
    aim = AnalyticalIsingModel("output/state_summary.csv", 10)
    print(aim.expected_epsilon)


