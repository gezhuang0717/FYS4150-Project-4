import pandas as pd
import numpy as np


class AnalyticalIsingModel:
    #  _N = 4
    #  _L = 2
    def __init__(self, filename: str, temperature: float, L=2):
        self.L = L
        self.N = L ** 2
        self._T = temperature
        self._beta = 1 / self._T  # Boltzman-constant?
        self._states = pd.read_csv(filename)
        #  self._Z = np.sum(np.exp(-self._beta * self._states["E(s)"]))
        #  self.expected_E = np.sum(self._states["E(s)"] * self.p(self._states["E(s)"]))
        #  self.expected_epsilon = self.expected_E / self.N

    @property
    def Z(self):
        return np.sum(np.exp(-self._beta * self._states["E(s)"]))

    @property
    def expected_E(self):
        return np.sum(self._states["E(s)"] * self.p(self._states["E(s)"]))

    @property
    def expected_M(self):
        return np.sum(np.abs(self._states["M(s)"]) * self.p(self._states["E(s)"]))

    @property
    def expected_epsilon(self):
        return self.expected_E / self.N

    @property
    def expected_epsilon(self):
        return self.expected_E / self.N

    def p(self, E_s):
        return 1 / self.Z * np.exp(-self._beta * E_s)


if __name__ == "__main__":
    aim = AnalyticalIsingModel("output/state_summary.csv", 10)
    print(aim.expected_E)
    print(aim.expected_M)
    print(aim.expected_epsilon)
