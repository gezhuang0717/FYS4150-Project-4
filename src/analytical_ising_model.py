import pandas as pd
import numpy as np


class AnalyticalIsingModel:
    def __init__(self, filename: str, temperature: float, L=2):
        self._L = L
        self._N = L ** 2
        self._T = temperature
        self._beta = 1 / self._T  # Boltzman-constant?
        self._states = pd.read_csv(filename)

    def p(self, E_s):
        return 1 / self.Z * np.exp(-self._beta * E_s)

    @property
    def Z(self):
        return np.sum(np.exp(-self._beta * self._states["E(s)"]))

    @property
    def expected_E(self):
        return np.sum(self._states["E(s)"] * self.p(self._states["E(s)"]))

    @property
    def expected_epsilon(self):
        return self.expected_E / self._N

    @property
    def expected_magnitization(self):
        return np.sum(self._states["M(s)"] * self.p(self._states["E(s)"]))

    @property
    def expected_abs_M(self):
        return np.sum(np.abs(self._states["M(s)"]) * self.p(self._states["E(s)"]))

    @property
    def expected_E_squared(self):
        return np.sum((self._states["E(s)"] ** 2) * self.p(self._states["E(s)"]))

    @property
    def expected_epsilon_squared(self):
        return self.expected_E_squared / (self._N * self._N)

    @property
    def expected_m_squared(self):
        return np.sum((self._states["M(s)"] ** 2) * self.p(self._states["E(s)"])) / (
            self._N * self._N
        )


class Superanalytical:
    def __init__(self, filename: str, temperature: float, L=2):
        self._L = L
        self._N = L ** 2
        self._T = temperature
        self._beta = 1 / self._T  # Boltzman-constant?
        self._states = pd.read_csv(filename)

    def p(self, E_s):
        return 1 / self.Z * np.exp(-self._beta * E_s)

    @property
    def Z(self):
        return 4 * np.cosh(8 * self._beta) + 12

    @property
    def expected_E(self):
        return -(8 * np.sinh(8 * self._beta)) / (np.cosh(8 * self._beta) + 3)

    @property
    def expected_abs_M(self):
        return (2 * np.exp(8 * self._beta) + 4) / (np.cosh(8 * self._beta) + 3)


if __name__ == "__main__":
    aim = AnalyticalIsingModel("output/state_summary.csv", temperature=1)
    sup = Superanalytical("output/state_summary.csv", temperature=1)
    print(aim.Z - sup.Z)
    print(aim.expected_E - sup.expected_E)
    print(aim.expected_abs_M - sup.expected_abs_M)
    #  print(aim.expected_E)
    #  print(aim.expected_epsilon)
    #  print(aim.expected_epsilon_squared)
    #
    #  print(aim.expected_magnitization)
    #
    #  print(aim.expected_absolute_m)
    #  print(aim.expected_m_squared)
