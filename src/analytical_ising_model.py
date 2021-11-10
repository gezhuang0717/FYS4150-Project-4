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
    def expected_absolute_m(self):
        return (
            np.sum(np.abs(self._states["M(s)"]) * self.p(self._states["E(s)"]))
            / self._N
        )

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


if __name__ == "__main__":
    aim = AnalyticalIsingModel("output/state_summary.csv", temperature=10)
    print(aim.expected_E)
    print(aim.expected_epsilon)
    print(aim.expected_epsilon_squared)

    print(aim.expected_absolute_m)
    print(aim.expected_m_squared)
