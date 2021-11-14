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
    def expected_M_squared(self):
        return np.sum(self._states["M(s)"] ** 2 * self.p(self._states["E(s)"]))

    @property
    def expected_abs_M(self):
        return np.sum(np.abs(self._states["M(s)"]) * self.p(self._states["E(s)"]))

    @property
    def expected_abs_m(self):
        return self.expected_abs_M / self._N

    @property
    def expected_E_squared(self):
        return np.sum((self._states["E(s)"] ** 2) * self.p(self._states["E(s)"]))

    @property
    def expected_epsilon_squared(self):
        return self.expected_E_squared / (self._N * self._N)

    @property
    def C_v(self):
        return (1 / self._N) * (1 / (self._T * self._T)) * (self.expected_E_squared - self.expected_E * self.expected_E) 
    
    @property
    def chi(self):
        return (1 / self._N) * (1 / self._T) * (self.expected_M_squared - self.expected_abs_M * self.expected_abs_M)


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

    with open("output/expected_epsilon_analytical_L=2.csv", "w") as epsilon_file, \
    open("output/expected_m_abs_analytical_L=2.csv", "w") as m_file, \
    open("output/C_v_analytical_L=2.csv", "w") as C_v_file, \
    open("output/chi_analytical_L=2.csv", "w") as chi_file:
        epsilon_file.write("T,<epsilon>\n")
        m_file.write("T,<|m|>\n")
        C_v_file.write("T,C_v\n")
        chi_file.write("T, chi\n")
        for T in np.arange(1, 2.5, .1):
            aim = AnalyticalIsingModel("output/state_summary.csv", temperature=T)
            epsilon_file.write(f"{T},{aim.expected_epsilon}\n")
            m_file.write(f"{T},{aim.expected_abs_m}\n")
            C_v_file.write(f"{T},{aim.C_v}\n")
            chi_file.write(f"{T},{aim.chi}\n")





    aim = AnalyticalIsingModel("output/state_summary.csv", temperature=1)
    sup = Superanalytical("output/state_summary.csv", temperature=1)
    print("")
    print(aim.expected_epsilon)
    print(aim.expected_abs_m)
    print(aim.C_v)
    print(aim.chi)
    #  print(aim.expected_E)
    #  print(aim.expected_epsilon)
    #  print(aim.expected_epsilon_squared)
    #
    #  print(aim.expected_magnitization)
    #
    #  print(aim.expected_absolute_m)
    #  print(aim.expected_m_squared)
