import pandas as pd
import numpy as np


class AnalyticalIsingModel:
    def __init__(self, filename: str, temperature: float, L=2):
        self._L = L
        self._N = L ** 2
        self._T = temperature
        self._beta = 1 / self._T  # Boltzman-constant?
        self._states = pd.read_csv(filename)

        np_state = self._states.to_numpy()
        np_state = np_state.transpose()

        self.positive_spins = np_state[0]
        self.E_s = np_state[1]
        self.M_s = np_state[2]
        self.degeneracy = np_state[3]

    def p(self, E_s):
        return 1 / self.Z * np.exp(-self._beta * E_s)

    @property
    def Z(self):
        return np.sum(np.exp(-self._beta * self.E_s) * self.degeneracy)

    @property
    def expected_E(self):
        return np.sum(self.degeneracy * self.E_s * self.p(self.E_s))

    @property
    def expected_epsilon(self):
        return self.expected_E / self._N

    @property
    def expected_E_squared(self):
        return np.sum((self.E_s ** 2) * self.degeneracy * self.p(self.E_s))

    @property
    def expected_epsilon_squared(self):
        return self.expected_E_squared / (self._N * self._N)

    @property
    def expected_abs_M(self):
        return np.sum(abs(self.M_s) * self.degeneracy * self.p(self.E_s))

    @property
    def expected_abs_m(self):
        return self.expected_abs_M / self._N

    @property
    def expected_M_squared(self):
        return np.sum((self.M_s ** 2) * self.degeneracy * self.p(self.E_s))
        #  return np.sum(self._states["M(s)"] ** 2 * self.p(self._states["E(s)"]))

    @property
    def expected_m_squared(self):
        return self.expected_M_squared / (self._N * self._N)

    @property
    def C_v(self):
        return (
            (1 / self._N)
            * (1 / (self._T * self._T))
            * (self.expected_E_squared - self.expected_E * self.expected_E)
        )

    @property
    def chi(self):
        return (
            (1 / self._N)
            * (1 / self._T)
            * (self.expected_M_squared - self.expected_abs_M * self.expected_abs_M)
        )


class Superanalytical:
    def __init__(self, filename: str, temperature: float, L=2):
        assert L == 2, "Only implemented for L=2 (N=4)"
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
    def expected_epsilon(self):
        return -(2 * np.sinh(8 * self._beta)) / (np.cosh(8 * self._beta) + 3)

    @property
    def expected_E(self):
        return self._N * self.expected_epsilon

    @property
    def expected_epsilon_squared(self):
        return 4 * np.cosh(8 * self._beta) / (np.cosh(8 * self._beta) + 3)

    @property
    def expected_E_squared(self):
        return self._N ** 2 * self.expected_epsilon_squared

    @property
    def expected_abs_m(self):
        return (np.exp(8 * self._beta) + 2) / (2 * np.cosh(8 * self._beta) + 6)

    @property
    def expected_abs_M(self):
        return self._N * self.expected_abs_m

    @property
    def expected_m_squared(self):
        return (np.exp(8 * self._beta) + 1) / (2 * np.cosh(8 * self._beta) + 6)

    @property
    def expected_M_squared(self):
        return self._N ** 2 * self.expected_m_squared

    @property
    def C_v(self):
        return (
            16
            / (self._T ** 2)
            * (1 + 3 * np.cosh(8 * self._beta))
            / ((np.cosh(8 * self._beta) + 3) ** 2)
        )

    @property
    def C_v_naive(self):
        return (
            (1 / self._N)
            * (1 / (self._T * self._T))
            * (self.expected_E_squared - self.expected_E * self.expected_E)
        )

    @property
    def chi(self):
        return (
            (1 / self._N)
            * (1 / self._T)
            * (
                (8 * np.exp(8 * self._beta) + 8) / (np.cosh(8 * self._beta) + 3)
                - (
                    4
                    * (np.exp(8 * self._beta) + 2) ** 2
                    / (np.cosh(8 * self._beta) + 3) ** 2
                )
            )
        )

    @property
    def chi_naive(self):
        return (
            (1 / self._N)
            * (1 / self._T)
            * (self.expected_M_squared - self.expected_abs_M * self.expected_abs_M)
        )


def test_analytical_ising_model():
    aim = AnalyticalIsingModel("output/state_summary.csv", temperature=1)
    sup = Superanalytical("output/state_summary.csv", temperature=1)

    assert aim.Z == sup.Z, "Z-values are not equal"
    assert (
        aim.expected_epsilon == sup.expected_epsilon
    ), "Expected epsilon-values are not equal"
    assert aim.expected_E == sup.expected_E, "Expected E-values are not equal"
    assert (
        aim.expected_epsilon_squared == sup.expected_epsilon_squared
    ), "Expected epsilon^2-values are not equal"
    assert (
        aim.expected_E_squared == sup.expected_E_squared
    ), "Expected E^2-values are not equal"

    assert aim.expected_abs_m == sup.expected_abs_m, "Expected |m|-values are not equal"
    assert aim.expected_abs_M == sup.expected_abs_M, "Expected |M|-values are not equal"

    assert (
        aim.expected_m_squared == sup.expected_m_squared
    ), "Expected m^2-values are not equal"
    assert (
        aim.expected_M_squared == sup.expected_M_squared
    ), "Expected M^2-values are not equal"
    assert abs(sup.chi - sup.chi_naive) < 1e-14, "Chi-values are not equal"

    assert (
        abs(sup.C_v - sup.C_v_naive) < 1e-14
    ), "C_v analytical value and formula value are not equal"
    assert abs(aim.C_v - sup.C_v) < 1e-14, "C_v-values are not equal"
    assert abs(aim.chi - sup.chi) < 1e-14, "chi-values are not equal"


def main():
    with open("output/analytical_L=2.csv", "w") as outfile:
        outfile.write("T,<epsilon>,<|m|>,C_v,chi\n")
        for T in np.arange(1, 2.5, 0.1):
            aim = AnalyticalIsingModel("output/state_summary.csv", temperature=T)
            outfile.write(
                f"{T},{aim.expected_epsilon},{aim.expected_abs_m},{aim.C_v},{aim.chi}\n"
            )

    test_analytical_ising_model()


if __name__ == "__main__":
    main()
