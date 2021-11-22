"""
Analytical solution of the IsingModel for the 2x2 case

Use to test against the nummerical implementation
"""

import pandas as pd
import numpy as np


class NumericalIsingModel:
    def __init__(self, filename: str, temperature: float, L=2):
        """Initialize the numerical ising model

        Parameters
        ----------
            filename : str
                The filename of the data
            temperature : float
                The temperature of the system
            L : int
                The size of the lattice
        """
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
        """Calculate the probability of the energy E(s)

        Parameters
        ----------
            E_s : float
                The energy of the state

        Returns
        -------
            p : float
                The probability of the energy E(s)
        """
        return 1 / self.Z * np.exp(-self._beta * E_s)

    @property
    def Z(self):
        """Calculate the partition function

        Returns
        -------
            Z : float
                The partition function
        """
        return np.sum(np.exp(-self._beta * self.E_s) * self.degeneracy)

    @property
    def expected_E(self):
        """Calculate the expected energy

        Returns
        -------
            E : float
                The expected energy
        """
        return np.sum(self.degeneracy * self.E_s * self.p(self.E_s))

    @property
    def expected_epsilon(self):
        """Calculate the expected energy per spin

        Returns
        -------
            epsilon : float
                The expected energy
        """
        return self.expected_E / self._N

    @property
    def expected_E_squared(self):
        """Calculate the expected energy squared

        Returns
        -------
            E_squared : float
                The expected energy squared
        """
        return np.sum((self.E_s ** 2) * self.degeneracy * self.p(self.E_s))

    @property
    def expected_epsilon_squared(self):
        """Calculate the expected energy squared per spin

        Returns
        -------
            epsilon_squared : float
                The expected energy squared
        """
        return self.expected_E_squared / (self._N * self._N)

    @property
    def expected_abs_M(self):
        """Calculate the expected absolute magnetization

        Returns
        -------
            M : float
                The expected absolute magnetization
        """
        return np.sum(abs(self.M_s) * self.degeneracy * self.p(self.E_s))

    @property
    def expected_abs_m(self):
        """Calculate the expected absolute magnetization per spin

        Returns
        -------
            m : float
                The expected absolute magnetization per spin
        """
        return self.expected_abs_M / self._N

    @property
    def expected_M_squared(self):
        """Calculate the expected magnetization squared

        Returns
        -------
            M_squared : float
                The expected magnetization squared
        """
        return np.sum((self.M_s ** 2) * self.degeneracy * self.p(self.E_s))

    @property
    def expected_m_squared(self):
        """Calculate the expected magnetization squared per spin

        Returns
        -------
            m_squared : float
                The expected magnetization squared per spin
        """
        return self.expected_M_squared / (self._N * self._N)

    @property
    def C_v(self):
        """Calculate the specific heat capacity

        Returns
        -------
            C_v : float
                The specific heat capacity
        """
        return (
            (1 / self._N)
            * (1 / (self._T * self._T))
            * (self.expected_E_squared - self.expected_E * self.expected_E)
        )

    @property
    def chi(self):
        """Calculate the susceptibility

        Returns
        -------
            chi : float
                The susceptibility
        """
        return (
            (1 / self._N)
            * (1 / self._T)
            * (self.expected_M_squared - self.expected_abs_M * self.expected_abs_M)
        )


class AnalyticalIsingModel:
    def __init__(self, filename: str, temperature: float, L=2):
        """Initialize the analytical ising model

        Parameters
        ----------
            filename : str
                The filename of the data
            temperature : float
                The temperature of the system
            L : int
                The size of the lattice
        """
        assert L == 2, "Only implemented for L=2 (N=4)"
        self._L = L
        self._N = L ** 2
        self._T = temperature
        self._beta = 1 / self._T  # Boltzman-constant?
        self._states = pd.read_csv(filename)

    def p(self, E_s):
        """Calculate the probability of the energy E(s)

        Parameters
        ----------
            E_s : float
                The energy of the state

        Returns
        -------
            p : float
                The probability of the energy E(s)
        """
        return 1 / self.Z * np.exp(-self._beta * E_s)

    @property
    def Z(self):
        """Calculate the partition function

        Returns
        -------
            Z : float
                The partition function
        """
        return 4 * np.cosh(8 * self._beta) + 12

    @property
    def expected_epsilon(self):
        """Calculate the expected energy per spin

        Returns
        -------
            epsilon : float
                The expected energy
        """
        return -(2 * np.sinh(8 * self._beta)) / (np.cosh(8 * self._beta) + 3)

    @property
    def expected_E(self):
        """Calculate the expected energy

        Returns
        -------
            E : float
                The expected energy
        """
        return self._N * self.expected_epsilon

    @property
    def expected_epsilon_squared(self):
        """Calculate the expected energy squared per spin

        Returns
        -------
            epsilon_squared : float
                The expected energy squared
        """
        return 4 * np.cosh(8 * self._beta) / (np.cosh(8 * self._beta) + 3)

    @property
    def expected_E_squared(self):
        """Calculate the expected energy squared

        Returns
        -------
            E_squared : float
                The expected energy squared
        """
        return self._N ** 2 * self.expected_epsilon_squared

    @property
    def expected_abs_m(self):
        """Calculate the expected absolute magnetization per spin

        Returns
        -------
            m : float
                The expected absolute magnetization per spin
        """
        return (np.exp(8 * self._beta) + 2) / (2 * np.cosh(8 * self._beta) + 6)

    @property
    def expected_abs_M(self):
        """Calculate the expected absolute magnetization

        Returns
        -------
            M : float
                The expected absolute magnetization
        """
        return self._N * self.expected_abs_m

    @property
    def expected_m_squared(self):
        """Calculate the expected magnetization squared per spin

        Returns
        -------
            m_squared : float
                The expected magnetization squared per spin
        """
        return (np.exp(8 * self._beta) + 1) / (2 * np.cosh(8 * self._beta) + 6)

    @property
    def expected_M_squared(self):
        """Calculate the expected magnetization squared

        Returns
        -------
            M_squared : float
                The expected magnetization squared
        """
        return self._N ** 2 * self.expected_m_squared

    @property
    def C_v(self):
        """Calculate the specific heat capacity

        Returns
        -------
            C_v : float
                The specific heat capacity
        """
        return (
            16
            / (self._T ** 2)
            * (1 + 3 * np.cosh(8 * self._beta))
            / ((np.cosh(8 * self._beta) + 3) ** 2)
        )

    @property
    def C_v_naive(self):
        """Calculate the specific heat capacity, using the analytical expressions for expected_E_squared and expected_E

        Returns
        -------
            C_v : float
                The specific heat capacity
        """
        return (
            (1 / self._N)
            * (1 / (self._T * self._T))
            * (self.expected_E_squared - self.expected_E * self.expected_E)
        )

    @property
    def chi(self):
        """Calculate the susceptibility

        Returns
        -------
            chi : float
                The susceptibility
        """
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
        """Calculate the susceptibility, using the analytical expressions for expected_M_squared and expected_abs_M

        Returns
        -------
            chi : float
                The susceptibility
        """
        return (
            (1 / self._N)
            * (1 / self._T)
            * (self.expected_M_squared - self.expected_abs_M * self.expected_abs_M)
        )


def test_analytical_ising_model():
    """Test the analytical ising model against numerical ising model"""
    numerical_ising_model = NumericalIsingModel(
        "output/state_summary.csv", temperature=1
    )
    analytical_ising_model = AnalyticalIsingModel(
        "output/state_summary.csv", temperature=1
    )

    assert numerical_ising_model.Z == analytical_ising_model.Z, "Z-values are not equal"

    assert (
        numerical_ising_model.expected_epsilon
        == analytical_ising_model.expected_epsilon
    ), "Expected epsilon-values are not equal"
    assert (
        numerical_ising_model.expected_E == analytical_ising_model.expected_E
    ), "Expected E-values are not equal"

    assert (
        numerical_ising_model.expected_epsilon_squared
        == analytical_ising_model.expected_epsilon_squared
    ), "Expected epsilon^2-values are not equal"
    assert (
        numerical_ising_model.expected_E_squared
        == analytical_ising_model.expected_E_squared
    ), "Expected E^2-values are not equal"

    assert (
        numerical_ising_model.expected_abs_m == analytical_ising_model.expected_abs_m
    ), "Expected |m|-values are not equal"
    assert (
        numerical_ising_model.expected_abs_M == analytical_ising_model.expected_abs_M
    ), "Expected |M|-values are not equal"

    assert (
        numerical_ising_model.expected_m_squared
        == analytical_ising_model.expected_m_squared
    ), "Expected m^2-values are not equal"
    assert (
        numerical_ising_model.expected_M_squared
        == analytical_ising_model.expected_M_squared
    ), "Expected M^2-values are not equal"

    assert (
        abs(analytical_ising_model.C_v - analytical_ising_model.C_v_naive) < 1e-14
    ), "C_v analytical value and formula value are not equal"
    assert (
        abs(numerical_ising_model.C_v - analytical_ising_model.C_v) < 1e-14
    ), "C_v-values are not equal"

    assert (
        abs(numerical_ising_model.chi - analytical_ising_model.chi) < 1e-14
    ), "chi-values are not equal"
    assert (
        abs(analytical_ising_model.chi - analytical_ising_model.chi_naive) < 1e-14
    ), "Chi-values are not equal"


def main():
    """
    Verify that the numerial and analytical expressions give the same
    results and write the analytical values to a file
    """
    test_analytical_ising_model()

    with open("output/analytical_L=2.csv", "w") as outfile:
        outfile.write("T,<epsilon>,<|m|>,C_v,chi\n")
        for T in np.arange(1, 2.5, 0.1):
            analytical_ising_model = AnalyticalIsingModel(
                "output/state_summary.csv", temperature=T
            )
            outfile.write(
                f"{T},{analytical_ising_model.expected_epsilon},{analytical_ising_model.expected_abs_m},{analytical_ising_model.C_v},{analytical_ising_model.chi}\n"
            )


if __name__ == "__main__":
    main()
