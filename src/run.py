import argparse
import itertools
import numpy as np
from collections import defaultdict
import plot
import analytical_ising_model
import zoom

def get_all_states(L):
    """
    List all possible states of a LxL lattice (equivelant states grouped)
    """
    N = L ** 2

    spins = itertools.product([-1, 1], repeat=N)

    state_summary = []

    for spin in spins:
        spin = np.array(spin).reshape(L, L)

        positive_spins = np.sum(spin == 1)

        E_h = np.multiply(spin, np.roll(spin, 1, axis=0))
        E_v = np.multiply(spin, np.roll(spin, 1, axis=1))
        E_s = -(np.sum(E_h) + np.sum(E_v))

        M_s = np.sum(spin)

        state_summary.append([positive_spins, E_s, M_s])

    # extract the unique states and number of degeneracies
    unique_states = np.unique(state_summary, axis=0)
    degeneracies = defaultdict(int)
    for state in state_summary:
        degeneracies[str(state)] += 1

    # combine degeneracies and unique states
    states = []
    for state in unique_states:
        states.append([*state, degeneracies[str(list(state))]])

    # write state summary to csv file with header
    np.savetxt(
        "output/state_summary.csv",
        sorted(states, reverse=True),
        delimiter=",",
        fmt="%s",
        header="positive spins,E(s),M(s),degeneracy",
        comments="",
    )

    print("Successfully saved state summary to file")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="To run the python and julia scripts")
    parser.add_argument(
        "-s",
        "--states",
        help="To get information about different states",
        action="store_true",
    )
    parser.add_argument(
        "-p",
        "--plot",
        help="To generate plots",
        action="store_true",
    )
    parser.add_argument(
        "-an",
        "--analytical",
        help="To generate analytical values",
        action="store_true",
    )
    parser.add_argument(
        "-z",
        "--zoom",
        help="Find maximum values and zoom",
        action="store_true",
    )
    parser.add_argument(
        "-a",
        "--all",
        help="To run everything",
        action="store_true",
    )

    args = parser.parse_args()

    if not any(vars(args).values()):
        parser.print_help()
    if args.states or args.all:
        get_all_states(2)
    if args.plot or args.all:
        plot.main()
    if args.analytical or args.all:
        analytical_ising_model.main()
    if args.zoom or args.all:
        zoom.main()
