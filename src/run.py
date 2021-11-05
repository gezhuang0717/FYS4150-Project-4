import argparse


def get_all_states():
    import julia

    j = julia.Julia()
    j.include("src/states.jl")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="To run the python and julia scripts")
    parser.add_argument(
        "-s",
        "--states",
        help="To get information about different states",
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
        get_all_states()
