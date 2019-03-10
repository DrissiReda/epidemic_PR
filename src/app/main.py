import argparse
import os
import sys
import csv
import numpy as np
import datetime
import math
import time
from app.PR_sim import *


def human_time(ms):
    ms = int(ms*1000)
    x = ms%1000
    left = int(ms - x)
    seconds = int(left/1000) % 60
    left = int(left/1000) - seconds
    minutes = int(left/60) %60
    left = int(left/60) - minutes
    hours = int(left/60) % 60
    res=""
    if hours != 0:
        res+=str(hours) + " hours "
    if minutes != 0:
        res+=str(minutes)+ " minutes "
    if seconds != 0:
        res+=str(seconds)+ " seconds "
    if x!=0:
        res+=str(x) + " miliseconds"
    return  res



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simple epidemic modeling example using PageRank"
    )
    parser.add_argument(
        "-g",
        "--input-graph",
        dest="filename",
        type=str,
        action="store",
        help="File containing the edge list of the matrix to process",
    )
    parser.add_argument(
        "-i",
        "--infected",
        dest="infected",
        type=float,
        action="store",
        default=0.05,
        help="Ratio of initially randomly individuals",
    )
    parser.add_argument(
        "-v",
        "--vaccinate",
        dest="random_vaccination",
        type=float,
        action="store",
        default=0.12,
        help="Ratio of initially randomly vaccinated individuals",
    )
    parser.add_argument(
        "-H",
        "--heal",
        dest="heal",
        type=float,
        action="store",
        default=0.26,
        help="Probability of randomly healing when infected",
    )
    parser.add_argument(
        "-c",
        "--contamination",
        dest="contamination",
        type=float,
        action="store",
        default=0.2,
        help="Probability of contaminating each neighbor when infected",
    )
    parser.add_argument(
        "-t",
        "--iteration",
        dest="iter",
        type=int,
        action="store",
        default=150,
        help="Number of iterations of the simulation",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    info, density, adjmat = open_graph(args.filename)
    print(
        str(info)
        + "\nDensity:"
        + str(density)
        + "\nRatio of initially infected individuals:"
        + str(args.infected)
        + "\nRatio of initially randomly vaccinated individuals:"
        + str(args.random_vaccination)
        + "\nProbability of contaminating each neighbor when infected:"
        + str(args.contamination)
        + "\nProbability of randomly healing when infected:"
        + str(args.heal)
        + "\nNumber of iterations:"
        + str(args.iter)
    )
    start = time.time()
    NoVac = zip(
        *epidemic_spread(
            adjmat,
            iter=args.iter,
            infec=args.infected,
            ngb_contam=args.contamination,
            heal=args.heal,
            vac=0,
        )
    )
    NoVac_time = time.time() - start
    start = time.time()
    RandomVac = zip(
        *epidemic_spread(
            adjmat,
            iter=args.iter,
            infec=args.infected,
            ngb_contam=args.contamination,
            heal=args.heal,
            vac=args.random_vaccination,
        )
    )
    RandomVac_time = time.time() - start
    start = time.time()
    PRVac = zip(
        *epidemic_spread(
            adjmat,
            type="pagerank",
            iter=args.iter,
            infec=args.infected,
            ngb_contam=args.contamination,
            heal=args.heal,
        )
    )
    PRVac_time = time.time() - start
    print(
        "Time spent to run each simulation:"
        + "\nNo Vaccination: "
        + str(human_time(NoVac_time))
        + "\nRandom Vaccination: "
        + str(human_time(RandomVac_time))
        + "\nPageRank Vaccination: "
        + str(human_time(PRVac_time))
    )
    plt.plot(*NoVac, label="No vaccination")
    plt.plot(*RandomVac, label="Random Vaccination")
    plt.plot(*PRVac, label="PageRank Vaccination")

    plt.xlabel("Steps")
    plt.ylabel("Percent of infected individuals")
    plt.title(
        "Epidemic spread simulation using infection vector\n with different parameters\n"
        +"Simulating "+basename(args.filename), loc="left"
    )
    plt.legend()
    plt.savefig("../plots/"+splitext(basename(args.filename))[0]+".png", bbox_inches='tight')
