import numpy as np
from scipy.sparse import csc_matrix
import networkx as nx
from networkx import NetworkXError
import matplotlib.pyplot as plt
import random as rand
from random import randint
from os.path import basename, splitext

def open_graph(name):
    Grap = nx.read_edgelist(
        name, delimiter=",", create_using=nx.DiGraph(), nodetype=int
    )
    try:
        adjmat = nx.adjacency_matrix(Grap).todense()
    except NetworkXError:
        try:
            Grap = nx.read_edgelist(
                name, delimiter=" ", create_using=nx.DiGraph(), nodetype=int
            )
            adjmat = nx.adjacency_matrix(Grap).todense()
        except NetworkXError:
            try:
                Grap = nx.read_edgelist(
                    name, delimiter=":", create_using=nx.DiGraph(), nodetype=int
                )
                adjmat = nx.adjacency_matrix(Grap).todense()
            except NetworkXError:
                Grap = nx.read_edgelist(
                    name, delimiter="\t", create_using=nx.DiGraph(), nodetype=int
                )
                adjmat = nx.adjacency_matrix(Grap).todense()
    nodes = len(Grap.nodes())
    density = Grap.size() / (nodes * (nodes - 1) * 0.5)
    return nx.info(Grap), density, adjmat


def adj2trans(adjmat):
    """
    Converts adjacency matrix into transition matrix

    Parameters
    ----------
    adjmat: adjacency matrix
    """
    n = adjmat.shape[0]
    P_mat = np.zeros((n, n))
    # converting adjmat to transmat
    for i in range(0, n):
        if np.sum(adjmat[i])!=0:
            P_mat[i] = adjmat[i] / np.sum(adjmat[i])
    return P_mat


def page_rank(P, s=0.8, maxerr=0.0001):
    """
    Computes the pagerank for each of the n states

    Parameters
    ----------
    P: matrix representing state transitions
       Gij is a binary value representing a transition from state i to j.

    s: probability of following a transition. 1-s probability of teleporting
       to another state.

    maxerr: if the sum of pageranks between iterations is bellow this we will
            have converged.
    """
    P = adj2trans(P)
    n = P.shape[0]
    r = np.ones(n) / n
    ro = np.ones(n)
    G = np.ones((n, n))
    while np.sum(np.abs(r - ro)) > maxerr:
        ro = r.copy()
        r = (
            s * r.dot(P)
            + (1 - s) * r.dot(G) / n
)
    # return normalized pagerank
    return np.argsort(-r)


def epidemic_spread(
    adjmat,
    type="random",
    iter=150,
    infec=0.05,
    ngb_contam=0.2,
    heal=0.02,
    vac=0.22,
    alpha=0.8,
):
    """
    Simulates the epidemic spread using either a random infection

    Parameters:
    -----------
    adjmat: input graph
    type: no vaccination, random, or pagerank optimized simulation
    iter: Number of iterations of the simulation
    infec: percent of randomly infected people in the initial state
    ngb_contam: Probability of contaminating each neighbor when infected
    heal: Probability of randomly healing when infected
    vac: Ratio of initially vaccinated individuals
    """

    n = adjmat.shape[0]
    # Grap = get_infection_vect(Grap,infec)
    if type == "pagerank":
        # outputs which nodes to vaccinate using pagerank optimization
        init_vac = page_rank(adjmat)[: int(n * vac)]
    else:  # random vaccination
        init_vac = np.random.choice(n, int(n * vac), replace=False)

    # delete vaccinated people
    adjmat = np.delete(adjmat, init_vac, 0)
    adjmat = np.delete(adjmat, init_vac, 1)
    n = adjmat.shape[0]
    contam_vector = np.random.choice(n, int(n * infec))
    nb_infect = []  # to be plotted with enum
    t = []
    # initially the infection ratio is the parameter given
    res = [infec]
    for i in range(iter):

        # get neighbors of contaminated people
        ngbs = adjmat[contam_vector, :].nonzero()[1]
        # infect neighbors
        contam_vector = np.unique(
            np.concatenate(
                (
                    contam_vector,
                    np.extract(np.random.rand(ngbs.shape[0]) < ngb_contam, ngbs),
                )
            )
        )
        # non neighbors get infected with probability 1-alpha
    # non_ngbs = np.array(list(set(np.arange(0, n)) - set(contam_vector) - set(ngbs)))
    # if False and len(non_ngbs) > 0:
    #     contam_vector = np.unique(
    #         np.concatenate(
    #             (
    #                 contam_vector,
    #                 np.extract(
    #                     np.random.rand(non_ngbs.shape[0]) < 1 - alpha, non_ngbs
    #                 ),
    #             )
    #         )
    #     )

        # probability to heal
        contam_vector = np.delete(
            contam_vector, np.where(np.random.rand(contam_vector.shape[0]) < heal)[0]
        )
        res.append(float(len(contam_vector) / n))
    return list(enumerate(res))


if __name__ == "__main__":
    print(epidemic_spread("tt", type="pagerank"))
    plt.plot(*zip(*epidemic_spread("tt.csv", vac=0)), label="No vaccination")
    plt.plot(*zip(*epidemic_spread("tt.csv")), label="Random Vaccination")
    plt.plot(
        *zip(*epidemic_spread("tt.csv", type="pagerank")), label="PageRank Vaccination"
    )
    plt.xlabel("Steps")
    plt.ylabel("Percent of infected individuals")
    plt.title(
        "Epidemic spread simulation using infection vector with different parameters"
    )
    plt.title("Sim for "+args.filename, loc=left)
    plt.legend()
    plt.savefig("../plots/"+basename(args.filename)+".png")
