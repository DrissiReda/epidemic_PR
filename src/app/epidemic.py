import numpy as np
import csv
import matplotlib.pyplot as plt


def page_rank(transition_matrix: np.ndarray, alpha: float = 0.85, error: float = 0.001):
    size = transition_matrix.shape[0]
    x1 = np.ones(size) / size
    x2 = np.ones(size)
    google_matrix = np.ones((size, size))
    while np.sum(np.abs(x1 - x2)) > error:
        x2 = x1.copy()
        x1 = (
            alpha * x1.dot(transition_matrix)
            + (1 - alpha) * x1.dot(google_matrix) / size
        )
    return np.argsort(-x1)


def transition_from_adjacency(adjacency_matrix: np.ndarray):
    ri, ci = adjacency_matrix.nonzero()
    transition_matrix = np.zeros(adjacency_matrix.shape)
    transition_matrix[ri, ci] = adjacency_matrix[ri, ci] / adjacency_matrix.sum(1)[ri]
    return transition_matrix


def adjacency_from_csv(path: str, delimiter: str = " ", header: bool = False):
    with open(path) as file:
        reader = csv.reader(file, delimiter=delimiter)
        if header:
            next(reader, None)
        edges = [(int(row[0]), int(row[1])) for row in reader]
    size = max(max(v1, v2) for v1, v2 in edges) + 1
    adjacency_matrix = np.zeros((size, size), dtype=np.int8)
    for (v1, v2) in edges:
        adjacency_matrix[v1, v2] = adjacency_matrix[v2, v1] = 1
    return adjacency_matrix


def iterate(
    adjacency_matrix: np.ndarray,
    initially_infected_ratio=0.05,
    contamination_probability=0.2,
    cure_probability=0.24,
    vaccinated_ratio=0.2,
    optimize_vaccination=False,
    iteration_count=150,
):
    original_size = adjacency_matrix.shape[0]
    if optimize_vaccination:
        vaccinated = page_rank(transition_from_adjacency(adjacency_matrix))[
            : int(original_size * vaccinated_ratio)
        ]
    else:
        vaccinated = np.random.choice(
            original_size, int(original_size * vaccinated_ratio), replace=False
        )
    adjacency_matrix = np.delete(adjacency_matrix, vaccinated, 0)
    adjacency_matrix = np.delete(adjacency_matrix, vaccinated, 1)
    contaminated = np.random.choice(
        adjacency_matrix.shape[0],
        int(original_size * initially_infected_ratio),
        replace=False,
    )
    result = [(0, len(contaminated) / original_size)]
    for i in range(1, iteration_count):
        # print(i)
        contaminated = np.delete(
            contaminated,
            np.where(np.random.rand(contaminated.shape[0]) < cure_probability)[0],
        )
        adjacent_to_contaminated = adjacency_matrix[contaminated, :].nonzero()[1]
        newly_contaminated = np.extract(
            np.random.rand(adjacent_to_contaminated.shape[0])
            < contamination_probability,
            adjacent_to_contaminated,
        )
        contaminated = np.concatenate((contaminated, newly_contaminated))
        contaminated = np.unique(contaminated)
        result += [(i, len(contaminated) / original_size)]
    return result


if __name__ == "__main__":
    A = adjacency_from_csv(
        "tt.csv",
        # "data/graph.csv",
        # "data/facebook_combined.txt",
        delimiter=",",
        header=True,
    )
    P = transition_from_adjacency(A)
    print(iterate(A, optimize_vaccination=True))
    plt.plot(*zip(*iterate(A, vaccinated_ratio=0)), label="No Vaccination")
    plt.plot(*zip(*iterate(A)), label="Random Vaccination")
    plt.plot(*zip(*iterate(A, optimize_vaccination=True)), label="PageRank Vaccination")
    plt.xlabel("Time (iteration)")
    plt.ylabel("Proportion of infected individuals")
    plt.title("Stochastic simulation using the infection vector")
    plt.legend()
    plt.show()
