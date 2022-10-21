import numpy as np

def normalize(scores):
    norm = np.linalg.norm(scores, ord=1)
    return scores / norm