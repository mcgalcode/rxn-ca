from pylattica.core import SimulationResult


class ReactionResult(SimulationResult):
    """A class that stores the result of running a simulation. Keeps track of all
    the steps that the simulation proceeded through, and the set of reactions that
    was used in the simulation.

    Compression (live frame storage) is configured by the runner via
    SimulationResult.configure_compression, driven by the compress_freq /
    num_frames arguments passed to runner.run(). This class adds no compression
    behavior of its own.
    """

    # Inherits __init__ and from_dict from SimulationResult (from_dict uses
    # cls(), so it produces a ReactionResult).