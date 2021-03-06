from abc import ABC, abstractmethod
from rxn_ca.core.basic_simulation_result import BasicSimulationResult

from rxn_ca.core.neighborhoods import NeighborhoodView


class BasicController(ABC):

    def instantiate_result(self):
        pass

    @abstractmethod
    def get_new_state(self, nb_view: NeighborhoodView):
        pass

    def instantiate_result(self):
        return BasicSimulationResult()
