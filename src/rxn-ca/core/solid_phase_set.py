from pylattica.discrete.phase_set import PhaseSet

class SolidPhaseSet(PhaseSet):

    FREE_SPACE = "Free Space"

    def __init__(self, phases):
        phases = phases + [SolidPhaseSet.FREE_SPACE]
        super().__init__(phases)

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "phases": self.phases,
        }

