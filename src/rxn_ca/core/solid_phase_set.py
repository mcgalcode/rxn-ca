from pylattica.discrete.phase_set import PhaseSet

class SolidPhaseSet(PhaseSet):

    FREE_SPACE = "Free Space"

    @classmethod
    def from_dict(cls, set_dict):
        return cls(
            set_dict["phases"],
            set_dict["volumes"]
        )

    def __init__(self, phases, volumes):
        phases = phases + [SolidPhaseSet.FREE_SPACE]
        self.volumes = volumes
        super().__init__(phases)

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "phases": self.phases,
            "volumes": self.volumes
        }

