from monty.json import MSONable
from typing import List
from .reaction_library import ReactionLibrary
from .scored_reaction import ScoredReaction
from .scored_reaction_set import ScoredReactionSet

class LibrarySet(MSONable):
    """Contains a set of ReactionLibrary objects."""

    def __init__(self, libraries: list[ReactionLibrary]):
        
        lib_lenghts = [len(lib.get_rxns_at_temp(t)) for lib in libraries for t in lib.temps ]
        if len(set(lib_lenghts)) != 1:
            raise ValueError("All libraries must have the same number of reactions at each temperature")
        
        temps = [t for lib in libraries for t in lib.temps]
        rxns_id_to_rxns = {i : rxn for i, rxn in enumerate(libraries[0].get_rxns_at_temp(temps[0]).reactions)}
        
        score_map = {}
        energy_per_atom_map = {}
        for l, lib in enumerate(libraries):
            score_map[l]= {t: {i: lib.get_rxns_at_temp(t).reactions[i].competitiveness for i in rxns_id_to_rxns.keys()} for t in temps}
            energy_per_atom_map[l]= {t: {i: lib.get_rxns_at_temp(t).reactions[i].energy_per_atom for i in rxns_id_to_rxns.keys()} for t in temps}
            
        self.score_map = score_map
        self.energy_per_atom_map = energy_per_atom_map
        self.rxns_id_to_rxns = rxns_id_to_rxns
        self.temps = temps
        self.phases = libraries[0].phases   
        
    def get_rxn_by_id(self, id: int) -> ScoredReaction:
        return self.rxns_id_to_rxns[id]
    
    def get_rxn_score(self, id: int, temp: int) -> float:
        return self.score_map[id][temp]
    
    def get_rxns_at_temp(self, temp: int) -> List[ScoredReaction]:
        return [self.rxns_id_to_rxns[i] for i in self.rxns_id_to_rxns.keys() if self.get_rxn_score(i, temp) is not None]
    
    def get_libraries(self) -> List[ReactionLibrary]:
        libraries = []
        for l in range(len(self.score_map)):
            libraries.append(ReactionLibrary(self.phases))
            for t in self.temps:
                scores = self.score_map[l][t]
                energy_per_atom = self.energy_per_atom_map[l][t]
                rxns = [self.rxns_id_to_rxns[i] for i in self.rxns_id_to_rxns.keys() if scores[i] is not None]
                rxns_rescored = [ScoredReaction(rxn._reactants, rxn._products, scores[i], energy_per_atom=energy_per_atom[i]) for i, rxn in enumerate(rxns)]
                libraries[l].add_rxns_at_temp(ScoredReactionSet(rxns_rescored, self.phases), t)
        return libraries