"""Objective function for reaction optimization."""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Callable
import hashlib
import json

from .base import OptimizationResult
from .optimizable_recipe import OptimizableRecipe
from .utilities import (
    get_result_analysis,
    AnalyzedResult,
    MaximumProductScorer,
    FinalProductScorer,
)


class ScorerType(str, Enum):
    """Types of scoring functions for evaluating reaction outcomes."""
    FINAL = "final"  # Final amount of target phase
    MAXIMUM = "maximum"  # Peak amount of target phase during reaction


_SCORER_MAP = {
    ScorerType.FINAL: FinalProductScorer,
    ScorerType.MAXIMUM: MaximumProductScorer,
}


@dataclass
class ObjectiveConfig:
    """Configuration for the objective function.

    Args:
        target_phase: Chemical formula of the target product phase (e.g., "BaTiO3")
        scorer_type: How to score the reaction (final amount vs peak amount)
        simulation_size: Grid size for the CA simulation
        num_realizations: Number of simulation runs to average
        cache_results: Whether to cache results to avoid re-evaluation
        compress_freq: If set, store a full state snapshot every compress_freq
            steps instead of per-step diffs. This avoids slow O(n) reconstruction
            during analysis. Cannot be combined with num_frames.
        num_frames: If set, store roughly this many full state snapshots over the
            run. Cannot be combined with compress_freq.
    """
    target_phase: str
    scorer_type: ScorerType = ScorerType.FINAL
    simulation_size: int = 10
    num_realizations: int = 3
    cache_results: bool = True
    compress_freq: int = 50
    num_frames: Optional[int] = None


class ObjectiveFunction:
    """Wrapper for running simulations and scoring results.

    This class handles:
    - Converting optimizer parameters to OptimizableRecipe
    - Running the CA simulation
    - Analyzing and scoring the results
    - Caching to avoid redundant evaluations
    """

    def __init__(
        self,
        config: ObjectiveConfig,
        base_reactions: Any = None,  # ReactionSet
        reaction_lib: Any = None,  # ReactionLibrary
        phase_set: Any = None,  # SolidPhaseSet
        precursor_slot_map: Optional[Dict[str, str]] = None,
    ):
        """Initialize the objective function.

        Args:
            config: ObjectiveConfig with scoring settings
            base_reactions: ReactionSet for generating reaction library
            reaction_lib: Pre-computed ReactionLibrary (optional, saves computation)
            phase_set: SolidPhaseSet for the system
            precursor_slot_map: Optional fixed mapping of slot names to precursors.
                If not provided, slots are extracted from optimizer parameters.
        """
        self.config = config
        self.base_reactions = base_reactions
        self.reaction_lib = reaction_lib
        self.phase_set = phase_set
        self.precursor_slot_map = precursor_slot_map

        # Initialize scorer
        scorer_class = _SCORER_MAP[config.scorer_type]
        self.scorer = scorer_class(config.target_phase)

        # Cache for avoiding re-evaluation
        self._cache: Dict[str, OptimizationResult] = {}

        # Cumulative reaction library - reuses scored reactions across evaluations
        # This avoids re-scoring reactions at the same temperature multiple times
        self._cumulative_lib: Any = None  # ReactionLibrary, initialized lazily

    def _params_to_cache_key(self, params: Dict[str, Any]) -> str:
        """Generate a cache key for a parameter configuration."""
        # Sort and serialize params for consistent hashing
        sorted_params = sorted(params.items())
        param_str = json.dumps(sorted_params, sort_keys=True)
        return hashlib.md5(param_str.encode()).hexdigest()

    def _extract_precursor_slots(self, params: Dict[str, Any]) -> Dict[str, str]:
        """Extract precursor slot selections from parameters.

        Looks for parameters that are precursor selections (string values
        that aren't standard parameter names).
        """
        # Start with fixed mapping if provided
        slot_map = dict(self.precursor_slot_map) if self.precursor_slot_map else {}

        # Extract dynamic slots from params (these may override fixed slots)
        standard_params = {"hold_temp", "hold_time", "ramp_rate"}

        for key, value in params.items():
            if key in standard_params:
                continue
            if key.endswith("_ratio"):
                continue
            # This is likely a precursor slot selection
            if isinstance(value, str):
                slot_map[key] = value

        return slot_map

    def _params_to_recipe(self, params: Dict[str, Any]) -> OptimizableRecipe:
        """Convert optimizer parameters to an OptimizableRecipe."""
        precursor_slots = self._extract_precursor_slots(params)

        return OptimizableRecipe.from_params(
            params=params,
            precursor_slots=precursor_slots,
            simulation_size=self.config.simulation_size,
            num_simulations=self.config.num_realizations,
        )

    def _run_simulation(self, recipe: OptimizableRecipe) -> Any:
        """Run the CA simulation and return the result document.

        Uses cumulative library caching to avoid re-scoring reactions at
        temperatures that have already been scored in previous evaluations.
        """
        # Import here to avoid circular imports
        from ..utilities.parallel_sim import run_sim_parallel
        from ..utilities.single_sim import run_single_sim

        reaction_recipe = recipe.to_recipe()

        # Use pre-computed library if provided, otherwise use cumulative caching
        if self.reaction_lib is not None:
            existing_lib = self.reaction_lib
        else:
            existing_lib = self._cumulative_lib

        if recipe.num_simulations > 1:
            result_doc = run_sim_parallel(
                reaction_recipe,
                base_reactions=self.base_reactions,
                reaction_lib=self.reaction_lib,
                phase_set=self.phase_set,
                existing_lib=existing_lib,
                compress_freq=self.config.compress_freq,
                num_frames=self.config.num_frames,
            )
        else:
            result_doc = run_single_sim(
                reaction_recipe,
                base_reactions=self.base_reactions,
                reaction_lib=self.reaction_lib,
                phase_set=self.phase_set,
                existing_lib=existing_lib,
                compress_freq=self.config.compress_freq,
                num_frames=self.config.num_frames,
            )

        # Update cumulative library with newly scored temps
        if self.reaction_lib is None and result_doc is not None:
            self._cumulative_lib = result_doc.reaction_library

        return result_doc

    def _score_result(self, result_doc: Any) -> float:
        """Score a simulation result."""
        try:
            traces = get_result_analysis(result_doc)
            analyzed = AnalyzedResult(traces)
            return self.scorer.score(analyzed)
        except (KeyError, IndexError):
            # Target phase not present in results or empty trace
            return 0.0
        except Exception as e:
            print(f"Warning: Scoring failed with error: {e}")
            return 0.0

    def evaluate(self, params: Dict[str, Any]) -> OptimizationResult:
        """Evaluate a single parameter configuration.

        Args:
            params: Dictionary of parameter values

        Returns:
            OptimizationResult with score and metadata
        """
        # Check cache
        if self.config.cache_results:
            cache_key = self._params_to_cache_key(params)
            if cache_key in self._cache:
                return self._cache[cache_key]

        try:
            # Convert to recipe
            recipe = self._params_to_recipe(params)

            # Run simulation
            result_doc = self._run_simulation(recipe)

            # Score
            score = self._score_result(result_doc)

        except Exception as e:
            print(f"Warning: Evaluation failed with error: {e}")
            import traceback
            traceback.print_exc()
            score = 0.0
            recipe = None
            result_doc = None

        # Create result
        result = OptimizationResult(
            parameters=params,
            score=score,
            recipe=recipe,
            result_doc=result_doc,
        )

        # Cache
        if self.config.cache_results and recipe is not None:
            self._cache[cache_key] = result

        return result

    def evaluate_batch(
        self, params_list: List[Dict[str, Any]]
    ) -> List[OptimizationResult]:
        """Evaluate multiple parameter configurations.

        Runs sequentially but benefits from reaction library caching -
        temperatures scored in earlier evaluations are reused in later ones.

        Args:
            params_list: List of parameter dictionaries

        Returns:
            List of OptimizationResult objects
        """
        results = []
        for params in params_list:
            result = self.evaluate(params)
            results.append(result)
        return results

    def clear_cache(self) -> None:
        """Clear the result cache."""
        self._cache.clear()

    def clear_reaction_lib_cache(self) -> None:
        """Clear the cumulative reaction library cache."""
        self._cumulative_lib = None

    @property
    def cache_size(self) -> int:
        """Number of cached results."""
        return len(self._cache)

    @property
    def cached_temps(self) -> List[int]:
        """List of temperatures with cached scored reactions."""
        if self._cumulative_lib is None:
            return []
        return self._cumulative_lib.temps


class MockObjectiveFunction(ObjectiveFunction):
    """A mock objective function for testing without running simulations.

    Uses a provided callable to compute scores based on parameters.
    """

    def __init__(
        self,
        config: ObjectiveConfig,
        score_fn: Callable[[Dict[str, Any]], float],
    ):
        """Initialize the mock objective function.

        Args:
            config: ObjectiveConfig (only used for metadata)
            score_fn: Function that takes params dict and returns a score
        """
        self.config = config
        self.score_fn = score_fn
        self._cache: Dict[str, OptimizationResult] = {}

    def evaluate(self, params: Dict[str, Any]) -> OptimizationResult:
        """Evaluate using the mock scoring function."""
        # Check cache
        if self.config.cache_results:
            cache_key = self._params_to_cache_key(params)
            if cache_key in self._cache:
                return self._cache[cache_key]

        score = self.score_fn(params)

        result = OptimizationResult(
            parameters=params,
            score=score,
            recipe=None,
            result_doc=None,
        )

        # Cache
        if self.config.cache_results:
            self._cache[cache_key] = result

        return result
