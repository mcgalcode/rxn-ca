"""Base classes and data structures for optimization module."""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional, Union
import json
import time

import pandas as pd


class ParameterType(str, Enum):
    """Types of optimization parameters."""
    CONTINUOUS = "continuous"
    DISCRETE = "discrete"
    CATEGORICAL = "categorical"
    PRECURSOR_SLOT = "precursor_slot"


@dataclass
class Parameter:
    """Base class for optimization parameters."""
    name: str
    param_type: ParameterType

    def validate(self, value: Any) -> bool:
        """Validate that a value is valid for this parameter."""
        raise NotImplementedError


@dataclass
class ContinuousParameter(Parameter):
    """A continuous parameter with bounds."""
    low: float
    high: float
    param_type: ParameterType = field(default=ParameterType.CONTINUOUS, init=False)

    def __post_init__(self):
        if self.low >= self.high:
            raise ValueError(f"low ({self.low}) must be less than high ({self.high})")

    def validate(self, value: float) -> bool:
        return self.low <= value <= self.high


@dataclass
class DiscreteParameter(Parameter):
    """A discrete parameter with step size."""
    low: float
    high: float
    step: float
    param_type: ParameterType = field(default=ParameterType.DISCRETE, init=False)

    def __post_init__(self):
        if self.low >= self.high:
            raise ValueError(f"low ({self.low}) must be less than high ({self.high})")
        if self.step <= 0:
            raise ValueError(f"step ({self.step}) must be positive")

    @property
    def values(self) -> List[float]:
        """Get all valid discrete values."""
        import numpy as np
        return list(np.arange(self.low, self.high + self.step / 2, self.step))

    def validate(self, value: float) -> bool:
        return value in self.values


@dataclass
class CategoricalParameter(Parameter):
    """A categorical parameter with discrete choices."""
    choices: List[str]
    param_type: ParameterType = field(default=ParameterType.CATEGORICAL, init=False)

    def __post_init__(self):
        if not self.choices:
            raise ValueError("choices must not be empty")

    def validate(self, value: str) -> bool:
        return value in self.choices


@dataclass
class PrecursorSlotParameter(Parameter):
    """A parameter for selecting precursors with chemical encodings.

    This is a special categorical parameter that uses chemical formulas
    as choices, enabling chemically-informed optimization via molecular
    fingerprints (e.g., MORDRED descriptors in BayBE).
    """
    candidates: List[str]  # Chemical formulas like "BaCO3", "BaO"
    param_type: ParameterType = field(default=ParameterType.PRECURSOR_SLOT, init=False)

    def __post_init__(self):
        if not self.candidates:
            raise ValueError("candidates must not be empty")

    @property
    def choices(self) -> List[str]:
        """Alias for candidates to match CategoricalParameter interface."""
        return self.candidates

    def validate(self, value: str) -> bool:
        return value in self.candidates


@dataclass
class OptimizationResult:
    """Result from a single optimization evaluation."""
    parameters: Dict[str, Any]
    score: float
    recipe: Any = None  # OptimizableRecipe
    result_doc: Any = None  # RxnCAResultDoc
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __repr__(self) -> str:
        param_str = ", ".join(f"{k}={v}" for k, v in self.parameters.items())
        return f"OptimizationResult(score={self.score:.4f}, {param_str})"


class OptimizationHistory:
    """History of optimization results with utility methods."""

    def __init__(self):
        self._results: List[OptimizationResult] = []

    def add(self, result: OptimizationResult) -> None:
        """Add a result to the history."""
        self._results.append(result)

    def add_batch(self, results: List[OptimizationResult]) -> None:
        """Add multiple results to the history."""
        self._results.extend(results)

    @property
    def results(self) -> List[OptimizationResult]:
        """Get all results."""
        return self._results

    @property
    def best_result(self) -> Optional[OptimizationResult]:
        """Get the result with the highest score."""
        if not self._results:
            return None
        return max(self._results, key=lambda r: r.score)

    @property
    def scores(self) -> List[float]:
        """Get all scores."""
        return [r.score for r in self._results]

    def __len__(self) -> int:
        return len(self._results)

    def __iter__(self):
        return iter(self._results)

    def __getitem__(self, idx: int) -> OptimizationResult:
        return self._results[idx]

    def to_dataframe(self) -> pd.DataFrame:
        """Convert history to a pandas DataFrame."""
        if not self._results:
            return pd.DataFrame()

        rows = []
        for i, result in enumerate(self._results):
            row = {"iteration": i, "score": result.score}
            row.update(result.parameters)
            rows.append(row)

        return pd.DataFrame(rows)

    def get_best_n(self, n: int) -> List[OptimizationResult]:
        """Get the top N results by score."""
        return sorted(self._results, key=lambda r: r.score, reverse=True)[:n]


class BaseOptimizer(ABC):
    """Abstract base class for optimizers.

    All optimizers implement a common interface with suggest/tell pattern
    and a convenience optimize() method for running the full loop.
    """

    def __init__(
        self,
        search_space: "SearchSpace",
        objective: "ObjectiveFunction",
        n_initial: int = 5,
        n_iterations: int = 20,
    ):
        """Initialize the optimizer.

        Args:
            search_space: The parameter space to search
            objective: The objective function to optimize
            n_initial: Number of initial random samples
            n_iterations: Number of optimization iterations after initial sampling
        """
        self.search_space = search_space
        self.objective = objective
        self.n_initial = n_initial
        self.n_iterations = n_iterations
        self.history = OptimizationHistory()

    @abstractmethod
    def suggest(self, n_suggestions: int = 1) -> List[Dict[str, Any]]:
        """Suggest parameter configurations to evaluate.

        Args:
            n_suggestions: Number of suggestions to return

        Returns:
            List of parameter dictionaries
        """
        pass

    @abstractmethod
    def tell(self, parameters: Dict[str, Any], score: float) -> None:
        """Report the result of an evaluation to the optimizer.

        Args:
            parameters: The evaluated parameter configuration
            score: The resulting score (higher is better)
        """
        pass

    def tell_batch(self, results: List[OptimizationResult]) -> None:
        """Report multiple results to the optimizer.

        Args:
            results: List of OptimizationResult objects
        """
        for result in results:
            self.tell(result.parameters, result.score)

    def optimize(
        self,
        verbose: bool = True,
        output_dir: Optional[str] = None,
    ) -> OptimizationHistory:
        """Run the full optimization loop.

        Args:
            verbose: Whether to print progress
            output_dir: Optional directory to save results. If provided, creates:
                - manifest.json: Optimization configuration and search space
                - results.json: Final optimization results
                - simulations/: Directory with individual simulation results

        Returns:
            OptimizationHistory containing all results
        """
        # Set up output directory if requested
        sim_output_dir = None
        if output_dir is not None:
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            # Create simulations subdirectory
            sim_output_dir = output_path / "simulations"
            sim_output_dir.mkdir(exist_ok=True)

            # Write manifest
            self._write_manifest(output_path)

        total_iterations = self.n_initial + self.n_iterations
        optimization_start = time.time()

        for i in range(total_iterations):
            if verbose:
                phase = "initial" if i < self.n_initial else "optimization"
                print(f"Iteration {i + 1}/{total_iterations} ({phase})")

            # Get suggestion
            suggestions = self.suggest(n_suggestions=1)
            params = suggestions[0]

            # Evaluate with timing
            iter_start = time.time()
            result = self.objective.evaluate(params)
            iter_duration = time.time() - iter_start

            # Store timing in result metadata
            result.metadata["duration_seconds"] = iter_duration
            result.metadata["iteration"] = i

            self.history.add(result)

            # Save simulation result if output_dir specified
            if sim_output_dir is not None:
                self._save_simulation_result(sim_output_dir, i, result)

            # Update optimizer
            self.tell(params, result.score)

            if verbose:
                print(f"  Score: {result.score:.4f} ({iter_duration:.1f}s)")
                if self.history.best_result:
                    print(f"  Best so far: {self.history.best_result.score:.4f}")

        self._total_duration = time.time() - optimization_start

        # Write final results
        if output_dir is not None:
            self._write_final_results(Path(output_dir))

        return self.history

    def _write_manifest(self, output_path: Path) -> None:
        """Write optimization manifest to output directory."""
        manifest = {
            "created_at": datetime.now().isoformat(),
            "optimizer_type": self.__class__.__name__,
            "n_initial": self.n_initial,
            "n_iterations": self.n_iterations,
            "total_iterations": self.n_initial + self.n_iterations,
            "search_space": self._serialize_search_space(),
            "objective_config": self._serialize_objective_config(),
        }

        manifest_path = output_path / "manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(manifest, f, indent=2)

    def _serialize_search_space(self) -> Dict[str, Any]:
        """Serialize search space for manifest."""
        params = []
        for p in self.search_space.parameters:
            param_dict = {"name": p.name, "type": p.param_type.value}
            if hasattr(p, "low"):
                param_dict["low"] = p.low
            if hasattr(p, "high"):
                param_dict["high"] = p.high
            if hasattr(p, "step"):
                param_dict["step"] = p.step
            if hasattr(p, "choices"):
                param_dict["choices"] = p.choices
            if hasattr(p, "candidates"):
                param_dict["candidates"] = p.candidates
            params.append(param_dict)
        return {"parameters": params}

    def _serialize_objective_config(self) -> Dict[str, Any]:
        """Serialize objective config for manifest."""
        config = self.objective.config
        return {
            "target_phase": config.target_phase,
            "scorer_type": config.scorer_type.value,
            "simulation_size": config.simulation_size,
            "num_realizations": config.num_realizations,
            "compress_freq": config.compress_freq,
            "num_frames": config.num_frames,
        }

    def _save_simulation_result(
        self, sim_dir: Path, iteration: int, result: "OptimizationResult"
    ) -> None:
        """Save a single simulation result and generate phase plot."""
        if result.result_doc is None:
            return

        # Save the result doc
        result_path = sim_dir / f"iteration_{iteration:03d}.json"
        result.result_doc.to_file(str(result_path))

        # Generate mass fraction plot
        try:
            self._generate_phase_plot(sim_dir, iteration, result)
        except Exception as e:
            print(f"Warning: Could not generate phase plot for iteration {iteration}: {e}")

    def _generate_phase_plot(
        self, sim_dir: Path, iteration: int, result: "OptimizationResult"
    ) -> None:
        """Generate mass fraction vs reaction coordinate plot."""
        from ..analysis.bulk_reaction_analyzer import BulkReactionAnalyzer
        from ..analysis.visualization.reaction_plotter import ReactionPlotter
        from ..analysis.visualization.phase_trace_calculator import PhaseTraceConfig

        # Create analyzer and plotter with 1% threshold to show all significant phases
        analyzer = BulkReactionAnalyzer.from_result_doc(result.result_doc)
        trace_config = PhaseTraceConfig(minimum_required_prevalence=0.01)
        plotter = ReactionPlotter(analyzer, trace_config=trace_config, include_heating_trace=True)

        # Generate mass fraction plot
        fig = plotter.plot_mass_fractions()

        # Add title with params
        params = result.parameters
        title_parts = [f"Iteration {iteration}"]
        if "hold_temp" in params:
            title_parts.append(f"{params['hold_temp']}K")
        title_parts.append(f"Score: {result.score:.4f}")
        fig.update_layout(title=" | ".join(title_parts))

        # Save as PNG
        plot_path = sim_dir / f"iteration_{iteration:03d}_phases.png"
        fig.write_image(str(plot_path), scale=2)

    def _write_final_results(self, output_path: Path) -> None:
        """Write final optimization results and generate plots."""
        # Compute timing stats
        durations = [r.metadata.get("duration_seconds", 0) for r in self.history.results]
        total_duration = getattr(self, "_total_duration", sum(durations))

        results_data = {
            "best_score": self.history.best_result.score if self.history.best_result else None,
            "best_params": self.history.best_result.parameters if self.history.best_result else None,
            "timing": {
                "total_seconds": total_duration,
                "mean_iteration_seconds": sum(durations) / len(durations) if durations else 0,
                "min_iteration_seconds": min(durations) if durations else 0,
                "max_iteration_seconds": max(durations) if durations else 0,
            },
            "all_results": [
                {
                    "score": r.score,
                    "params": r.parameters,
                    "duration_seconds": r.metadata.get("duration_seconds"),
                }
                for r in self.history.results
            ],
        }

        results_path = output_path / "results.json"
        with open(results_path, "w") as f:
            json.dump(results_data, f, indent=2)

        # Generate plots
        self._generate_plots(output_path, results_data["all_results"])

    def _generate_plots(self, output_path: Path, results: List[Dict[str, Any]]) -> None:
        """Generate optimization plots."""
        from .plotting import (
            plot_optimization_summary,
            plot_parameter_grid,
        )

        target_phase = self.objective.config.target_phase

        # Summary plot
        try:
            fig = plot_optimization_summary(results, target_phase=target_phase)
            fig.write_image(str(output_path / "summary.png"), scale=2)
        except Exception as e:
            print(f"Warning: Could not generate summary plot: {e}")

        # Parameter grid
        try:
            fig = plot_parameter_grid(results)
            fig.write_image(str(output_path / "parameter_exploration.png"), scale=2)
        except Exception as e:
            print(f"Warning: Could not generate parameter plot: {e}")

    def optimize_batch(
        self, batch_size: int = 1, verbose: bool = True
    ) -> OptimizationHistory:
        """Run optimization with batch evaluation.

        Args:
            batch_size: Number of parallel evaluations per iteration
            verbose: Whether to print progress

        Returns:
            OptimizationHistory containing all results
        """
        total_iterations = self.n_initial + self.n_iterations
        iterations = (total_iterations + batch_size - 1) // batch_size

        for i in range(iterations):
            remaining = total_iterations - len(self.history)
            current_batch = min(batch_size, remaining)

            if current_batch <= 0:
                break

            if verbose:
                print(f"Batch {i + 1}: evaluating {current_batch} configurations")

            # Get suggestions
            suggestions = self.suggest(n_suggestions=current_batch)

            # Evaluate batch
            results = self.objective.evaluate_batch(suggestions)
            self.history.add_batch(results)

            # Update optimizer
            self.tell_batch(results)

            if verbose:
                best_in_batch = max(results, key=lambda r: r.score)
                print(f"  Best in batch: {best_in_batch.score:.4f}")
                if self.history.best_result:
                    print(f"  Best overall: {self.history.best_result.score:.4f}")

        return self.history
