from ..core.recipe import ReactionRecipe

from ..reactions import ReactionLibrary
from ..phases import SolidPhaseSet
from ..computing.schemas.ca_result_schema import RxnCAResultDoc, compress_doc, get_metadata_from_results

from rxn_network.reactions.reaction_set import ReactionSet
from pylattica.core import Simulation

import multiprocessing as mp
import ray

from .single_sim import run_single_sim
from .get_scored_rxns import get_scored_rxns

_reaction_lib = "reaction_lib"
_recipe = "recipe"
_initial_simulation = "initial_simulation"
_compress_freq = "compress_freq"
_num_frames = "num_frames"
_analysis_only = "analysis_only"


def _get_result(_):

    result: RxnCAResultDoc = run_single_sim(
        mp_globals[_recipe],
        reaction_lib=mp_globals.get(_reaction_lib),
        initial_simulation=mp_globals.get(_initial_simulation),
        compress_freq=mp_globals.get(_compress_freq, None),
        num_frames=mp_globals.get(_num_frames, None),
        analysis_only=mp_globals.get(_analysis_only, False),
    )
    # In analysis-only mode the per-realization payload is the ObservedResult;
    # otherwise it is the full ReactionResult.
    if result.observed_results:
        return (result.observed_results[0], result.final_simulation)
    return (result.results[0], result.final_simulation)


def run_sim_parallel(recipe: ReactionRecipe,
                     base_reactions: ReactionSet = None,
                     reaction_lib: ReactionLibrary = None,
                     initial_simulation: Simulation = None,
                     phase_set: SolidPhaseSet = None,
                     existing_lib: ReactionLibrary = None,
                     compress_freq: int = None,
                     num_frames: int = None,
                     analysis_only: bool = False):
    """Run simulation with multiple realizations in parallel.

    Args:
        recipe: The reaction recipe to simulate
        base_reactions: Base reactions for scoring (required if reaction_lib not provided)
        reaction_lib: Pre-computed complete reaction library (skips all scoring)
        initial_simulation: Initial simulation state (optional)
        phase_set: Phase set for the system
        existing_lib: Existing library with some temps already scored.
            New temps will be added to this library incrementally.
        compress_freq: If set, store a full state snapshot every compress_freq
            steps instead of per-step diffs. Cannot be combined with num_frames.
        num_frames: If set, store roughly this many full state snapshots over
            the run. Cannot be combined with compress_freq.
        analysis_only: If True, retain only the reduced per-phase-volume view of
            each realization (via a PhaseVolumeObserver) at the compress_freq /
            num_frames cadence, stored under RxnCAResultDoc.observed_results.

    Returns:
        RxnCAResultDoc with averaged results from all realizations
    """
    if base_reactions is None and reaction_lib is None:
        raise ValueError("Must provide either base_reactions or reaction_lib")

    if reaction_lib is None:
        # Get required temperatures
        required_temps = recipe.heating_schedule.all_temps
        cached_temps = existing_lib.temps if existing_lib else []
        new_temps = [t for t in required_temps if t not in cached_temps]

        if new_temps:
            print("================= RETRIEVING AND SCORING REACTIONS =================")
            if existing_lib:
                print(f"    (Reusing {len(cached_temps)} cached temps, scoring {len(new_temps)} new temps)")

        reaction_lib: ReactionLibrary = get_scored_rxns(
            base_reactions,
            heating_sched=recipe.heating_schedule,
            scorer_class=recipe.get_score_class(),
            phase_set=phase_set,
            existing_lib=existing_lib,
        )

    print()
    print()
    print()

    print(f'================= RUNNING SIMULATION w/ {recipe.num_realizations} REALIZATIONS =================')


    global mp_globals

    mp_globals = {
        _reaction_lib: reaction_lib,
        _recipe: recipe,
        _initial_simulation: initial_simulation,
        _compress_freq: compress_freq,
        _num_frames: num_frames,
        _analysis_only: analysis_only,
    }

    with mp.get_context("fork").Pool(recipe.num_realizations) as pool:
        sim_results = pool.map(_get_result, [_ for _ in range(recipe.num_realizations)])
        final_simulations = [res[1] for res in sim_results]
        results = [res[0] for res in sim_results]

    good_results = [res for res in results if res is not None]
    print(f'{len(good_results)} results achieved out of {len(sim_results)}')

    # In analysis-only mode each realization returns an ObservedResult, which
    # belongs in observed_results rather than results.
    if analysis_only:
        result_doc = RxnCAResultDoc(
            recipe=recipe,
            observed_results=good_results,
            reaction_library=reaction_lib,
            phases=reaction_lib.phases,
            final_simulation=final_simulations[-1]
        )
    else:
        result_doc = RxnCAResultDoc(
            recipe=recipe,
            results=good_results,
            reaction_library=reaction_lib,
            phases=reaction_lib.phases,
            final_simulation=final_simulations[-1]
        )

    return result_doc


@ray.remote
def _run_single_realization(
    recipe: ReactionRecipe,
    recipe_index: int,
    realization_id: int,
    reaction_lib: ReactionLibrary,
    initial_sim: Simulation = None,
    compress: bool = True,
    num_steps: int = 500,
    **sim_kwargs
) -> tuple[int, int, any]:
    """
    Run a single realization of a recipe using Ray.
    
    Args:
        recipe: The ReactionRecipe to run
        recipe_index: Index of the recipe in the original list (for grouping results)
        realization_id: ID of this realization (0 to num_realizations-1)
        reaction_lib: The ReactionLibrary to use
        initial_sim: Optional initial Simulation
        **sim_kwargs: Additional kwargs to pass to run_single_sim
    
    Returns:
        Tuple of (recipe_index, realization_id, result, final_simulation)
    """
    result_doc: RxnCAResultDoc = run_single_sim(
        recipe, 
        reaction_lib=reaction_lib,
        initial_simulation=initial_sim,
        **sim_kwargs
    )
    
    result_doc.metadata = get_metadata_from_results(result_doc.results)
    
    if compress:
        result_doc = compress_doc(result_doc, num_steps)
    
    return (recipe_index, realization_id, result_doc.results[0], result_doc.final_simulation)


def run_multi_recipe_parallel_ray(
    recipes: list[ReactionRecipe],
    reaction_libraries: list[ReactionLibrary],
    initial_simulations: list[Simulation] = None,
    num_workers: int = None,
    memory_per_task_gb: float = None,
    cpus_per_task: int = 1,
    **kwargs
) -> list[RxnCAResultDoc]:
    """
    Run multiple recipes in parallel using Ray, with all realizations for each recipe
    distributed across available cores. 
    For each recipe[i], uses reaction_libraries[i] and optionally initial_simulations[i]
    to run all realizations for that recipe. 
    
    Args:
        recipes: List of ReactionRecipe objects
        reaction_libraries: List of ReactionLibrary objects (one per recipe)
        initial_simulations: Optional list of Simulation objects (one per recipe)
        num_workers: Number of CPU cores to use (default: all available on the current node)
        memory_per_task_gb: Memory requirement per task in GB (optional, for memory-aware scheduling)
        cpus_per_task: Number of CPUs per task (default: 1, since simulations typically don't benefit from multiple CPUs)
        **kwargs: Additional kwargs to pass to run_single_sim (e.g., phase_set)
    
    Returns:
        List of RxnCAResultDoc objects, one per recipe
    """
    MAX_CPUS_PER_TASK = 4 # Any more does not seem to improve performance
    
    # Validate inputs
    if len(reaction_libraries) != len(recipes):
        raise ValueError(
            f"Number of reaction_libraries ({len(reaction_libraries)}) must match "
            f"number of recipes ({len(recipes)})"
        )
    
    if initial_simulations is not None and len(initial_simulations) != len(recipes):
        raise ValueError(
            f"Number of initial_simulations ({len(initial_simulations)}) must match "
            f"number of recipes ({len(recipes)})"
        )
    
    # Initialize Ray if not already initialized
    if not ray.is_initialized():
        init_kwargs = {}
        if num_workers:
            init_kwargs['num_cpus'] = num_workers
        ray.init(**init_kwargs)
    
    # Calculate total realizations
    total_realizations = sum(recipe.num_realizations for recipe in recipes)
    
    # Determine CPUs per task
    # Default to 1 CPU per task since simulations typically don't benefit from multiple CPUs
    # and the overhead can actually slow things down
    cpus_per_task = max(1, min(cpus_per_task, MAX_CPUS_PER_TASK))
    
    available_cpus = int(ray.cluster_resources().get('CPU', 1))
    print(f"Using {cpus_per_task} CPU(s) per task (total: {total_realizations} realizations, {available_cpus} CPUs available)")
    
    # Generate all tasks: for each recipe, create tasks for all its realizations
    task_refs = []
    
    for recipe_idx, recipe in enumerate(recipes):
        num_realizations = recipe.num_realizations
        
        task_options = {"num_cpus": cpus_per_task}
        if memory_per_task_gb:
            task_options["memory"] = int(memory_per_task_gb * 1_000_000_000)  # Convert GB to bytes
        
        run_task = _run_single_realization.options(**task_options)
        
        # Launch all realizations for this recipe
        for realization_id in range(num_realizations):
            task_ref = run_task.remote(
                recipe=recipe,
                recipe_index=recipe_idx,
                realization_id=realization_id,
                reaction_lib=reaction_libraries[recipe_idx],
                initial_sim=initial_simulations[recipe_idx] if initial_simulations else None,
                **kwargs
            )
            task_refs.append(task_ref)
    
    print(f"================= RUNNING {total_realizations} REALIZATIONS ACROSS {len(recipes)} RECIPES ON {ray.cluster_resources()['CPU']} CORES =================")
    print(f"Distributing tasks across available cores...")
    
    # Execute all tasks (Ray handles scheduling and load balancing)
    results = ray.get(task_refs)
    final_simulations = [res[3] for res in results]
    
    # Group results by recipe
    recipe_results = {}
    for recipe_idx, _, result, _ in results:
        if recipe_idx not in recipe_results:
            recipe_results[recipe_idx] = []
        recipe_results[recipe_idx].append(result)
    
    result_docs = []
    for recipe_idx, recipe in enumerate(recipes):
        good_results = [res for res in recipe_results.get(recipe_idx, []) if res is not None]
        num_good = len(good_results)
        num_total = recipe.num_realizations
        
        print(f"Recipe {recipe_idx}: {num_good}/{num_total} successful realizations")
        
        result_doc = RxnCAResultDoc(
            recipe=recipe,
            results=good_results,
            reaction_library=reaction_libraries[recipe_idx],
            phases=reaction_libraries[recipe_idx].phases,
            final_simulation=final_simulations[recipe_idx]
        )
        result_docs.append(result_doc)
    
    return result_docs