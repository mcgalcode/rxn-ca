#!/usr/bin/env python
"""Command-line interface entry points for rxn-ca."""

import argparse
import json
import os
import sys


def react():
    """Run the rxn-ca automaton simulation."""
    from rxn_ca.computing.schemas.ca_result_schema import get_metadata_from_results
    from rxn_ca.core.recipe import ReactionRecipe
    from rxn_ca.reactions import ReactionLibrary
    from rxn_ca.utilities.parallel_sim import run_sim_parallel
    from rxn_ca.utilities.prints import print_banner
    from rxn_ca.utilities.single_sim import run_single_sim
    from pylattica.core import Simulation

    parser = argparse.ArgumentParser(
        prog="react",
        description="Runs the rxn-ca automaton",
    )

    parser.add_argument('recipe_location')
    parser.add_argument('-d', '--input-dir')
    parser.add_argument('-o', '--output-file')
    parser.add_argument('-p', '--output-dir')
    parser.add_argument('--no-compress', action='store_true',
                        help='Disable live compression (keeps all frames)')
    parser.add_argument('--num-frames', type=int, default=None,
                        help='Store roughly this many snapshots over the run '
                             '(sets the frame/observation cadence). Mutually '
                             'exclusive with the default compress_freq.')
    parser.add_argument('--analysis-only', action='store_true',
                        help='Retain only the reduced per-phase-volume view of '
                             'each frame (no full states). Much smaller on disk, '
                             'but only volume-derived analyses are available.')
    parser.add_argument('--count-reactions', action='store_true',
                        help='Assemble reaction choice metadata (requires --no-compress)')
    parser.add_argument('-l', '--reaction-library-file')
    parser.add_argument('-i', '--initial-simulation-file')

    parser.add_argument('-s', '--single', default=False, action='store_true')
    parser.add_argument('--store-lib', default=False, action=argparse.BooleanOptionalAction)

    args = parser.parse_args()

    if args.count_reactions and not args.no_compress:
        print("Error: --count-reactions requires --no-compress.")
        print("       Live compression (enabled by default) only stores frames at intervals,")
        print("       losing most REACTION_CHOSEN values needed for --count-reactions.")
        sys.exit(1)

    if args.analysis_only and args.count_reactions:
        print("Error: --analysis-only cannot be combined with --count-reactions.")
        print("       Reaction-choice metadata requires full per-step states, which")
        print("       analysis-only mode does not retain.")
        sys.exit(1)

    if args.analysis_only and args.no_compress:
        print("Error: --analysis-only cannot be combined with --no-compress.")
        print("       --no-compress retains every full frame; --analysis-only retains none.")
        sys.exit(1)

    output_file_arg = args.output_file
    output_dir = args.output_dir
    recipe_location = args.recipe_location
    reaction_library_filename = args.reaction_library_file
    initial_simulation_filename = args.initial_simulation_file
    store_lib = args.store_lib

    print_banner()

    print(recipe_location)

    if reaction_library_filename is not None:
        print(f"Reading reaction library from {reaction_library_filename}...")
        rxn_lib = ReactionLibrary.from_file(reaction_library_filename)
        reaction_set = None
        phases = rxn_lib.phases
    else:
        print("--reaction-library-file is required to use the react script")
        sys.exit(1)

    if initial_simulation_filename is not None:
        initial_simulation = Simulation.from_file(initial_simulation_filename)
    else:
        initial_simulation = None

    if os.path.isdir(recipe_location):
        recipe_filenames = [os.path.join(recipe_location, fpath) for fpath in os.listdir(recipe_location)]
        recipe_filenames = [fpath for fpath in recipe_filenames if os.path.isfile(fpath)]
    else:
        recipe_filenames = [recipe_location]

    print(f"Identified the following recipes: {', '.join(recipe_filenames)}")

    for recipe_filename in recipe_filenames:
        print(f"Reading recipe from {recipe_filename}...")
        recipe = ReactionRecipe.from_file(recipe_filename)

        if output_file_arg is None:
            if recipe.name is None:
                output_fname = recipe_filename.split("/")[-1]
            else:
                output_fname = f"{recipe.name}.json"

            if output_dir is None:
                output_dir = os.path.dirname(recipe_filename)

            output_file = os.path.join(output_dir, output_fname)
        else:
            output_file = output_file_arg

        print(f"Choosing {output_file} as output location")

        # Live compression is enabled by default; --no-compress disables it
        use_live_compress = not args.no_compress

        # Determine the frame / observation cadence. An explicit --num-frames
        # wins; otherwise a default compress_freq is used whenever frames are
        # being retained (live compression or analysis-only observation).
        num_frames = args.num_frames
        compress_freq = None
        if num_frames is None and (use_live_compress or args.analysis_only):
            compress_freq = 500

        # When not embedding the library, reference it by its file path so it is
        # omitted from the artifact but still lazily recoverable.
        reaction_lib_path = None if store_lib else reaction_library_filename

        runner = run_single_sim if args.single else run_sim_parallel

        result_doc = runner(
            recipe,
            base_reactions=reaction_set,
            reaction_lib=rxn_lib,
            initial_simulation=initial_simulation,
            phase_set=phases,
            compress_freq=compress_freq,
            num_frames=num_frames,
            analysis_only=args.analysis_only,
            reaction_lib_path=reaction_lib_path,
        )

        if args.count_reactions:
            print("Assembling metadata from results...")
            result_doc.metadata = get_metadata_from_results(result_doc.results)

        print(f'================= SAVING RESULTS to {output_file} =================')

        if store_lib:
            print("Embedding reaction library in result doc...")
        else:
            print(f"Referencing reaction library by path: {reaction_library_filename}")

        if args.analysis_only:
            analysis_fpath = output_file.split(".")[0] + "_analysis.json"
            print(f"Saving analysis-only results to {analysis_fpath}")
            result_doc.to_file(analysis_fpath)
        elif use_live_compress:
            # With live_compress, result is already compressed - just save it
            compressed_fpath = output_file.split(".")[0] + "_compressed.json"
            print(f"Saving compressed results to {compressed_fpath}")
            result_doc.to_file(compressed_fpath)
        else:
            print(f"Saving full results to {output_file}...")
            result_doc.to_file(output_file)


def enumerate_rxns():
    """Enumerate reactions for a chemical system."""
    from rxn_ca.phases import SolidPhaseSet
    from rxn_ca.utilities.enumerate_rxns import enumerate_rxns
    from rxn_ca.utilities.get_entries import get_entries
    from rxn_network.entries.entry_set import GibbsEntrySet

    parser = argparse.ArgumentParser(
        prog="enumerate",
        description="Enumerates reactions for a chemical system",
    )

    parser.add_argument('-m', '--reaction-manifest')

    parser.add_argument('-s', '--chemical-system')
    parser.add_argument('-e', '--energy-cutoff')
    parser.add_argument('-f', '--formulas-to-include')
    parser.add_argument('-t', '--extra-entry-set-file')

    parser.add_argument('-o', '--output-file')

    args = parser.parse_args()

    manifest_filename = args.reaction_manifest
    extra_entry_set = None
    entry_metadata = None

    if manifest_filename is not None:
        with open(manifest_filename, 'r+') as f:
            manifest = json.load(f)

        chem_sys = manifest.get("chemical_system")
        formulas_to_include = manifest.get("formulas_to_include", [])
        energy_cutoff = manifest.get("energy_cutoff", 0.01)
    else:
        chem_sys = args.chemical_system
        energy_cutoff = float(args.energy_cutoff) if args.energy_cutoff is not None else 0.01
        formulas_to_include = args.formulas_to_include
        if formulas_to_include is None:
            formulas_to_include = []
        else:
            formulas_to_include = formulas_to_include.split(",")

    if chem_sys is None:
        print("chemical system must be provided! either by the -s flag, or as the chemical_system value in the reaction manifest")
        sys.exit(1)

    print(f"Enumerating rxns for {chem_sys} using energy cutoff {energy_cutoff} and ensuring formulas {formulas_to_include} are present")

    if args.extra_entry_set_file is not None:
        with open(args.extra_entry_set_file, "r+") as f:
            eset_dict = json.load(f)
            extra_entry_set = GibbsEntrySet.from_dict(eset_dict['entry_set'])
            entry_metadata = eset_dict['entry_metadata']

    custom_entries = extra_entry_set.entries_list if extra_entry_set is not None else None

    entries = get_entries(
        chem_sys,
        stability_cutoff=energy_cutoff,
        formulas_to_include=formulas_to_include,
        custom_entries=custom_entries
    )

    phase_set = SolidPhaseSet.from_entry_set(entries, entry_metadata=entry_metadata)

    result = enumerate_rxns(
        chem_sys=chem_sys,
        stability_cutoff=energy_cutoff,
        formulas_to_include=formulas_to_include
    )

    output_filename = args.output_file

    if output_filename is None:
        output_filename = f'{chem_sys}_reactions.json'

    result.to_file(output_filename)
    print(f"Saved enumerated reactions to {output_filename}")


def suggest_precursors():
    """Suggest precursor combinations for synthesizing a target phase."""
    from rxn_ca.optimization.precursor_selection import (
        get_expanded_elements,
        suggest_recipes,
        suggest_recipes_from_literature,
    )
    from rxn_ca.optimization.synthesis_data import SynthesisDataset
    from rxn_ca.utilities.get_entries import get_entries

    parser = argparse.ArgumentParser(
        prog="suggest-precursors",
        description="Suggest precursor combinations for a target phase",
    )

    parser.add_argument(
        'target',
        help='Target phase formula (e.g., BaTiO3)',
    )
    parser.add_argument(
        '-n', '--n-precursors',
        type=int,
        default=2,
        help='Number of precursors per template (default: 2)',
    )
    parser.add_argument(
        '-m', '--max-templates',
        type=int,
        default=10,
        help='Maximum number of templates to show (default: 10)',
    )
    parser.add_argument(
        '-e', '--energy-cutoff',
        type=float,
        default=0.1,
        help='Metastability cutoff for Materials Project query (default: 0.1 eV/atom)',
    )
    parser.add_argument(
        '-l', '--literature-data',
        type=str,
        default=None,
        help='Path to synthesis dataset JSON for literature-based ranking',
    )
    parser.add_argument(
        '--min-frequency',
        type=int,
        default=5,
        help='Minimum literature occurrences when using --literature-data (default: 5)',
    )
    parser.add_argument(
        '--json',
        action='store_true',
        help='Output as JSON instead of human-readable format',
    )

    args = parser.parse_args()

    target = args.target
    print(f"Finding precursor combinations for: {target}", file=sys.stderr)

    # Expand elements to include precursor anions (C, N, H, etc.)
    print("Expanding element set for precursor phases...", file=sys.stderr)
    elements = get_expanded_elements(target)
    print(f"  Elements: {', '.join(sorted(elements))}", file=sys.stderr)

    # Fetch entries from Materials Project
    print(f"Fetching phases from Materials Project (cutoff={args.energy_cutoff} eV/atom)...", file=sys.stderr)
    entries = get_entries(
        elements,
        metastability_cutoff=args.energy_cutoff,
    )
    available_phases = [e.composition.reduced_formula for e in entries]
    print(f"  Found {len(available_phases)} phases", file=sys.stderr)

    # Generate templates
    if args.literature_data:
        print(f"Loading literature data from {args.literature_data}...", file=sys.stderr)
        dataset = SynthesisDataset.from_json_file(args.literature_data)
        print(f"  Loaded {len(dataset.records)} synthesis records", file=sys.stderr)

        templates = suggest_recipes_from_literature(
            target_phase=target,
            available_phases=available_phases,
            synthesis_dataset=dataset,
            n_precursors=args.n_precursors,
            min_frequency=args.min_frequency,
            max_templates=args.max_templates,
        )
        score_key = "literature_score"
    else:
        templates = suggest_recipes(
            target_phase=target,
            available_phases=available_phases,
            n_precursors=args.n_precursors,
            practical_only=True,
            max_templates=args.max_templates,
        )
        score_key = "practicality_score"

    # Output results
    if args.json:
        output = []
        for t in templates:
            output.append({
                "precursors": t.precursors,
                "target": t.target_phase,
                "score": t.metadata.get(score_key, 0),
                "metadata": t.metadata,
            })
        print(json.dumps(output, indent=2))
    else:
        if not templates:
            print(f"\nNo valid precursor combinations found for {target}")
            sys.exit(0)

        print(f"\nPrecursor combinations for {target}:")
        print("-" * 60)
        for i, t in enumerate(templates, 1):
            precursor_str = " + ".join(t.precursors)
            score = t.metadata.get(score_key, 0)
            print(f"{i:2d}. {precursor_str}")
            print(f"    Score: {score:.2f} ({score_key.replace('_', ' ')})")
        print("-" * 60)
        print(f"Total: {len(templates)} combinations")


def build_library():
    """Build a reaction library from a reaction set."""
    from rxn_ca.computing.schemas.enumerated_rxns_schema import EnumeratedRxnsModel
    from rxn_ca.core.recipe import ReactionRecipe
    from rxn_ca.phases import DEFAULT_GASES, SolidPhaseSet
    from rxn_ca.utilities.get_scored_rxns import get_scored_rxns

    parser = argparse.ArgumentParser(
        prog="build-library",
        description="Build reaction library from reaction set",
    )

    parser.add_argument('-e', '--reaction-enumeration-file')
    parser.add_argument('-r', '--recipe-file')

    parser.add_argument('-o', '--output-file')

    args = parser.parse_args()

    reaction_enumeration_file = args.reaction_enumeration_file
    recipe_file = args.recipe_file

    recipe: ReactionRecipe = ReactionRecipe.from_file(recipe_file)
    enumeration: EnumeratedRxnsModel = EnumeratedRxnsModel.from_file(reaction_enumeration_file)

    gases = [*recipe.additional_gas_phases, *DEFAULT_GASES]

    print("Building phase set using gases ", gases)

    phase_set = SolidPhaseSet.from_rxn_set(enumeration.rxn_set, gas_phases=gases)

    lib = get_scored_rxns(
        enumeration.rxn_set,
        recipe.heating_schedule,
        exclude_theoretical=recipe.exclude_theoretical,
        exclude_phases=recipe.exclude_phases,
        exclude_pure_elements=recipe.exclude_pure_elements,
        phase_set=phase_set
    )

    output_filename = args.output_file

    if output_filename is None:
        output_filename = 'reaction_library.json'

    lib.to_file(output_filename)
    print(f"Saved reaction library to {output_filename}")
