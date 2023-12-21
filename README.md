# rxn-ca - A lattice model for simulating solid state reactions

## Getting Started

#### Setting up jobflow

`rxn-ca` uses `jobflow` to store reactions that are eventually used. Follow the jobflow setup. When you are finished, you should be able to retrieve a reference to your jobflow store in a python environment by running `JobflowSettings().JOB_STORE`.

### Enumerating Reactions

Use the `enumerate` script to enumerate reactions prior to running your simulation. It accepts a single positional argument which is a filepath for a reaction manifest file. This file specifies all the chemical systems, temperatures, and energy cutoffs for which reactions should be enumerated. You can see an example of this file in `data/rxn_manifest.json`.

Once `rxn-ca` is installed, and your manifest file is written, enumeration is performed like this:

```
enumerate my_reaction_manifest.json
```

### Writing your reaction recipe

Before running a reaction, you must specify the recipe that will run it. Here's an example python snippet that creates a simple recipe:

```
from rxn_ca.core.recipe import ReactionRecipe
from rxn_ca.core.heating import HeatingSchedule, HeatingStep

heating_schedule = HeatingSchedule(
    HeatingStep.sweep(600, 1100, stage_length=30000, step_size=100),
    HeatingStep.hold(1100, 60000),    
)

recipe = ReactionRecipe(
    chem_sys="Ba-Ti-O-S-Na",
    reactant_amounts={
        "BaS": 1,
        "Na2TiO3": 1,
    },
    particle_size=0.5,
    heating_schedule=heating_schedule,
    simulation_size=15,
    num_realizations=3
)

recipe.to_file("BaS_Na2TiO3_recipe.json")
```

### Running your reaction

Use the `react` script to run your reaction. In the simplest case, just pass the recipe file name to the script:

```
react my_recipe.json
```

The output of your reaction will be stored in a json file whose name will be printed at the end of the stdout.

### Analyzing your reaction outcome

We recommend analyzing your reaction in a jupyter notebook, where plotting capabilities are fully supported. All you will need is the filename of your reaction result produced by the `react` script.

The following snippet shows a basic analysis:

```
from rxn_ca.analysis.bulk_reaction_analyzer import BulkReactionAnalyzer

analyzer = BulkReactionAnalyzer.from_result_doc_file("my_result.json")
analyzer.plot_molar_phase_amounts()
```

## Debugging

### grpcio

On M1 Macs, importing code related to the reaction-network can cause breakages. The suggested action from the error message, reproduced below, should fix any issues.

```
Failed to import grpc on Apple Silicon. On Apple Silicon machines, try `pip uninstall grpcio; conda install grpcio`. Check out https://docs.ray.io/en/master/ray-overview/installation.html#m1-mac-apple-silicon-support for more details.
```
