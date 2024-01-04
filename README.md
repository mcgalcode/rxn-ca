# rxn-ca - A lattice model for simulating solid state reactions

## Getting Started

### Enumerating Reactions

The first step in using `rxn-ca` is to enumerate reactions. This process utilizes the `reaction-network` package to calculate the stoichiometrically viable reactions and their energies using data from the Materials Project. The result of this step is an output file containing a set of reactions that will be used in subsequent steps.

To perform this enumeration, use the `enumerate` script, which is made available upon installation of `rxn-ca`. It accepts a single positional argument which is a filepath for a reaction manifest file. This file specifies the chemical system, and energy stability cutoffs, and any specific formulas that should be included regardless of their stability. Here is an example file:

```
{
    "chemical_system": "Ba-Ti-O",
    "formulas_to_include": ["Ba2Ti2O5"],
    "energy_cutoff": 0.02,
}
```

These parameters can be specified as command line arguments as well.

Once `rxn-ca` is installed, and your manifest file is written, enumeration is performed like this:

```
enumerate -m my_reaction_manifest.json -o my_output_file.json
```

If successful, the enumerated reactions will be written to the `my_output_file.json` file.

### Writing a reaction recipe

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
    heating_schedule=heating_schedule
)

recipe.to_file("BaS_Na2TiO3_recipe.json")
```

This snippet will create yet another JSON file containing the details of the recipe.

### Building a library of reactions for your recipe

`rxn-ca` takes into account the temperatures specified in the reaction recipe. Since the free energy change in each reaction is a function of temperature, the energetics of every reaction must be recalculated at every temperature specified in the recipe. When this is performed, the resulting collection of reactions and their energies is called a `ReactionLibrary`. You can build one using the `build-library` script.

All you need to specify for this script is the location of the file containing your reaction recipe, and the location of the reaction enumeration you performed earlier.

You can generate your reaction library file like this:

```
build-library -s my_reaction_set.json -r my_recipe.json -o my_output_library.json
```

The resulting `ReactionLibrary` will be saved in the `my_output_library.json` file.

### Running your reaction

Use the `react` script to run your reaction. In the simplest case, just pass the recipe and reaction library file names to the script:

```
react my_recipe.json -l my_reaction_library.json
```

The output of your reaction will be stored in a JSON file whose name will be printed at the end of the stdout.

### Analyzing your reaction outcome

We recommend analyzing your reaction in a jupyter notebook, where plotting capabilities are fully supported. All you will need is the filename of your reaction result produced by the `react` script.

The following snippet shows a basic analysis:

```
from rxn_ca.analysis.bulk_reaction_analyzer import BulkReactionAnalyzer

analyzer = BulkReactionAnalyzer.from_result_doc_file("my_result.json")
analyzer.plot_molar_phase_amounts()
```