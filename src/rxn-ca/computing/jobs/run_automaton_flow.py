from .run_automaton import RunRxnAutomatonMaker

from jobflow import Flow

RUN_AUTOMATON_FLOW = "RUN_AUTOMATON_FLOW"

def run_automaton_flow(
                    chem_sys: str,
                    temp: int,
                    setup_style: str,
                    setup_args: dict,
                    num_steps: int,
                    scored_rxns_task_id: str = None,
                    db_connection_params: dict = {},
                    db_file: str = None,
                    dimensionality: int = 2,
                    inertia: float = 1,
                    open_species: dict = {},
                    free_species: list = [],
                    parallel: bool = True
    ):

    maker = RunRxnAutomatonMaker()
    job = maker.make(
        chem_sys = chem_sys,
        temp = temp,
        setup_style = setup_style,
        setup_args = setup_args,
        num_steps = num_steps,
        scored_rxns_task_id = scored_rxns_task_id,
        db_connection_params = db_connection_params,
        db_file = db_file,
        dimensionality = dimensionality,
        inertia = inertia,
        open_species = open_species,
        free_species = free_species,
        parallel = parallel,
    )

    return Flow([job], name=RUN_AUTOMATON_FLOW)
