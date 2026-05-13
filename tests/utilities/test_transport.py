import numpy as np


def test_transport_db_loads_and_has_docs():
    from rxn_ca.utilities.transport import TransportDatabase

    db = TransportDatabase.from_chemsys("Ba-Ti-O")
    assert len(db.docs) > 0


def test_transport_db_nearest_temperature_lookup():
    from rxn_ca.utilities.transport import TransportDatabase

    db = TransportDatabase.from_chemsys("Ba-Ti-O")

    # BaTiO3 is present at 1000.0K in the data; 1001 should choose nearest (1000)
    doc_exact = db.get_doc("BaTiO3", 1000.0)
    doc_near = db.get_doc("BaTiO3", 1001.0)
    assert doc_exact.reduced_formula == "BaTiO3"
    assert doc_near.reduced_formula == "BaTiO3"
    assert float(doc_near.temperature) == float(doc_exact.temperature)


def test_transport_doc_from_py_oats_serialized_dict_minimal_fields():
    from rxn_ca.utilities.transport import TransportDoc

    d = {
        "reduced_formula": "BaTiO3",
        "temperature": 1000.0,
        "species": ["Ba", "Ti", "O"],
        "mapping": {"Ba": 0, "Ti": 1, "O": 2},
        "L_tensor": {
            "@module": "numpy",
            "@class": "array",
            "dtype": "float64",
            "data": [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]],
        },
        "L_tensor_self": {
            "@module": "numpy",
            "@class": "array",
            "dtype": "float64",
            "data": [[1.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 6.0]],
        },
        # extra fields from py-oats that rxn_ca doesn't need should be ignored
        "times": [0.0, 1.0],
        "correlation_functions": {"Ba-Ba": {"total": [0.0, 1.0]}},
    }

    doc = TransportDoc.from_py_oats_dict(d)
    assert doc.reduced_formula == "BaTiO3"
    assert doc.temperature == 1000.0
    assert doc.species == ["Ba", "Ti", "O"]
    assert doc.mapping["Ti"] == 1
    assert isinstance(doc.L_tensor, np.ndarray)
    assert isinstance(doc.L_tensor_self, np.ndarray)
    assert doc.L_tensor_dis is not None
    np.testing.assert_allclose(doc.L_tensor_dis, doc.L_tensor - doc.L_tensor_self)


def test_compute_fluxes_returns_species_dict():
    from rxn_ca.utilities.transport import TransportDatabase

    db = TransportDatabase.from_chemsys("Ba-Ti-O")
    mu = {"Ba": 0.1, "Ti": -0.05, "O": 0.02}
    fluxes = db.compute_fluxes("BaTiO3", 1000.0, mu)
    assert isinstance(fluxes, dict)
    assert set(fluxes.keys()) == {"Ba", "Ti", "O"}


def test_compute_fluxes_identity_fallback_for_unknown_formula():
    from rxn_ca.utilities.transport import TransportDatabase

    db = TransportDatabase.from_chemsys("Ba-Ti-O")
    mu = {"Ba": 0.1, "Ti": -0.05, "O": 0.02}
    # Unknown formula should return identity-scaled flux without raising
    fluxes = db.compute_fluxes("XYZ", 1000.0, mu)
    assert isinstance(fluxes, dict)
    assert set(fluxes.keys()) == {"Ba", "Ti", "O"}
