"""Tests for workflow schema serialization."""

import pytest
from monty.json import jsanitize, MontyDecoder
import json

from rxn_ca.workflow.schemas import ReactionLibraryData, SimulationOutput


class TestReactionLibraryData:
    """Tests for ReactionLibraryData serialization."""

    def test_as_dict_contains_class_info(self):
        """Test that as_dict includes @module and @class for MSONable."""
        data = ReactionLibraryData(
            phase_set_dict={"phases": ["Fe", "Si"]},
            reaction_library_dict={"reactions": []},
            chemical_system="Fe-Si",
            temperatures=[300, 400, 500],
            phases_available=["Fe", "Si", "FeSi"],
            reaction_library_path="/path/to/lib.json",
        )

        d = data.as_dict()

        assert "@module" in d
        assert "@class" in d
        assert d["@class"] == "ReactionLibraryData"

    def test_round_trip_serialization(self):
        """Test that from_dict(as_dict()) returns equivalent object."""
        original = ReactionLibraryData(
            phase_set_dict={"phases": ["Fe", "Si"], "volumes": {"Fe": 1.0}},
            reaction_library_dict={"reactions": [{"id": 1}]},
            chemical_system="Fe-Si",
            temperatures=[300, 400, 500],
            phases_available=["Fe", "Si", "FeSi"],
            reaction_library_path="/path/to/lib.json",
        )

        reconstructed = ReactionLibraryData.from_dict(original.as_dict())

        assert reconstructed.phase_set_dict == original.phase_set_dict
        assert reconstructed.reaction_library_dict == original.reaction_library_dict
        assert reconstructed.chemical_system == original.chemical_system
        assert reconstructed.temperatures == original.temperatures
        assert reconstructed.phases_available == original.phases_available
        assert reconstructed.reaction_library_path == original.reaction_library_path

    def test_round_trip_with_none_path(self):
        """Test serialization when reaction_library_path is None."""
        original = ReactionLibraryData(
            phase_set_dict={},
            reaction_library_dict={},
            chemical_system="Fe-Si",
            temperatures=[300],
            phases_available=["Fe"],
            reaction_library_path=None,
        )

        reconstructed = ReactionLibraryData.from_dict(original.as_dict())

        assert reconstructed.reaction_library_path is None

    def test_round_trip_with_metadata(self):
        """Test that metadata survives serialization."""
        original = ReactionLibraryData(
            phase_set_dict={},
            reaction_library_dict={},
            chemical_system="Fe-Si",
            temperatures=[300],
            phases_available=["Fe"],
            reaction_library_path="/path/to/lib.json",
            metadata={"n_entries": 3, "thermo_types": ["R2SCAN"]},
        )

        reconstructed = ReactionLibraryData.from_dict(original.as_dict())

        assert reconstructed.metadata == {"n_entries": 3, "thermo_types": ["R2SCAN"]}

    def test_from_dict_without_metadata(self):
        """Test that a dict written before metadata existed still loads."""
        d = {
            "phase_set_dict": {},
            "reaction_library_dict": {},
            "chemical_system": "Fe-Si",
            "temperatures": [300],
            "phases_available": ["Fe"],
            "reaction_library_path": None,
        }

        assert ReactionLibraryData.from_dict(d).metadata is None

    def test_jsanitize_compatibility(self):
        """Test that schema works with monty's jsanitize (used by jobflow)."""
        data = ReactionLibraryData(
            phase_set_dict={"test": "value"},
            reaction_library_dict={},
            chemical_system="Fe-Si",
            temperatures=[300, 400],
            phases_available=["Fe", "Si"],
        )

        # jsanitize should not raise
        sanitized = jsanitize(data, strict=True)

        assert isinstance(sanitized, dict)
        assert sanitized["chemical_system"] == "Fe-Si"

    def test_json_round_trip(self):
        """Test full JSON serialization round trip."""
        original = ReactionLibraryData(
            phase_set_dict={"phases": ["A", "B"]},
            reaction_library_dict={"rxns": [1, 2, 3]},
            chemical_system="A-B",
            temperatures=[300, 400, 500],
            phases_available=["A", "B", "AB"],
            reaction_library_path="/tmp/test.json",
        )

        # Serialize to JSON string
        json_str = json.dumps(original.as_dict())

        # Deserialize
        d = json.loads(json_str)
        reconstructed = ReactionLibraryData.from_dict(d)

        assert reconstructed.chemical_system == original.chemical_system
        assert reconstructed.temperatures == original.temperatures


class TestSimulationOutput:
    """Tests for SimulationOutput serialization."""

    def test_as_dict_contains_class_info(self):
        """Test that as_dict includes @module and @class for MSONable."""
        output = SimulationOutput(
            final_molar_amounts={"Fe": 0.5, "Si": 0.5},
            molar_amounts_trajectory={"Fe": [1.0, 0.5], "Si": [1.0, 0.5]},
            temperature_trajectory=[300, 400],
            step_indices=[0, 1],
            phase_set_dict={"phases": ["Fe", "Si"]},
        )

        d = output.as_dict()

        assert "@module" in d
        assert "@class" in d
        assert d["@class"] == "SimulationOutput"

    def test_round_trip_serialization(self):
        """Test that from_dict(as_dict()) returns equivalent object."""
        original = SimulationOutput(
            final_molar_amounts={"Fe": 0.3, "Si": 0.2, "FeSi": 0.5},
            molar_amounts_trajectory={
                "Fe": [1.0, 0.5, 0.3],
                "Si": [1.0, 0.5, 0.2],
                "FeSi": [0.0, 0.5, 0.5],
            },
            temperature_trajectory=[300, 400, 500],
            step_indices=[0, 10, 20],
            phase_set_dict={"phases": ["Fe", "Si", "FeSi"]},
            reaction_library_path="/path/to/lib.json",
            result_doc_path="/path/to/result.json",
            recipe_dict={"reactants": {"Fe": 1.0, "Si": 1.0}},
            chemical_system="Fe-Si",
            num_realizations=3,
            metadata={"experiment_id": "exp001"},
        )

        reconstructed = SimulationOutput.from_dict(original.as_dict())

        assert reconstructed.final_molar_amounts == original.final_molar_amounts
        assert reconstructed.molar_amounts_trajectory == original.molar_amounts_trajectory
        assert reconstructed.temperature_trajectory == original.temperature_trajectory
        assert reconstructed.step_indices == original.step_indices
        assert reconstructed.phase_set_dict == original.phase_set_dict
        assert reconstructed.reaction_library_path == original.reaction_library_path
        assert reconstructed.result_doc_path == original.result_doc_path
        assert reconstructed.recipe_dict == original.recipe_dict
        assert reconstructed.chemical_system == original.chemical_system
        assert reconstructed.num_realizations == original.num_realizations
        assert reconstructed.metadata == original.metadata

    def test_round_trip_with_defaults(self):
        """Test serialization with default/None values."""
        original = SimulationOutput(
            final_molar_amounts={},
            molar_amounts_trajectory={},
            temperature_trajectory=[],
            step_indices=[],
            phase_set_dict={},
        )

        reconstructed = SimulationOutput.from_dict(original.as_dict())

        assert reconstructed.reaction_library_path is None
        assert reconstructed.result_doc_path is None
        assert reconstructed.recipe_dict == {}
        assert reconstructed.chemical_system is None
        assert reconstructed.num_realizations == 1
        assert reconstructed.metadata is None

    def test_jsanitize_compatibility(self):
        """Test that schema works with monty's jsanitize (used by jobflow)."""
        output = SimulationOutput(
            final_molar_amounts={"Fe": 0.5},
            molar_amounts_trajectory={"Fe": [1.0, 0.5]},
            temperature_trajectory=[300, 400],
            step_indices=[0, 1],
            phase_set_dict={},
            chemical_system="Fe-Si",
            num_realizations=2,
        )

        # jsanitize should not raise
        sanitized = jsanitize(output, strict=True)

        assert isinstance(sanitized, dict)
        assert sanitized["chemical_system"] == "Fe-Si"
        assert sanitized["num_realizations"] == 2

    def test_json_round_trip(self):
        """Test full JSON serialization round trip."""
        original = SimulationOutput(
            final_molar_amounts={"A": 0.5, "B": 0.5},
            molar_amounts_trajectory={"A": [1.0, 0.5], "B": [0.0, 0.5]},
            temperature_trajectory=[300, 400],
            step_indices=[0, 10],
            phase_set_dict={"test": True},
            metadata={"key": "value"},
        )

        # Serialize to JSON string
        json_str = json.dumps(original.as_dict())

        # Deserialize
        d = json.loads(json_str)
        reconstructed = SimulationOutput.from_dict(d)

        assert reconstructed.final_molar_amounts == original.final_molar_amounts
        assert reconstructed.metadata == original.metadata
