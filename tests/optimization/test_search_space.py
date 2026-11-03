"""Tests for SearchSpace builder."""

import pytest
from rxn_ca.optimization import SearchSpace
from rxn_ca.optimization.base import (
    ContinuousParameter,
    DiscreteParameter,
    CategoricalParameter,
    PrecursorSlotParameter,
)


class TestSearchSpace:
    """Tests for the SearchSpace fluent builder API."""

    def test_empty_search_space(self):
        """Empty search space has no parameters."""
        ss = SearchSpace()
        assert len(ss.parameters) == 0

    def test_add_temperature_range(self):
        """Temperature range creates a discrete parameter."""
        ss = SearchSpace().add_temperature_range(800, 1200, step=100)

        assert len(ss.parameters) == 1
        param = ss.parameters[0]
        assert param.name == "hold_temp"
        assert isinstance(param, DiscreteParameter)
        assert param.low == 800
        assert param.high == 1200
        assert param.step == 100
        assert 900 in param.values
        assert 1000 in param.values

    def test_add_hold_time_range(self):
        """Hold time creates a continuous parameter."""
        ss = SearchSpace().add_hold_time_range(1, 10)

        assert len(ss.parameters) == 1
        param = ss.parameters[0]
        assert param.name == "hold_time"
        assert isinstance(param, DiscreteParameter)
        assert param.low == 1
        assert param.high == 10

    def test_add_ramp_step_time_range(self):
        """Ramp step time creates a continuous parameter."""
        ss = SearchSpace().add_ramp_step_time_range(5, 20)

        assert len(ss.parameters) == 1
        param = ss.parameters[0]
        assert param.name == "ramp_step_time"
        assert isinstance(param, DiscreteParameter)
        assert param.low == 5.0
        assert param.high == 20.0

    def test_add_precursor_slot(self):
        """Precursor slot creates a PrecursorSlotParameter."""
        ss = SearchSpace().add_precursor_slot("Ba_source", ["BaCO3", "BaO"])

        assert len(ss.parameters) == 1
        param = ss.parameters[0]
        assert param.name == "Ba_source"
        assert isinstance(param, PrecursorSlotParameter)
        assert "BaCO3" in param.candidates
        assert "BaO" in param.candidates

    def test_add_precursor_ratio(self):
        """Precursor ratio creates a continuous parameter."""
        ss = SearchSpace().add_precursor_ratio("Ba_source", 0.4, 0.6)

        assert len(ss.parameters) == 1
        param = ss.parameters[0]
        assert param.name == "Ba_source_ratio"
        assert isinstance(param, ContinuousParameter)
        assert param.low == 0.4
        assert param.high == 0.6

    def test_fluent_chaining(self):
        """Builder methods can be chained."""
        ss = (
            SearchSpace()
            .add_temperature_range(800, 1400, step=50)
            .add_hold_time_range(1, 15)
            .add_ramp_step_time_range(5, 20)
            .add_precursor_slot("Ba_source", ["BaCO3", "BaO"])
            .add_precursor_ratio("Ba_source", 0.4, 0.6)
        )

        assert len(ss.parameters) == 5
        param_names = [p.name for p in ss.parameters]
        assert "hold_temp" in param_names
        assert "hold_time" in param_names
        assert "ramp_step_time" in param_names
        assert "Ba_source" in param_names
        assert "Ba_source_ratio" in param_names

    def test_get_parameter_names(self):
        """Can get list of parameter names."""
        ss = (
            SearchSpace()
            .add_temperature_range(800, 1200)
            .add_hold_time_range(1, 10)
        )

        names = ss.parameter_names
        assert "hold_temp" in names
        assert "hold_time" in names

    def test_str_representation(self):
        """Search space has readable string representation."""
        ss = SearchSpace().add_temperature_range(800, 1200, step=100)

        s = str(ss)
        assert "hold_temp" in s
        assert "800" in s
        assert "1200" in s


class TestParameterValidation:
    """Tests for parameter validation."""

    def test_continuous_parameter_validation(self):
        """Continuous parameter validates bounds."""
        param = ContinuousParameter(name="test", low=0.0, high=1.0)

        assert param.validate(0.5) is True
        assert param.validate(0.0) is True
        assert param.validate(1.0) is True
        assert param.validate(-0.1) is False
        assert param.validate(1.1) is False

    def test_continuous_parameter_invalid_bounds(self):
        """Continuous parameter rejects low >= high."""
        with pytest.raises(ValueError):
            ContinuousParameter(name="test", low=1.0, high=0.0)

    def test_discrete_parameter_values(self):
        """Discrete parameter generates correct values."""
        param = DiscreteParameter(name="test", low=0, high=10, step=2)

        values = param.values
        assert 0 in values
        assert 2 in values
        assert 4 in values
        assert 10 in values
        assert 3 not in values

    def test_categorical_parameter_validation(self):
        """Categorical parameter validates choices."""
        param = CategoricalParameter(name="test", choices=["a", "b", "c"])

        assert param.validate("a") is True
        assert param.validate("d") is False

    def test_precursor_slot_parameter(self):
        """Precursor slot parameter works like categorical."""
        param = PrecursorSlotParameter(name="Ba_source", candidates=["BaCO3", "BaO"])

        assert param.validate("BaCO3") is True
        assert param.validate("BaO") is True
        assert param.validate("Ba") is False
        assert param.choices == ["BaCO3", "BaO"]
