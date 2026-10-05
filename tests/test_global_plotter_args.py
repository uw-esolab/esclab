import copy

import pytest

from esclab.plotting import OnlinePlotter
from esclab.simulate import Model


@pytest.fixture
def model(monkeypatch):
    # Exercise argument resolution and construction without opening a window.
    monkeypatch.setattr(OnlinePlotter, "_build_widgets", lambda self: None)
    instance = Model()
    instance.initialize()
    return instance


def test_all_options_are_inherited_by_regular_and_stats_plotters(model):
    options = {
        "y1lim": (0, 10), "y2lim": (-1, 1),
        "y1label": "Primary", "y2label": "Secondary",
        "nmax_points": 25, "update_every": 5, "tab_title": "Global",
        "show_live": False, "time_units": "hr", "plotter_size": (0.5, 0.6),
    }
    model.settings.global_plotter_args = options.copy()
    model.add_plotter(1.0, 2.0)
    model.add_plotter_simstats()

    for plotter in model._plotters:
        assert plotter._y1lim == (0, 10)
        assert plotter._y2lim == (-1, 1)
        assert plotter.y1label == "Primary"
        assert plotter.y2label == "Secondary"
        assert plotter.nmax_points == 25
        assert plotter.update_every == 5
        assert plotter._tab_title == "Global"
        assert plotter.show_live is False
        assert plotter.time_units == "hr"
        assert plotter._plotter_size == (0.5, 0.6)
    assert model._plotters[0].y1_items[0].v == 1.0
    assert model._plotters[0].y2_items[0].v == 2.0
    assert model._plotters[1].y1_items[0]._series is model.historian["timestep"]
    assert model._plotters[1].y2_items[0]._series is model.historian["iterations"]
    assert model.settings.global_plotter_args == options


@pytest.mark.parametrize("stats", [False, True])
def test_explicit_none_and_builtin_default_values_override_globals(model, stats):
    model.settings.global_plotter_args = {
        "y1lim": (0, 10), "y2lim": (-1, 1), "tab_title": "Global",
        "show_live": False, "time_units": "hr", "update_every": 5,
        "y1label": "Global label", "plotter_size": (0.5, 0.6),
    }
    options = dict(y1lim=None, y2lim=None, tab_title=None, show_live=True,
                   time_units="sec", update_every=1, y1label="",
                   plotter_size=(0.9, 0.9))
    if stats:
        model.add_plotter_simstats(**options)
    else:
        model.add_plotter(1.0, **options)
    plotter = model._plotters[-1]
    assert plotter._y1lim is None
    assert plotter._y2lim is None
    assert plotter._tab_title is None
    assert plotter.show_live is True
    assert plotter.time_units == "sec"
    assert plotter.update_every == 1
    assert plotter.y1label == ""
    assert plotter._plotter_size == (0.9, 0.9)


def test_missing_globals_preserve_regular_and_stats_defaults(model):
    model.add_plotter(1.0)
    model.add_plotter_simstats()
    regular, stats = model._plotters
    for plotter in (regular, stats):
        assert plotter.show_live is True
        assert plotter.time_units == "sec"
        assert plotter.nmax_points == 1000
        assert plotter.update_every == 1
        assert plotter._y1lim is None
        assert plotter._y2lim is None
        assert plotter._plotter_size == (0.9, 0.9)
    assert regular.y1label == regular.y2label == ""
    assert regular._tab_title is None
    assert stats.y1label == "Timestep"
    assert stats.y2label == "Iterations"
    assert stats._tab_title == "Simulation Stats"


def test_class_defaults_are_copied_and_model_settings_are_independent(monkeypatch):
    defaults = {"show_live": False, "y1lim": [0, 10]}
    monkeypatch.setattr(Model.Settings, "global_plotter_args", defaults)
    first, second = Model(), Model()
    first.settings.global_plotter_args["show_live"] = True
    first.settings.global_plotter_args["y1lim"][0] = -10
    assert second.settings.global_plotter_args == defaults
    assert defaults == {"show_live": False, "y1lim": [0, 10]}
    second.add_plotter(1.0)
    assert second._plotters[-1].show_live is False


def test_updates_apply_only_to_new_plotters_and_do_not_mutate_settings(model):
    model.settings.global_plotter_args = {"show_live": False, "time_units": "hr"}
    model.add_plotter(1.0, time_units="min")
    assert model.settings.global_plotter_args["time_units"] == "hr"
    model.settings.global_plotter_args["time_units"] = "sec"
    model.add_plotter(2.0)
    assert [p.time_units for p in model._plotters] == ["min", "sec"]


@pytest.mark.parametrize("stats", [False, True])
@pytest.mark.parametrize("options,error,match", [
    ({"y1": [1.0]}, ValueError, "cannot contain 'y1'"),
    ({"y2": [2.0]}, ValueError, "cannot contain 'y2'"),
    ({"typo": True}, TypeError, "Unsupported plotter argument 'typo'"),
    (None, TypeError, "must be a dict"),
    ([], TypeError, "must be a dict"),
])
def test_invalid_globals_are_rejected_without_adding_plotters(model, stats, options, error, match):
    model.settings.global_plotter_args = copy.deepcopy(options)
    with pytest.raises(error, match=match):
        if stats:
            model.add_plotter_simstats()
        else:
            model.add_plotter(1.0)
    assert model._plotters == []


@pytest.mark.parametrize("options,error,match", [
    ({"y1": [1.0]}, ValueError, "cannot contain 'y1'"),
    ({"y2": [2.0]}, ValueError, "cannot contain 'y2'"),
    ({"typo": True}, TypeError, "Unsupported plotter argument 'typo'"),
])
def test_stats_reject_invalid_explicit_options(model, options, error, match):
    with pytest.raises(error, match=match):
        model.add_plotter_simstats(**options)
    assert model._plotters == []


def test_existing_positional_arguments_override_globals(model):
    model.settings.global_plotter_args = {
        "show_live": True, "time_units": "hr", "tab_title": "Global",
        "nmax_points": 25, "update_every": 5,
    }
    model.add_plotter(1.0, 2.0, None, None, "", "", 1000, 1, None, False, "sec")
    model.add_plotter_simstats(None, 1000, 1, False, "sec")
    for plotter in model._plotters:
        assert plotter.show_live is False
        assert plotter.time_units == "sec"
        assert plotter._tab_title is None
        assert plotter.nmax_points == 1000
        assert plotter.update_every == 1
