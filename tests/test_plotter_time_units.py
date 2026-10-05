import numpy as np
import pytest
from pyqtgraph.Qt import QtCore

from esclab.plotting import OnlinePlotter
from esclab.simulate import Model


def _plotter(**kwargs):
    return OnlinePlotter(
        [1.0], [2.0], None, None, "Primary", "Secondary", 2, 1,
        **kwargs,
    )


@pytest.mark.parametrize("units,seconds", [
    ("nsec", 1e-9), ("msec", 1e-3), ("sec", 1.0), ("min", 60.0),
    ("hr", 3600.0), ("days", 86400.0), ("weeks", 604800.0),
    ("years", 365 * 86400.0),
])
def test_time_and_variable_step_widths_are_converted(units, seconds):
    plotter = _plotter(show_live=False, time_units=units)
    plotter.preallocate(n_steps=1, timestep=seconds / 2)
    plotter.log_step(seconds, 0.25, 0.5)
    plotter.log_step(3 * seconds, 0.5, 0.75, timestep=2 * seconds)

    np.testing.assert_allclose(plotter.x_data[:2], [1, 3])
    np.testing.assert_allclose(plotter._step_timesteps[:2], [0.5, 2])
    np.testing.assert_allclose(plotter.y1_data[0][:2], [1, 1])
    np.testing.assert_allclose(plotter.y2_data[0][:2], [2, 2])
    np.testing.assert_allclose(plotter.conv_data[:2], [0.25, 0.5])
    np.testing.assert_allclose(plotter.iter_data[:2], [0.5, 0.75])


@pytest.mark.parametrize("units", ["", "seconds", "invalid", None, 60, []])
def test_invalid_units_fail_before_creating_widgets(units, monkeypatch):
    def unexpected_widgets(self):
        pytest.fail("Invalid units should be rejected before creating widgets")

    monkeypatch.setattr(OnlinePlotter, "_build_widgets", unexpected_widgets)
    with pytest.raises(ValueError, match="Unsupported time_units"):
        _plotter(time_units=units)


def test_model_plotters_forward_units_and_default_to_seconds():
    model = Model()
    model.initialize()
    model.add_plotter([1.0], show_live=False, time_units="hr")
    model.add_plotter_simstats(show_live=False, time_units="min")
    model.add_plotter([1.0], show_live=False)
    model.add_plotter_simstats(show_live=False)

    assert [p.time_units for p in model._plotters] == ["hr", "min", "sec", "sec"]
    for plotter in model._plotters:
        plotter.log_step(3600, timestep=60)
    assert [p.x_data[0] for p in model._plotters] == [1, 60, 3600, 3600]
    assert model.time == 0


@pytest.mark.parametrize("show_live", [True, False])
def test_rendering_ranges_strips_and_tooltips_use_selected_units(show_live):
    plotter = _plotter(show_live=show_live, time_units="hr")
    try:
        plotter.preallocate(timestep=1800)
        for time in (3600, 7200, 10800):
            plotter.log_step(time, conv_fraction=0.25, iter_fraction=0.5)
        if not show_live:
            plotter._finalize()
        OnlinePlotter.app.processEvents()

        plotter.apply_font_size()
        assert plotter.ax1.getAxis("bottom").labelText == "Time (hr)"
        np.testing.assert_allclose(plotter.y1_lines[0].xData, [1, 2, 3])
        np.testing.assert_allclose(plotter.y2_lines[0].xData, [1, 2, 3])
        for bar in (plotter._conv_bar_item, plotter._iter_bar_item):
            np.testing.assert_allclose(bar.opts["x0"], [0.75, 1.75, 2.75])
            np.testing.assert_allclose(bar.opts["x1"], [1.25, 2.25, 3.25])

        # Follow shows the most recent two points; Fit All includes all three.
        start, end = plotter.ax1.viewRange()[0]
        assert 1 < start < 2 < 3 < end
        plotter._fit_to_full_data()
        start, end = plotter.ax1.viewRange()[0]
        assert start < 1 < 3 < end

        # Manual zoom must search the buffer in displayed units as well.
        plotter.ax1.setXRange(2.4, 2.6, padding=0)
        plotter._render_visible_window()
        np.testing.assert_allclose(plotter.y1_lines[0].xData, [2, 3])
        plotter._fit_to_full_data()
        plotter._plot_tooltip_enabled = True
        for axis in (plotter.ax1, plotter.ax_conv, plotter.ax_iter):
            plotter._tooltip_last_idx = -1
            scene_pos = axis.vb.mapViewToScene(QtCore.QPointF(2, 0.5))
            plotter._on_strip_mouse_move((scene_pos,))
            assert plotter._tip_label.text().startswith("time | 2 hr\n")
            assert plotter._tooltip_vline.value() == pytest.approx(2)
    finally:
        if hasattr(plotter, "win"):
            OnlinePlotter.tab_widget.removeTab(OnlinePlotter.tab_widget.indexOf(plotter.win))
            OnlinePlotter.instances.remove(plotter)
            plotter._tip_label.close()
            plotter.win.close()
