from esclab.simulate import Component, Model


class _Source(Component):
    out = Component.Output()

    def calculate(self):
        self.out.v = 1.0


class _Pipe(Component):
    inlet = Component.Input(0.0)
    outlet = Component.Output()

    def calculate(self):
        self.outlet.v = self.inlet.v


class _Sink(Component):
    inlet = Component.Input(0.0)

    def calculate(self):
        pass


def _build_linear_model():
    model = Model()
    model.src = _Source()
    model.pipe = _Pipe()
    model.sink = _Sink()

    model.connect(model.src.out, model.pipe.inlet)
    model.connect(model.pipe.outlet, model.sink.inlet)
    return model


def test_add_network_graph_without_gui_tab():
    model = _build_linear_model()

    result = model.add_network_graph(show_tab=False)

    assert result["view_created"] is False
    assert result["n_components"] == 3
    assert result["n_edges"] == 2
    assert result["exported_paths"] == ()


def test_add_network_graph_requires_path_for_export():
    model = _build_linear_model()

    try:
        model.add_network_graph(show_tab=False, save_png=True)
    except ValueError as exc:
        assert "path_base" in str(exc)
    else:
        raise AssertionError("Expected ValueError when save_png=True without path_base")
