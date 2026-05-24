import math

from pyqtgraph.Qt import QtCore, QtGui, QtWidgets

from esclab.plotting.online_plotter import OnlinePlotter


class NetworkTopologyView:
    def __init__(self, model, tab_title="Network", include_subnetworks=True):
        self.model = model
        self.include_subnetworks = include_subnetworks
        OnlinePlotter.ensure_window()

        self.view = QtWidgets.QGraphicsView()
        self.scene = QtWidgets.QGraphicsScene(self.view)
        self.view.setScene(self.scene)
        self.view.setRenderHint(QtGui.QPainter.Antialiasing)
        OnlinePlotter.tab_widget.addTab(self.view, tab_title)

        self._draw_topology()

    def _draw_topology(self):
        self.scene.clear()
        analysis = self.model._network_analysis
        if analysis is None:
            return

        components = sorted(self.model._components, key=lambda c: c.name)
        if not components:
            return

        positions = self._component_positions(components)
        edge_set = sorted({(src, dst) for src, dst, _, _ in analysis["edges"]}, key=lambda pair: (pair[0].name, pair[1].name))

        plan_colors = self._plan_color_map(analysis) if self.include_subnetworks else {}

        for src, dst in edge_set:
            self._draw_edge(positions[src], positions[dst])

        label_map = self._component_label_map(components)
        for comp in components:
            color = plan_colors.get(comp, QtGui.QColor(60, 120, 180))
            self._draw_node(positions[comp], label_map[comp], color)

        self.view.setSceneRect(self.scene.itemsBoundingRect().adjusted(-40, -40, 40, 40))
        self.view.fitInView(self.view.sceneRect(), QtCore.Qt.KeepAspectRatio)

    def _component_label_map(self, components):
        """Return stable, non-empty labels for all components.

        Components often have an empty `.name`; fall back to class name and add
        indices when there are duplicates.
        """
        base_names = []
        for comp in components:
            base = str(getattr(comp, "name", "") or "").strip()
            if not base:
                base = type(comp).__name__
            base_names.append(base)

        counts = {}
        for base in base_names:
            counts[base] = counts.get(base, 0) + 1

        seen = {}
        labels = {}
        for comp, base in zip(components, base_names):
            if counts[base] == 1:
                labels[comp] = base
                continue
            seen[base] = seen.get(base, 0) + 1
            labels[comp] = f"{base} {seen[base]}"
        return labels

    def _plan_color_map(self, analysis):
        color_map = {}
        palette = [
            QtGui.QColor(60, 120, 180),
            QtGui.QColor(230, 140, 30),
            QtGui.QColor(100, 170, 90),
            QtGui.QColor(170, 90, 170),
            QtGui.QColor(190, 80, 80),
            QtGui.QColor(90, 160, 160),
        ]

        for plan_index, plan in enumerate(analysis.get("plans", [])):
            base = palette[plan_index % len(palette)]
            if plan.get("mode") == "coupled":
                base = base.darker(105)
            for comp in plan.get("components", []):
                color_map[comp] = base
        return color_map

    def _component_positions(self, components):
        n = len(components)
        radius = 220 + 10 * max(n - 8, 0)
        center_x = 0
        center_y = 0
        positions = {}

        for i, comp in enumerate(components):
            angle = (2.0 * math.pi * i) / max(n, 1)
            x = center_x + radius * math.cos(angle)
            y = center_y + radius * math.sin(angle)
            positions[comp] = QtCore.QPointF(x, y)
        return positions

    def _draw_node(self, center, label, color):
        width = 150
        height = 56
        rect = QtCore.QRectF(center.x() - width / 2, center.y() - height / 2, width, height)
        pen = QtGui.QPen(QtGui.QColor(30, 30, 30))
        pen.setWidth(2)
        brush = QtGui.QBrush(color)
        path = QtGui.QPainterPath()
        path.addRoundedRect(rect, 10, 10)
        self.scene.addPath(path, pen, brush)

        text = self.scene.addText(label)
        text.setDefaultTextColor(QtGui.QColor(255, 255, 255))
        text_rect = text.boundingRect()
        text.setPos(center.x() - text_rect.width() / 2, center.y() - text_rect.height() / 2)

    def _draw_edge(self, src, dst):
        line_pen = QtGui.QPen(QtGui.QColor(70, 70, 70))
        line_pen.setWidth(2)

        line = QtCore.QLineF(src, dst)
        if line.length() < 1.0:
            return

        shrink = 82
        line.setLength(max(line.length() - shrink, 0.0))
        start = src
        end = line.p2()

        self.scene.addLine(QtCore.QLineF(start, end), line_pen)
        self._draw_arrowhead(start, end)

    def _draw_arrowhead(self, start, end):
        arrow_size = 10
        direction = QtCore.QLineF(start, end)
        angle = math.atan2(-direction.dy(), direction.dx())

        p1 = end + QtCore.QPointF(
            math.sin(angle - math.pi / 3) * arrow_size,
            math.cos(angle - math.pi / 3) * arrow_size,
        )
        p2 = end + QtCore.QPointF(
            math.sin(angle - math.pi + math.pi / 3) * arrow_size,
            math.cos(angle - math.pi + math.pi / 3) * arrow_size,
        )

        polygon = QtGui.QPolygonF([end, p1, p2])
        brush = QtGui.QBrush(QtGui.QColor(70, 70, 70))
        pen = QtGui.QPen(QtGui.QColor(70, 70, 70))
        self.scene.addPolygon(polygon, pen, brush)

    def export_png(self, path):
        image = QtGui.QImage(self.view.viewport().size(), QtGui.QImage.Format_ARGB32)
        image.fill(QtGui.QColor("white"))
        painter = QtGui.QPainter(image)
        self.view.render(painter)
        painter.end()
        image.save(path)

    def export_svg(self, path):
        try:
            from PyQt5.QtSvg import QSvgGenerator
        except Exception as exc:
            raise RuntimeError("SVG export requires PyQt5 QtSvg support.") from exc

        generator = QSvgGenerator()
        generator.setFileName(path)
        rect = self.view.sceneRect().toRect()
        generator.setSize(rect.size())
        generator.setViewBox(rect)

        painter = QtGui.QPainter(generator)
        self.scene.render(painter)
        painter.end()
