import math

from pyqtgraph.Qt import QtCore, QtGui, QtWidgets

from esclab.plotting.online_plotter import OnlinePlotter


class NetworkTopologyView:
    NODE_WIDTH = 150
    NODE_HEIGHT = 56

    def __init__(self, model, tab_title="Network", include_subnetworks=True, show_connection_labels=True):
        self.model = model
        self.include_subnetworks = include_subnetworks
        self.show_connection_labels = show_connection_labels
        OnlinePlotter.ensure_window()

        self.view = QtWidgets.QGraphicsView()
        self.scene = QtWidgets.QGraphicsScene(self.view)
        self._label_items = []
        self._node_rects = {}
        self.view.setScene(self.scene)
        self.view.setRenderHint(QtGui.QPainter.Antialiasing)
        OnlinePlotter.tab_widget.addTab(self.view, tab_title)
        OnlinePlotter.register_font_target(self)

        self._draw_topology()

    def _draw_topology(self):
        self.scene.clear()
        self._label_items = []
        self._node_rects = {}
        analysis = self.model._network_analysis
        if analysis is None:
            return

        components = sorted(self.model._components, key=lambda c: c.name)
        if not components:
            return

        positions = self._component_positions(components)
        edge_labels_by_pair = self._edge_label_groups(analysis)
        edge_set = sorted(edge_labels_by_pair.keys(), key=lambda pair: (pair[0].name, pair[1].name))

        plan_colors = self._plan_color_map(analysis) if self.include_subnetworks else {}

        label_map = self._component_label_map(components)
        for comp in components:
            color = plan_colors.get(comp, QtGui.QColor(60, 120, 180))
            self._draw_node(comp, positions[comp], label_map[comp], color)

        for src, dst in edge_set:
            self._draw_edge(src, dst, positions[src], positions[dst], edge_labels_by_pair.get((src, dst), ()))

        self.view.setSceneRect(self.scene.itemsBoundingRect().adjusted(-40, -40, 40, 40))
        self.view.fitInView(self.view.sceneRect(), QtCore.Qt.KeepAspectRatio)
        self.apply_font_size()

    def apply_font_size(self):
        font = QtGui.QFont()
        font.setPointSize(OnlinePlotter.font_size_pt)
        for text_item in self._label_items:
            text_item.setFont(font)

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

    def _edge_label_groups(self, analysis):
        grouped = {}
        for src, dst, source_output, destination_input in analysis.get("edges", []):
            label = f"{self._short_port_name(source_output.name)} -> {self._short_port_name(destination_input.name)}"
            key = (src, dst)
            if key not in grouped:
                grouped[key] = []
            if label not in grouped[key]:
                grouped[key].append(label)
        return {key: tuple(labels) for key, labels in grouped.items()}

    @staticmethod
    def _short_port_name(full_name):
        name = str(full_name or "")
        if "." in name:
            return name.split(".", 1)[1]
        return name

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

    def _draw_node(self, component, center, label, color):
        width = self.NODE_WIDTH
        height = self.NODE_HEIGHT
        rect = QtCore.QRectF(center.x() - width / 2, center.y() - height / 2, width, height)
        pen = QtGui.QPen(QtGui.QColor(30, 30, 30))
        pen.setWidth(2)
        brush = QtGui.QBrush(color)
        path = QtGui.QPainterPath()
        path.addRoundedRect(rect, 10, 10)
        self.scene.addPath(path, pen, brush)
        self._node_rects[component] = rect

        text = self.scene.addText(label)
        text.setDefaultTextColor(QtGui.QColor(255, 255, 255))
        self._label_items.append(text)
        text_rect = text.boundingRect()
        text.setPos(center.x() - text_rect.width() / 2, center.y() - text_rect.height() / 2)

    def _draw_edge(self, src_component, dst_component, src, dst, edge_labels):
        line_pen = QtGui.QPen(QtGui.QColor(70, 70, 70))
        line_pen.setWidth(2)

        line = QtCore.QLineF(src, dst)
        if line.length() < 1.0:
            return

        src_rect = self._node_rects.get(src_component)
        dst_rect = self._node_rects.get(dst_component)
        if src_rect is None or dst_rect is None:
            return

        start = self._ray_exit_rect(src_rect, src, dst)
        end = self._ray_exit_rect(dst_rect, dst, src)
        if QtCore.QLineF(start, end).length() < 1.0:
            return

        self.scene.addLine(QtCore.QLineF(start, end), line_pen)
        self._draw_arrowhead(start, end)
        self._draw_edge_label(start, end, edge_labels)

    def _draw_edge_label(self, start, end, edge_labels):
        if not self.show_connection_labels or not edge_labels:
            return

        text = "\n".join(edge_labels)
        text_item = self.scene.addText(text)
        text_item.setDefaultTextColor(QtGui.QColor(20, 20, 20))
        self._label_items.append(text_item)

        line = QtCore.QLineF(start, end)
        length = line.length()
        if length < 1e-9:
            return

        mid = QtCore.QPointF((start.x() + end.x()) * 0.5, (start.y() + end.y()) * 0.5)
        nx = -line.dy() / length
        ny = line.dx() / length
        offset = 12.0
        center = QtCore.QPointF(mid.x() + nx * offset, mid.y() + ny * offset)

        text_rect = text_item.boundingRect()
        text_item.setPos(center.x() - text_rect.width() / 2, center.y() - text_rect.height() / 2)

    @staticmethod
    def _ray_exit_rect(rect, center, toward):
        """Return the point where a ray from center toward toward exits rect.

        Assumes center is inside rect and rect is axis-aligned.
        """
        x0 = center.x()
        y0 = center.y()
        dx = toward.x() - x0
        dy = toward.y() - y0
        if abs(dx) < 1e-12 and abs(dy) < 1e-12:
            return QtCore.QPointF(x0, y0)

        candidates = []

        if dx > 1e-12:
            t = (rect.right() - x0) / dx
            candidates.append(t)
        elif dx < -1e-12:
            t = (rect.left() - x0) / dx
            candidates.append(t)

        if dy > 1e-12:
            t = (rect.bottom() - y0) / dy
            candidates.append(t)
        elif dy < -1e-12:
            t = (rect.top() - y0) / dy
            candidates.append(t)

        t_pos = [t for t in candidates if t > 0.0]
        if not t_pos:
            return QtCore.QPointF(x0, y0)
        t_exit = min(t_pos)
        return QtCore.QPointF(x0 + dx * t_exit, y0 + dy * t_exit)

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
