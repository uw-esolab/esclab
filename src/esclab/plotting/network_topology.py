import math

from pyqtgraph.Qt import QtCore, QtGui, QtWidgets

from esclab.plotting.online_plotter import OnlinePlotter


class _EdgeLabelItem(QtWidgets.QGraphicsTextItem):
    """Edge label text item with built-in hover tooltip."""

    def __init__(self, display_text, tooltip_text):
        super().__init__(display_text)
        self.setAcceptHoverEvents(True)
        self.setToolTip(tooltip_text)
        self.setFlag(QtWidgets.QGraphicsItem.ItemIgnoresTransformations, True)


class _TopologyGraphicsView(QtWidgets.QGraphicsView):
    """Interactive graphics view for topology: wheel zoom + drag pan."""

    ZOOM_IN_FACTOR = 1.15
    ZOOM_OUT_FACTOR = 1.0 / ZOOM_IN_FACTOR
    MIN_SCALE = 0.02
    MAX_SCALE = 50.0

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.transform_changed_callback = None
        self.setTransformationAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
        self.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
        self.setViewportUpdateMode(QtWidgets.QGraphicsView.SmartViewportUpdate)

    def wheelEvent(self, event):
        delta = event.angleDelta().y()
        if delta == 0:
            event.ignore()
            return

        factor = self.ZOOM_IN_FACTOR if delta > 0 else self.ZOOM_OUT_FACTOR
        current_scale = self.transform().m11()
        if current_scale <= 0.0:
            current_scale = 1.0
        target_scale = max(self.MIN_SCALE, min(self.MAX_SCALE, current_scale * factor))
        scale_factor = target_scale / current_scale

        self.scale(scale_factor, scale_factor)
        if self.transform_changed_callback is not None:
            self.transform_changed_callback()
        event.accept()


class NetworkTopologyView:
    NODE_WIDTH = 150
    NODE_HEIGHT = 56
    CONNECTION_LABEL_SEPARATOR = "→"
    # Alternative if preferred: "▶"
    LEADER_MID_BAND_MIN = 0.40
    LEADER_MID_BAND_MAX = 0.60
    LEADER_BIAS_T_LOW = 0.35
    LEADER_BIAS_T_HIGH = 0.65

    def __init__(self, model, tab_title="Network", include_subnetworks=True, show_connection_labels=True):
        self.model = model
        self.include_subnetworks = include_subnetworks
        self.show_connection_labels = show_connection_labels
        OnlinePlotter.ensure_window()

        self.view = _TopologyGraphicsView()
        self.scene = QtWidgets.QGraphicsScene(self.view)
        self._label_items = []
        self._node_rects = {}
        self._label_occupied_rects = []
        self._edge_label_links = []
        self.view.setScene(self.scene)
        self.view.setRenderHint(QtGui.QPainter.Antialiasing)
        OnlinePlotter.tab_widget.addTab(self.view, tab_title)
        OnlinePlotter.register_font_target(self)
        self.view.transform_changed_callback = self._update_edge_leaders

        self._draw_topology()

    def _draw_topology(self):
        self.scene.clear()
        self._label_items = []
        self._node_rects = {}
        self._label_occupied_rects = []
        self._edge_label_links = []
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
        self._refresh_node_label_fonts()
        self._update_edge_leaders()

    def apply_font_size(self):
        # Re-layout edge labels so collision avoidance stays valid after font changes.
        self._draw_topology()

    def _refresh_node_label_fonts(self):
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
            label = (
                f"{self._short_port_name(source_output.name)} "
                f"{self.CONNECTION_LABEL_SEPARATOR} "
                f"{self._short_port_name(destination_input.name)}"
            )
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
        text.setFlag(QtWidgets.QGraphicsItem.ItemIgnoresTransformations, True)
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

        display_text, tooltip_text = self._build_smart_label_text(edge_labels)
        text_item = _EdgeLabelItem(display_text, tooltip_text)
        self.scene.addItem(text_item)
        text_item.setDefaultTextColor(QtGui.QColor(20, 20, 20))

        font = QtGui.QFont()
        font.setPointSize(OnlinePlotter.font_size_pt)
        text_item.setFont(font)

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
        placed_rect = self._place_label_rect(center, text_rect, nx, ny)
        text_item.setPos(placed_rect.left(), placed_rect.top())
        self._register_edge_leader(text_item, start, end)
        self._label_occupied_rects.append(self._expand_rect(placed_rect, 4.0))

    def _register_edge_leader(self, text_item, line_start, line_end):
        pen = QtGui.QPen(QtGui.QColor(150, 150, 150, 210))
        pen.setWidth(1)
        leader_item = self.scene.addLine(QtCore.QLineF(line_start, line_start), pen)
        self._edge_label_links.append(
            {
                "text_item": text_item,
                "leader_item": leader_item,
                "line_start": QtCore.QPointF(line_start.x(), line_start.y()),
                "line_end": QtCore.QPointF(line_end.x(), line_end.y()),
            }
        )

    def _update_edge_leaders(self):
        for link in self._edge_label_links:
            self._update_one_edge_leader(link)

    def _update_one_edge_leader(self, link):
        text_item = link["text_item"]
        line_start = link["line_start"]
        line_end = link["line_end"]

        label_rect = self._label_visual_rect_in_scene(text_item)
        label_center = label_rect.center()
        anchor = self._leader_anchor_point(label_center, line_start, line_end)
        leader_start = self._ray_exit_rect(label_rect, label_center, anchor)

        line_item = link["leader_item"]
        if QtCore.QLineF(leader_start, anchor).length() < 3.0:
            line_item.setVisible(False)
            return

        line_item.setVisible(True)
        line_item.setLine(QtCore.QLineF(leader_start, anchor))

    def _label_visual_rect_in_scene(self, text_item):
        """Return the label's rendered scene rect accounting for ignore-transform mode."""
        rect = text_item.boundingRect()
        sx = abs(self.view.transform().m11())
        sy = abs(self.view.transform().m22())
        if sx < 1e-9:
            sx = 1.0
        if sy < 1e-9:
            sy = 1.0

        pos = text_item.scenePos()
        return QtCore.QRectF(
            pos.x(),
            pos.y(),
            rect.width() / sx,
            rect.height() / sy,
        )

    @classmethod
    def _leader_anchor_point(cls, label_center, seg_start, seg_end):
        """Pick a connector anchor point that avoids ambiguous center hotspots.

        Start from the orthogonal projection on the segment. If that projection lands
        near the segment midpoint, bias the anchor toward one side using the label's
        side-of-line sign so opposing labels on crossing lines do not both point to
        the same intersection point.
        """
        vx = seg_end.x() - seg_start.x()
        vy = seg_end.y() - seg_start.y()
        denom = vx * vx + vy * vy
        if denom <= 1e-12:
            return QtCore.QPointF(seg_start.x(), seg_start.y())

        t = ((label_center.x() - seg_start.x()) * vx + (label_center.y() - seg_start.y()) * vy) / denom
        t = max(0.0, min(1.0, t))

        if cls.LEADER_MID_BAND_MIN <= t <= cls.LEADER_MID_BAND_MAX:
            mid_x = (seg_start.x() + seg_end.x()) * 0.5
            mid_y = (seg_start.y() + seg_end.y()) * 0.5
            wx = label_center.x() - mid_x
            wy = label_center.y() - mid_y
            cross = vx * wy - vy * wx
            t = cls.LEADER_BIAS_T_LOW if cross >= 0.0 else cls.LEADER_BIAS_T_HIGH

        return QtCore.QPointF(seg_start.x() + t * vx, seg_start.y() + t * vy)

    @staticmethod
    def _build_smart_label_text(edge_labels, max_lines=2):
        full_text = "\n".join(edge_labels)
        if len(edge_labels) <= max_lines:
            return full_text, full_text

        kept = list(edge_labels[:max_lines])
        omitted = len(edge_labels) - max_lines
        if omitted >= 2:
            kept.append(f"... (+{omitted} more)")
        else:
            kept.append(edge_labels[max_lines])
        display_text = "\n".join(kept)
        return display_text, full_text

    def _place_label_rect(self, center, text_rect, nx, ny):
        base_offset = 12.0
        step = 12.0
        max_tries = 12

        best_rect = None
        best_penalty = float("inf")

        candidates = [0.0]
        for i in range(1, max_tries + 1):
            delta = i * step
            candidates.append(delta)
            candidates.append(-delta)

        for delta in candidates:
            offset = base_offset + delta
            cx = center.x() + nx * offset
            cy = center.y() + ny * offset
            rect = QtCore.QRectF(
                cx - text_rect.width() / 2,
                cy - text_rect.height() / 2,
                text_rect.width(),
                text_rect.height(),
            )
            penalty = self._placement_penalty(rect)
            if penalty < best_penalty:
                best_penalty = penalty
                best_rect = rect
                if penalty == 0.0:
                    break

        return best_rect

    def _placement_penalty(self, rect):
        padded_rect = self._expand_rect(rect, 2.0)
        penalty = 0.0

        for node_rect in self._node_rects.values():
            if padded_rect.intersects(node_rect):
                penalty += self._intersection_area(padded_rect, node_rect) + 1000.0

        for used_rect in self._label_occupied_rects:
            if padded_rect.intersects(used_rect):
                penalty += self._intersection_area(padded_rect, used_rect) + 500.0

        return penalty

    @staticmethod
    def _expand_rect(rect, margin):
        return QtCore.QRectF(
            rect.left() - margin,
            rect.top() - margin,
            rect.width() + 2 * margin,
            rect.height() + 2 * margin,
        )

    @staticmethod
    def _intersection_area(a, b):
        inter = a.intersected(b)
        if inter.isEmpty():
            return 0.0
        return inter.width() * inter.height()

    @staticmethod
    def _project_point_to_segment(point, seg_start, seg_end):
        """Project a point onto a line segment and clamp to segment bounds."""
        vx = seg_end.x() - seg_start.x()
        vy = seg_end.y() - seg_start.y()
        denom = vx * vx + vy * vy
        if denom <= 1e-12:
            return QtCore.QPointF(seg_start.x(), seg_start.y())

        t = ((point.x() - seg_start.x()) * vx + (point.y() - seg_start.y()) * vy) / denom
        t = max(0.0, min(1.0, t))
        return QtCore.QPointF(seg_start.x() + t * vx, seg_start.y() + t * vy)

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
