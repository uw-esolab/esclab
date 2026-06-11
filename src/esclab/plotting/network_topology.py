import json
import math
from html import escape

import numpy as np
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets

from esclab.plotting.online_plotter import OnlinePlotter


class _EdgeLabelItem(QtWidgets.QGraphicsTextItem):
    """Edge label text item with built-in hover tooltip."""

    def __init__(self, display_text, tooltip_text=None, tooltip_builder=None):
        super().__init__(display_text)
        self.setAcceptHoverEvents(True)
        self._tooltip_builder = tooltip_builder
        if tooltip_text is not None:
            self.setToolTip(tooltip_text)
        self.setFlag(QtWidgets.QGraphicsItem.ItemIgnoresTransformations, True)

    def hoverEnterEvent(self, event):
        if self._tooltip_builder is not None:
            self.setToolTip(self._tooltip_builder())
        super().hoverEnterEvent(event)


class _DraggableNodeItem(QtWidgets.QGraphicsPathItem):
    """A node shape that can be dragged to reposition it in the topology view."""

    def __init__(self, path, pen, brush, comp, topology_view, text_item):
        super().__init__(path)
        self.setPen(pen)
        self.setBrush(brush)
        self._comp = comp
        self._topology_view = topology_view
        self._text_item = text_item
        self.setFlags(
            QtWidgets.QGraphicsItem.ItemIsMovable
            | QtWidgets.QGraphicsItem.ItemSendsGeometryChanges
        )
        self.setCursor(QtCore.Qt.SizeAllCursor)

    def itemChange(self, change, value):
        if change == QtWidgets.QGraphicsItem.ItemPositionHasChanged:
            center = self.sceneBoundingRect().center()
            self._topology_view._center_text_item(self._text_item, center)
        return super().itemChange(change, value)

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)
        center = self.sceneBoundingRect().center()
        self._topology_view._on_node_drag_finished(self._comp, center)


class _DraggableEdgeSegmentItem(QtWidgets.QGraphicsLineItem):
    """Invisible draggable hit target for an edge segment."""

    def __init__(self, line, edge_key, segment_index, orientation, topology_view):
        super().__init__(line)
        self._edge_key = edge_key
        self._segment_index = segment_index
        self._orientation = orientation
        self._topology_view = topology_view

        pen = QtGui.QPen(QtGui.QColor(0, 0, 0, 0))
        pen.setWidth(12)
        pen.setCosmetic(True)
        self.setPen(pen)
        self.setFlags(
            QtWidgets.QGraphicsItem.ItemIsMovable
            | QtWidgets.QGraphicsItem.ItemSendsGeometryChanges
        )
        if self._orientation == "vertical":
            self.setCursor(QtCore.Qt.SizeHorCursor)
        else:
            self.setCursor(QtCore.Qt.SizeVerCursor)
        self.setAcceptedMouseButtons(QtCore.Qt.LeftButton)
        self.setZValue(1.5)

    def itemChange(self, change, value):
        if change == QtWidgets.QGraphicsItem.ItemPositionChange:
            pos = QtCore.QPointF(value)
            if self._orientation == "vertical":
                pos.setY(0.0)
            else:
                pos.setX(0.0)
            return pos
        return super().itemChange(change, value)

    def mouseReleaseEvent(self, event):
        super().mouseReleaseEvent(event)
        pos = self.pos()
        delta = pos.x() if self._orientation == "vertical" else pos.y()
        self._topology_view._on_edge_segment_drag_finished(
            self._edge_key,
            self._segment_index,
            delta,
        )


class _TopologyGraphicsView(QtWidgets.QGraphicsView):
    """Interactive graphics view for topology: wheel zoom + drag pan."""

    ZOOM_IN_FACTOR = 1.15
    ZOOM_OUT_FACTOR = 1.0 / ZOOM_IN_FACTOR
    MIN_SCALE = 0.02
    MAX_SCALE = 50.0

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.transform_changed_callback = None
        self._panning = False
        self._pan_start = QtCore.QPoint()
        self._pan_hbar_start = 0
        self._pan_vbar_start = 0
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
        self.setResizeAnchor(QtWidgets.QGraphicsView.AnchorUnderMouse)
        self.setDragMode(QtWidgets.QGraphicsView.NoDrag)
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

    def mousePressEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton:
            scene_pos = self.mapToScene(event.pos())
            for item in self.scene().items(scene_pos):
                if isinstance(item, (_DraggableNodeItem, _DraggableEdgeSegmentItem)):
                    super().mousePressEvent(event)
                    return
            self._panning = True
            self._pan_start = event.pos()
            self._pan_hbar_start = self.horizontalScrollBar().value()
            self._pan_vbar_start = self.verticalScrollBar().value()
            self.setCursor(QtCore.Qt.ClosedHandCursor)
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if self._panning:
            delta = event.pos() - self._pan_start
            self.horizontalScrollBar().setValue(self._pan_hbar_start - delta.x())
            self.verticalScrollBar().setValue(self._pan_vbar_start - delta.y())
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if event.button() == QtCore.Qt.LeftButton and self._panning:
            self._panning = False
            self.setCursor(QtCore.Qt.ArrowCursor)
            event.accept()
            return
        super().mouseReleaseEvent(event)


class NetworkTopologyView:
    NODE_WIDTH = 140
    NODE_HEIGHT = 35
    SCENE_RECT_PADDING = 40.0
    RECT_COLOR = (60, 180, 200)
    INITIAL_VIEW_FIT_MARGIN = 18.0
    INITIAL_VIEW_SCALE_BOOST = 1.0
    LAYOUT_ASPECT_STRETCH_TRIGGER = 1.00
    MAX_VERTICAL_LAYOUT_STRETCH = 1.50
    CONNECTION_LABEL_SEPARATOR = "→"
    # Alternative if preferred: "▶"
    LEADER_MID_BAND_MIN = 0.40
    LEADER_MID_BAND_MAX = 0.60
    LEADER_BIAS_T_LOW = 0.35
    LEADER_BIAS_T_HIGH = 1-LEADER_BIAS_T_LOW
    LABEL_LINE_OVERLAP_PENALTY = 900.0
    PID_LAYER_X_SPACING = 250.0
    PID_LAYER_Y_SPACING = 120.0
    PID_BARYCENTRIC_SWEEPS = 4
    EDGE_ORTH_CLEARANCE = 18.0
    EDGE_ROUTE_LANE_SPACING = 12.0
    EDGE_ROUTE_LAYER_CHANNEL_SPACING = 10.0
    EDGE_NODE_AVOID_MARGIN = 6.0
    EDGE_NODE_INTERSECTION_PENALTY = 10000.0
    SNAP_GRID_SIZE_PX = 20.0

    def __init__(
        self,
        model,
        tab_title="Network",
        include_subnetworks=True,
        show_connection_labels=True,
        layout_file=None,
    ):
        self.model = model
        self.include_subnetworks = include_subnetworks
        self.show_connection_labels = show_connection_labels
        OnlinePlotter.ensure_window()

        self.view = _TopologyGraphicsView()
        self.scene = QtWidgets.QGraphicsScene(self.view)
        self._label_items = []
        self._node_label_items = []
        self._node_rects = {}
        self._node_items = {}
        self._component_layout_keys = {}
        self._override_positions = {}
        self._edge_path_overrides = {}
        self._rendered_edge_paths = {}
        self._label_occupied_rects = []
        self._edge_label_links = []
        self._connector_segments = []
        self._layout_layers = {}
        self._edge_route_offsets = {}
        self._feedback_edge_pairs = set()
        self._layout_loopback_y = 0.0
        self._snap_to_grid_enabled = False
        self._snap_grid_size_px = float(self.SNAP_GRID_SIZE_PX)
        self._snap_grid_spinbox = None
        self.view.setScene(self.scene)
        self.view.setRenderHint(QtGui.QPainter.Antialiasing)

        # Build a container widget with a slim toolbar above the graphics view.
        container = QtWidgets.QWidget()
        vlayout = QtWidgets.QVBoxLayout(container)
        vlayout.setContentsMargins(0, 2, 0, 0)
        vlayout.setSpacing(2)
        toolbar = QtWidgets.QWidget()
        hlayout = QtWidgets.QHBoxLayout(toolbar)
        hlayout.setContentsMargins(4, 2, 4, 2)
        hlayout.setSpacing(6)
        btn_save = QtWidgets.QPushButton("Save Layout")
        btn_load = QtWidgets.QPushButton("Load Layout")
        btn_reset = QtWidgets.QPushButton("Reset Layout")
        btn_snap = QtWidgets.QPushButton("Snap to Grid")
        lbl_grid = QtWidgets.QLabel("Grid")
        spin_grid = QtWidgets.QSpinBox()
        btn_save.setFixedHeight(22)
        btn_load.setFixedHeight(22)
        btn_reset.setFixedHeight(22)
        btn_snap.setFixedHeight(22)
        lbl_grid.setFixedHeight(22)
        spin_grid.setFixedHeight(22)
        spin_grid.setRange(1, 1000)
        spin_grid.setSingleStep(5)
        spin_grid.setSuffix(" px")
        spin_grid.setValue(int(round(self._snap_grid_size_px)))
        spin_grid.setEnabled(self._snap_to_grid_enabled)
        spin_grid.valueChanged.connect(self._set_snap_grid_size)
        btn_snap.setCheckable(True)
        btn_snap.setChecked(self._snap_to_grid_enabled)
        btn_snap.toggled.connect(self._toggle_snap_to_grid)
        self._snap_grid_spinbox = spin_grid
        btn_save.clicked.connect(self._toolbar_save)
        btn_load.clicked.connect(self._toolbar_load)
        btn_reset.clicked.connect(self._toolbar_reset)
        hlayout.addWidget(btn_save)
        hlayout.addWidget(btn_load)
        hlayout.addWidget(btn_reset)
        hlayout.addWidget(btn_snap)
        hlayout.addWidget(lbl_grid)
        hlayout.addWidget(spin_grid)
        hlayout.addStretch()
        toolbar.setLayout(hlayout)
        vlayout.addWidget(toolbar)
        vlayout.addWidget(self.view)
        container.setLayout(vlayout)
        OnlinePlotter.tab_widget.addTab(container, tab_title)
        OnlinePlotter.register_font_target(self)
        self.view.transform_changed_callback = self._on_view_transform_changed

        if layout_file is not None:
            self._load_layout_file(layout_file)
        self._draw_topology()

    def _draw_topology(self):
        self.scene.clear()
        self._label_items = []
        self._node_label_items = []
        self._node_rects = {}
        self._node_items = {}
        self._component_layout_keys = {}
        self._rendered_edge_paths = {}
        self._label_occupied_rects = []
        self._edge_label_links = []
        self._connector_segments = []
        self._layout_layers = {}
        self._edge_route_offsets = {}
        self._feedback_edge_pairs = set()
        self._layout_loopback_y = 0.0
        analysis = self.model._network_analysis
        if analysis is None:
            return

        components = sorted(self.model._components, key=lambda c: c.name)
        if not components:
            return

        self._component_layout_keys = self._component_layout_key_map(components)

        positions = self._component_positions_pid_flow(components, analysis)
        positions = self._stretch_layout_for_viewport(positions)
        for comp in components:
            key = self._component_layout_keys.get(comp)
            if key in self._override_positions:
                x, y = self._override_positions[key]
                positions[comp] = QtCore.QPointF(x, y)
        if positions:
            y_vals = [p.y() for p in positions.values()]
            self._layout_loopback_y = max(y_vals) + self.PID_LAYER_Y_SPACING * 1.5
        edge_connections_by_pair = self._edge_connection_groups(analysis)
        edge_set = sorted(edge_connections_by_pair.keys(), key=lambda pair: (pair[0].name, pair[1].name))
        self._edge_route_offsets = self._compute_edge_route_offsets(edge_set)

        plan_colors = self._plan_color_map(analysis) if self.include_subnetworks else {}

        label_map = self._component_label_map(components)
        for comp in components:
            color = plan_colors.get(comp, QtGui.QColor(*self.RECT_COLOR))
            self._draw_node(comp, positions[comp], label_map[comp], color)

        edge_draw_data = []
        for src, dst in edge_set:
            is_feedback = (src, dst) in self._feedback_edge_pairs
            auto_path_points = self._compute_edge_path(
                src,
                dst,
                positions[src],
                positions[dst],
                route_offset=self._edge_route_offsets.get((src, dst), 0.0),
                is_feedback=is_feedback,
            )
            if auto_path_points is None or len(auto_path_points) < 2:
                continue
            edge_key = self._edge_layout_key(src, dst)
            path_points = self._apply_edge_path_override(edge_key, auto_path_points)
            if path_points is None or len(path_points) < 2:
                continue
            self._rendered_edge_paths[edge_key] = self._copy_path_points(path_points)
            label_segment = self._label_segment_for_path(path_points)
            segments = self._path_segments(path_points)
            edge_draw_data.append(
                (
                    edge_key,
                    src,
                    dst,
                    path_points,
                    label_segment,
                    segments,
                    edge_connections_by_pair.get((src, dst), ()),
                    label_map[src],
                    label_map[dst],
                    is_feedback,
                )
            )
            # Feedback loopback arcs run outside the node field — don't include
            # them in the segment registry used for forward-edge label avoidance.
            if not is_feedback:
                self._connector_segments.extend(segments)

        for edge_key, src, dst, path_points, label_segment, segments, edge_connections, src_label, dst_label, is_feedback in edge_draw_data:
            self._draw_edge(
                edge_key,
                src,
                dst,
                path_points,
                label_segment,
                segments,
                edge_connections,
                src_label,
                dst_label,
                is_feedback=is_feedback,
            )

        # Initial framing: fit to core graph geometry (nodes + routed edges),
        # not to label/leader extents, so the process network fills the view.
        graph_fit_rect = self._graph_fit_rect(edge_draw_data)
        scene_padding = self.SCENE_RECT_PADDING
        self.view.setSceneRect(self.scene.itemsBoundingRect().adjusted(-scene_padding, -scene_padding, scene_padding, scene_padding))
        if self.INITIAL_VIEW_SCALE_BOOST > 1.0:
            boost = self.INITIAL_VIEW_SCALE_BOOST
            current_scale = self.view.transform().m11()
            max_scale = self.view.MAX_SCALE if self.view.MAX_SCALE > 0.0 else current_scale * boost
            if current_scale > 0.0:
                capped_boost = min(boost, max_scale / current_scale)
                if capped_boost > 1.0:
                    self.view.scale(capped_boost, capped_boost)
        self._refresh_node_label_fonts()

        # Place edge labels after fitInView so label scene extents are computed
        # with the final transform, keeping label offsets stable across zoom.
        for edge_key, src, dst, path_points, label_segment, segments, edge_connections, src_label, dst_label, is_feedback in edge_draw_data:
            if label_segment is None:
                continue
            self._draw_edge_label(
                label_segment[0],
                label_segment[1],
                edge_connections,
                src_label,
                dst_label,
                own_segments=segments,
            )

        self._update_edge_leaders()

        # Expand scene bounds after labels/leaders are in place so panning can
        # still reach all annotations even though initial fit ignores them.
        self.view.setSceneRect(self.scene.itemsBoundingRect().adjusted(-scene_padding, -scene_padding, scene_padding, scene_padding))

    def _graph_fit_rect(self, edge_draw_data):
        """Build a fit rect from nodes and edge paths only.

        Labels and leader lines are intentionally excluded so they do not force
        an overly zoomed-out initial framing. Feedback loopback arcs are also
        excluded because they are decorative return paths that can be very long
        and would otherwise dominate the initial zoom.
        """
        fit_rect = QtCore.QRectF()
        has_geometry = False

        for node_rect in self._node_rects.values():
            fit_rect = node_rect if not has_geometry else fit_rect.united(node_rect)
            has_geometry = True

        for _, _, _, path_points, _, _, _, _, _, is_feedback in edge_draw_data:
            if is_feedback:
                continue
            for point in path_points:
                point_rect = QtCore.QRectF(point.x(), point.y(), 0.0, 0.0)
                fit_rect = point_rect if not has_geometry else fit_rect.united(point_rect)
                has_geometry = True

        if not has_geometry:
            fallback = self.scene.itemsBoundingRect()
            if fallback.isNull():
                return QtCore.QRectF(-1.0, -1.0, 2.0, 2.0)
            fit_rect = fallback

        margin = self.INITIAL_VIEW_FIT_MARGIN
        return fit_rect.adjusted(-margin, -margin, margin, margin)

    def _stretch_layout_for_viewport(self, positions):
        """Vertically stretch wide layouts so they better fill the viewport.

        The layered process layout can be much wider than tall, which leaves
        large unused vertical space after fit-in-view. Stretching y around the
        layout centroid improves initial screen utilization while preserving left
        to right ordering and relative vertical ranks.
        """
        if not positions or len(positions) < 2:
            return positions

        xs = [point.x() for point in positions.values()]
        ys = [point.y() for point in positions.values()]
        width = (max(xs) - min(xs)) + self.NODE_WIDTH
        height = (max(ys) - min(ys)) + self.NODE_HEIGHT
        if width <= 1e-9 or height <= 1e-9:
            return positions

        graph_aspect = width / height
        viewport = self.view.viewport().size()
        viewport_width = max(1, viewport.width())
        viewport_height = max(1, viewport.height())
        target_aspect = viewport_width / viewport_height
        if target_aspect <= 1e-9:
            target_aspect = 16.0 / 9.0

        if graph_aspect <= target_aspect * self.LAYOUT_ASPECT_STRETCH_TRIGGER:
            return positions

        stretch = min(self.MAX_VERTICAL_LAYOUT_STRETCH, graph_aspect / target_aspect)
        if stretch <= 1.0:
            return positions

        y_center = 0.5 * (max(ys) + min(ys))
        return {
            comp: QtCore.QPointF(point.x(), y_center + (point.y() - y_center) * stretch)
            for comp, point in positions.items()
        }

    def apply_font_size(self):
        # Re-layout edge labels so collision avoidance stays valid after font changes.
        self._draw_topology()

    def _on_view_transform_changed(self):
        self._recenter_node_labels()
        self._update_edge_leaders()

    def _refresh_node_label_fonts(self):
        font = QtGui.QFont()
        font.setPointSize(OnlinePlotter.font_size_pt)
        for text_item in self._label_items:
            text_item.setFont(font)
        self._recenter_node_labels()

    def _recenter_node_labels(self):
        for text_item, node_item in self._node_label_items:
            self._center_text_item(text_item, node_item.sceneBoundingRect().center())

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

    def _component_layout_key_map(self, components):
        """Return stable, unique keys for JSON layout persistence.

        Use component names when available; otherwise fall back to class names,
        and append indices when duplicates exist.
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
        keys = {}
        for comp, base in zip(components, base_names):
            if counts[base] == 1:
                keys[comp] = base
                continue
            seen[base] = seen.get(base, 0) + 1
            keys[comp] = f"{base}__{seen[base]}"
        return keys

    def _edge_layout_key(self, src, dst):
        src_key = self._component_layout_keys.get(src)
        dst_key = self._component_layout_keys.get(dst)
        if src_key is None:
            src_key = str(getattr(src, "name", "") or "").strip() or type(src).__name__
        if dst_key is None:
            dst_key = str(getattr(dst, "name", "") or "").strip() or type(dst).__name__
        return f"{src_key}__TO__{dst_key}"

    @staticmethod
    def _copy_path_points(path_points):
        return [QtCore.QPointF(point.x(), point.y()) for point in path_points]

    def _apply_edge_path_override(self, edge_key, auto_path_points):
        override_points = self._edge_path_overrides.get(edge_key)
        if not override_points:
            return auto_path_points

        path_points = self._copy_path_points(override_points)
        if len(path_points) < 2:
            return auto_path_points

        path_points[0] = QtCore.QPointF(auto_path_points[0].x(), auto_path_points[0].y())
        path_points[-1] = QtCore.QPointF(auto_path_points[-1].x(), auto_path_points[-1].y())
        path_points = self._simplify_path(path_points)
        if len(path_points) < 2:
            return auto_path_points
        return path_points

    @staticmethod
    def _path_points_to_json(path_points):
        return [{"x": point.x(), "y": point.y()} for point in path_points]

    @staticmethod
    def _path_points_from_json(points_json):
        points = []
        for item in points_json:
            if not isinstance(item, dict):
                continue
            if "x" not in item or "y" not in item:
                continue
            points.append(QtCore.QPointF(float(item["x"]), float(item["y"])))
        return points

    def _plan_color_map(self, analysis):
        color_map = {}
        palette = [
            QtGui.QColor(*self.RECT_COLOR),
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

    def _edge_connection_groups(self, analysis):
        grouped = {}
        grouped_seen = {}
        for src, dst, source_output, destination_input in analysis.get("edges", []):
            key = (src, dst)
            if key not in grouped:
                grouped[key] = []
                grouped_seen[key] = set()

            full_pair = (str(source_output.name), str(destination_input.name))
            if full_pair in grouped_seen[key]:
                continue
            grouped_seen[key].add(full_pair)

            connection = destination_input.connection
            grouped[key].append(
                {
                    "source_full": str(source_output.name),
                    "dest_full": str(destination_input.name),
                    "source_short": self._short_port_name(source_output.name),
                    "dest_short": self._short_port_name(destination_input.name),
                    "connection": connection,
                }
            )

        return {key: tuple(connections) for key, connections in grouped.items()}

    @staticmethod
    def _short_port_name(full_name):
        name = str(full_name or "")
        if "." in name:
            return name.split(".", 1)[1]
        return name

    def _component_positions_pid_flow(self, components, analysis):
        if not components:
            return {}

        order_index = {component: idx for idx, component in enumerate(components)}
        nodes = list(components)
        adjacency, reverse_adjacency = self._build_unique_adjacency(nodes, analysis, order_index)
        if all(len(neighbors) == 0 for neighbors in adjacency.values()):
            # No edges: arrange in a single row.
            return {comp: QtCore.QPointF(i * self.PID_LAYER_X_SPACING, 0.0) for i, comp in enumerate(components)}

        feedback_edges = self._find_feedback_edges(nodes, adjacency, order_index)
        self._feedback_edge_pairs = set(feedback_edges)
        reduced_adjacency, reduced_reverse = self._reduce_adjacency(adjacency, reverse_adjacency, feedback_edges)
        layers = self._assign_layers(nodes, reduced_reverse, reduced_adjacency, order_index)
        self._layout_layers = dict(layers)
        layer_buckets = self._initial_layer_buckets(nodes, layers, order_index)
        self._barycentric_reorder(layer_buckets, layers, reduced_adjacency, reduced_reverse, order_index)

        max_height = max(len(layer_nodes) for layer_nodes in layer_buckets.values())
        y_spacing = self.PID_LAYER_Y_SPACING
        if max_height >= 8:
            y_spacing *= 0.9
        if max_height >= 12:
            y_spacing *= 0.85

        layer_y_positions = self._compute_layered_y_positions(
            layer_buckets,
            reduced_reverse,
            reduced_adjacency,
            layers,
            y_spacing,
            order_index,
        )

        positions = {}
        for layer_id in sorted(layer_buckets):
            layer_nodes = layer_buckets[layer_id]
            for rank, component in enumerate(layer_nodes):
                x = layer_id * self.PID_LAYER_X_SPACING
                y = layer_y_positions.get(component)
                if y is None:
                    n_layer = len(layer_nodes)
                    y = (rank - (n_layer - 1) * 0.5) * y_spacing
                positions[component] = QtCore.QPointF(x, y)

        return positions

    @staticmethod
    def _compute_layered_y_positions(layer_buckets, reverse_adjacency, adjacency, layers, y_spacing, order_index):
        y_positions = {}

        for layer_id in sorted(layer_buckets):
            layer_nodes = list(layer_buckets[layer_id])
            if not layer_nodes:
                continue

            n_layer = len(layer_nodes)
            default_targets = {
                node: (idx - (n_layer - 1) * 0.5) * y_spacing
                for idx, node in enumerate(layer_nodes)
            }

            targets = []
            for rank_index, node in enumerate(layer_nodes):
                preds = [pred for pred in reverse_adjacency[node] if pred in y_positions]
                if preds:
                    pred_mean = sum(y_positions[pred] for pred in preds) / float(len(preds))

                    branch_targets = []
                    for pred in sorted(preds, key=lambda item: order_index[item]):
                        pred_y = y_positions[pred]
                        siblings = [
                            child
                            for child in adjacency[pred]
                            if layers.get(child) == layer_id
                        ]
                        siblings.sort(key=lambda item: order_index[item])
                        if len(siblings) > 1 and node in siblings:
                            sibling_rank = siblings.index(node) - (len(siblings) - 1) * 0.5
                            branch_targets.append(pred_y + sibling_rank * y_spacing)
                        else:
                            branch_targets.append(pred_y)

                    branch_mean = sum(branch_targets) / float(len(branch_targets))
                    target = 0.65 * branch_mean + 0.20 * pred_mean + 0.15 * default_targets[node]
                else:
                    target = default_targets[node]
                targets.append((target, rank_index, node))

            targets.sort(key=lambda item: (item[0], item[1]))

            assigned = []
            assigned_targets = []
            for target, _rank_index, node in targets:
                y_val = target
                if assigned and y_val - assigned[-1][0] < y_spacing:
                    y_val = assigned[-1][0] + y_spacing
                assigned.append((y_val, node))
                assigned_targets.append(target)

            center_target = sum(assigned_targets) / float(len(assigned_targets))
            center_assigned = sum(y for y, _node in assigned) / float(len(assigned))
            shift = center_target - center_assigned

            for y_val, node in assigned:
                y_positions[node] = y_val + shift

        if len(y_positions) > 1:
            y_vals = list(y_positions.values())
            y_min = min(y_vals)
            y_max = max(y_vals)
            current_span = y_max - y_min
            target_span = y_spacing * max(2.0, min(12.0, 1.5 + 0.75 * math.sqrt(len(y_positions))))

            if current_span > 1e-9 and current_span < target_span:
                center = 0.5 * (y_min + y_max)
                scale = target_span / current_span
                for node in list(y_positions.keys()):
                    y_positions[node] = center + (y_positions[node] - center) * scale
            elif current_span <= 1e-9:
                # Degenerate case: all nodes collapsed into one y-level.
                for node in y_positions:
                    layer_id = layers.get(node, 0)
                    y_positions[node] = ((layer_id % 5) - 2.0) * (0.5 * y_spacing)

        return y_positions

    def _compute_edge_route_offsets(self, edge_set):
        outgoing = {}
        incoming = {}
        for src, dst in edge_set:
            outgoing.setdefault(src, []).append(dst)
            incoming.setdefault(dst, []).append(src)

        out_rank = {
            src: self._centered_rank_map(sorted(dsts, key=lambda comp: comp.name))
            for src, dsts in outgoing.items()
        }
        in_rank = {
            dst: self._centered_rank_map(sorted(srcs, key=lambda comp: comp.name))
            for dst, srcs in incoming.items()
        }

        layer_pair_edges = {}
        for src, dst in edge_set:
            src_layer = self._layout_layers.get(src, 0)
            dst_layer = self._layout_layers.get(dst, 0)
            layer_pair_edges.setdefault((src_layer, dst_layer), []).append((src, dst))

        pair_rank = {}
        for pair_key, edges in layer_pair_edges.items():
            sorted_edges = sorted(edges, key=lambda edge: (edge[0].name, edge[1].name))
            n_edges = len(sorted_edges)
            for idx, edge in enumerate(sorted_edges):
                pair_rank[edge] = idx - (n_edges - 1) * 0.5

        offsets = {}
        for src, dst in edge_set:
            src_rank = out_rank.get(src, {}).get(dst, 0.0)
            dst_rank = in_rank.get(dst, {}).get(src, 0.0)
            base = (src_rank + dst_rank) * self.EDGE_ROUTE_LANE_SPACING
            base += pair_rank.get((src, dst), 0.0) * self.EDGE_ROUTE_LAYER_CHANNEL_SPACING

            src_layer = self._layout_layers.get(src, 0)
            dst_layer = self._layout_layers.get(dst, 0)
            if dst_layer <= src_layer:
                # Backward/cycle edges get a larger offset lane to avoid stacking.
                if abs(base) < 1e-9:
                    base = self.EDGE_ROUTE_LANE_SPACING
                base *= 1.4

            offsets[(src, dst)] = base

        return offsets

    @staticmethod
    def _centered_rank_map(items):
        n_items = len(items)
        return {item: idx - (n_items - 1) * 0.5 for idx, item in enumerate(items)}

    @staticmethod
    def _build_unique_adjacency(nodes, analysis, order_index):
        adjacency = {node: set() for node in nodes}
        reverse_adjacency = {node: set() for node in nodes}
        dedupe = set()

        for src, dst, _source_output, _destination_input in analysis.get("edges", []):
            if src not in adjacency or dst not in adjacency:
                continue
            key = (order_index[src], order_index[dst])
            if key in dedupe:
                continue
            dedupe.add(key)
            adjacency[src].add(dst)
            reverse_adjacency[dst].add(src)

        return adjacency, reverse_adjacency

    @staticmethod
    def _reduce_adjacency(adjacency, reverse_adjacency, feedback_edges):
        reduced_adjacency = {node: set() for node in adjacency}
        reduced_reverse = {node: set() for node in reverse_adjacency}

        for src, neighbors in adjacency.items():
            for dst in neighbors:
                if (src, dst) in feedback_edges:
                    continue
                reduced_adjacency[src].add(dst)
                reduced_reverse[dst].add(src)

        return reduced_adjacency, reduced_reverse

    @staticmethod
    def _find_feedback_edges(nodes, adjacency, order_index):
        # Score each edge: in_degree(src) - out_degree(dst).  Higher = more
        # likely a return/recycle path (merge node flowing back to source node).
        # Build the forward graph greedily in score-ascending order; any edge
        # that would close a cycle becomes a feedback edge instead.
        in_count = {node: 0 for node in nodes}
        out_count = {node: len(adjacency[node]) for node in nodes}
        for src in nodes:
            for dst in adjacency[src]:
                in_count[dst] += 1

        all_edges = sorted(
            ((src, dst) for src in nodes for dst in adjacency[src]),
            key=lambda e: (in_count[e[0]] - out_count[e[1]], order_index[e[0]], order_index[e[1]]),
        )

        def reachable(start, target, adj):
            if start is target:
                return True
            visited, stack = {start}, [start]
            while stack:
                node = stack.pop()
                for nb in adj[node]:
                    if nb is target:
                        return True
                    if nb not in visited:
                        visited.add(nb)
                        stack.append(nb)
            return False

        fwd_adjacency = {node: set() for node in nodes}
        feedback_edges = set()
        for src, dst in all_edges:
            if reachable(dst, src, fwd_adjacency):
                feedback_edges.add((src, dst))
            else:
                fwd_adjacency[src].add(dst)

        return feedback_edges

    @staticmethod
    def _assign_layers(nodes, reverse_adjacency, adjacency, order_index):
        in_degree = {node: len(reverse_adjacency[node]) for node in nodes}
        queue = [node for node in nodes if in_degree[node] == 0]
        queue.sort(key=lambda item: order_index[item])

        topo_order = []
        while queue:
            node = queue.pop(0)
            topo_order.append(node)
            for neighbor in sorted(adjacency[node], key=lambda item: order_index[item]):
                in_degree[neighbor] -= 1
                if in_degree[neighbor] == 0:
                    queue.append(neighbor)
            queue.sort(key=lambda item: order_index[item])

        if len(topo_order) < len(nodes):
            emitted = set(topo_order)
            topo_order.extend(node for node in nodes if node not in emitted)

        layers = {node: 0 for node in nodes}
        for node in topo_order:
            preds = reverse_adjacency[node]
            if preds:
                layers[node] = max(layers[pred] + 1 for pred in preds)
        return layers

    @staticmethod
    def _initial_layer_buckets(nodes, layers, order_index):
        buckets = {}
        for node in nodes:
            buckets.setdefault(layers[node], []).append(node)
        for layer_id in buckets:
            buckets[layer_id].sort(key=lambda item: order_index[item])
        return buckets

    def _barycentric_reorder(self, layer_buckets, layers, adjacency, reverse_adjacency, order_index):
        max_layer = max(layer_buckets) if layer_buckets else -1
        if max_layer <= 0:
            return

        for _ in range(self.PID_BARYCENTRIC_SWEEPS):
            for layer_id in range(1, max_layer + 1):
                previous_layer = layer_buckets.get(layer_id - 1, [])
                if not previous_layer:
                    continue
                previous_rank = {node: idx for idx, node in enumerate(previous_layer)}
                self._reorder_layer(
                    layer_buckets,
                    layer_id,
                    previous_rank,
                    reverse_adjacency,
                    layers,
                    layer_id - 1,
                    order_index,
                )

            for layer_id in range(max_layer - 1, -1, -1):
                next_layer = layer_buckets.get(layer_id + 1, [])
                if not next_layer:
                    continue
                next_rank = {node: idx for idx, node in enumerate(next_layer)}
                self._reorder_layer(
                    layer_buckets,
                    layer_id,
                    next_rank,
                    adjacency,
                    layers,
                    layer_id + 1,
                    order_index,
                )

            self._local_swap_improvement(layer_buckets, layers, adjacency, reverse_adjacency)

    @staticmethod
    def _reorder_layer(layer_buckets, layer_id, neighbor_rank, connectivity, layers, expected_neighbor_layer, order_index):
        layer_nodes = list(layer_buckets.get(layer_id, []))
        if len(layer_nodes) <= 1:
            return

        old_rank = {node: idx for idx, node in enumerate(layer_nodes)}

        def barycenter(node):
            neighbors = [
                neighbor
                for neighbor in connectivity[node]
                if layers.get(neighbor) == expected_neighbor_layer and neighbor in neighbor_rank
            ]
            if not neighbors:
                return float(old_rank[node])
            return float(sum(neighbor_rank[neighbor] for neighbor in neighbors)) / float(len(neighbors))

        layer_nodes.sort(key=lambda node: (barycenter(node), old_rank[node], order_index[node]))
        layer_buckets[layer_id] = layer_nodes

    def _local_swap_improvement(self, layer_buckets, layers, adjacency, reverse_adjacency):
        if not layer_buckets:
            return

        for layer_id in sorted(layer_buckets):
            layer_nodes = list(layer_buckets[layer_id])
            if len(layer_nodes) <= 2:
                continue

            improved = True
            while improved:
                improved = False
                idx = 0
                while idx < len(layer_nodes) - 1:
                    base_score = self._layer_crossing_score(
                        layer_id,
                        layer_nodes,
                        layer_buckets,
                        layers,
                        adjacency,
                        reverse_adjacency,
                    )

                    candidate_nodes = list(layer_nodes)
                    candidate_nodes[idx], candidate_nodes[idx + 1] = candidate_nodes[idx + 1], candidate_nodes[idx]
                    candidate_score = self._layer_crossing_score(
                        layer_id,
                        candidate_nodes,
                        layer_buckets,
                        layers,
                        adjacency,
                        reverse_adjacency,
                    )

                    if candidate_score < base_score:
                        layer_nodes = candidate_nodes
                        improved = True
                    idx += 1

            layer_buckets[layer_id] = layer_nodes

    @staticmethod
    def _layer_crossing_score(layer_id, layer_nodes, layer_buckets, layers, adjacency, reverse_adjacency):
        rank_current = {node: idx for idx, node in enumerate(layer_nodes)}
        score = 0

        prev_nodes = layer_buckets.get(layer_id - 1, [])
        if prev_nodes:
            rank_prev = {node: idx for idx, node in enumerate(prev_nodes)}
            edge_pairs = []
            for src in prev_nodes:
                for dst in adjacency[src]:
                    if layers.get(dst) == layer_id and dst in rank_current:
                        edge_pairs.append((rank_prev[src], rank_current[dst]))
            score += NetworkTopologyView._count_pair_crossings(edge_pairs)

        next_nodes = layer_buckets.get(layer_id + 1, [])
        if next_nodes:
            rank_next = {node: idx for idx, node in enumerate(next_nodes)}
            edge_pairs = []
            for src in layer_nodes:
                for dst in adjacency[src]:
                    if layers.get(dst) == layer_id + 1 and dst in rank_next:
                        edge_pairs.append((rank_current[src], rank_next[dst]))
            score += NetworkTopologyView._count_pair_crossings(edge_pairs)

            # Also consider reverse links from next layer to current for coupled graphs.
            edge_pairs = []
            for dst in layer_nodes:
                for src in reverse_adjacency[dst]:
                    if layers.get(src) == layer_id + 1 and src in rank_next:
                        edge_pairs.append((rank_current[dst], rank_next[src]))
            score += NetworkTopologyView._count_pair_crossings(edge_pairs)

        return score

    @staticmethod
    def _count_pair_crossings(edge_pairs):
        crossings = 0
        n = len(edge_pairs)
        for i in range(n):
            a1, b1 = edge_pairs[i]
            for j in range(i + 1, n):
                a2, b2 = edge_pairs[j]
                if (a1 - a2) * (b1 - b2) < 0:
                    crossings += 1
        return crossings

    def _draw_node(self, component, center, label, color):
        width = self.NODE_WIDTH
        height = self.NODE_HEIGHT
        rect = QtCore.QRectF(center.x() - width / 2, center.y() - height / 2, width, height)
        pen = QtGui.QPen(QtGui.QColor(30, 30, 30))
        pen.setWidth(2)
        brush = QtGui.QBrush(color)
        path = QtGui.QPainterPath()
        path.addRoundedRect(rect, 10, 10)
        self._node_rects[component] = rect

        text = QtWidgets.QGraphicsTextItem(label)
        text.setDefaultTextColor(QtGui.QColor(255, 255, 255))
        text.setFlag(QtWidgets.QGraphicsItem.ItemIgnoresTransformations, True)
        text.setZValue(2.0)
        font = QtGui.QFont()
        font.setPointSize(OnlinePlotter.font_size_pt)
        text.setFont(font)
        self.scene.addItem(text)

        node_item = _DraggableNodeItem(path, pen, brush, component, self, text)
        node_item.setZValue(1.0)
        self.scene.addItem(node_item)
        self._node_items[component] = node_item
        self._label_items.append(text)
        self._node_label_items.append((text, node_item))
        self._center_text_item(text, center)

    def _center_text_item(self, text_item, center):
        # ItemIgnoresTransformations keeps text in device-space size; convert the
        # local text extents to scene-space by dividing by current view scale.
        text_rect = text_item.boundingRect()
        sx = abs(self.view.transform().m11())
        sy = abs(self.view.transform().m22())
        if sx < 1e-9:
            sx = 1.0
        if sy < 1e-9:
            sy = 1.0

        dx = 0.5 * text_rect.width() / sx
        dy = 0.5 * text_rect.height() / sy
        text_item.setPos(center.x() - dx, center.y() - dy)

    def _compute_edge_segment(self, src_component, dst_component, src, dst):
        line = QtCore.QLineF(src, dst)
        if line.length() < 1.0:
            return None

        src_rect = self._node_rects.get(src_component)
        dst_rect = self._node_rects.get(dst_component)
        if src_rect is None or dst_rect is None:
            return None

        start = self._ray_exit_rect(src_rect, src, dst)
        end = self._ray_exit_rect(dst_rect, dst, src)
        if QtCore.QLineF(start, end).length() < 1.0:
            return None

        return start, end

    def _compute_edge_path(self, src_component, dst_component, src, dst, route_offset=0.0, is_feedback=False):
        segment = self._compute_edge_segment(src_component, dst_component, src, dst)
        if segment is None:
            return None

        start, end = segment
        src_rect = self._node_rects.get(src_component)
        dst_rect = self._node_rects.get(dst_component)
        if src_rect is None or dst_rect is None:
            return [start, end]

        if is_feedback:
            return self._loopback_path(src_rect, dst_rect, src, dst)

        start, end = self._preferred_horizontal_anchors(src_rect, dst_rect, src, dst)
        path_points = self._orthogonal_path(start, end, route_offset=route_offset)
        path_points = self._reroute_path_around_nodes(
            path_points,
            start,
            end,
            src_component,
            dst_component,
            route_offset=route_offset,
        )
        if len(path_points) < 2:
            return [start, end]
        return path_points

    def _reroute_path_around_nodes(self, path_points, start, end, src_component, dst_component, route_offset=0.0):
        """Reroute H-V-H paths when their jog cuts through neighboring nodes."""
        if len(path_points) < 2:
            return path_points
        if abs(start.y() - end.y()) < 1e-9:
            return path_points

        obstacle_rects = []
        margin = self.EDGE_NODE_AVOID_MARGIN
        for component, rect in self._node_rects.items():
            if component is src_component or component is dst_component:
                continue
            obstacle_rects.append(self._expand_rect(rect, margin))

        if not obstacle_rects:
            return path_points
        if not self._path_intersects_rects(path_points, obstacle_rects):
            return path_points

        x_base = (start.x() + end.x()) * 0.5 + route_offset
        if end.x() >= start.x():
            x_min, x_max = start.x(), end.x()
            x_base = max(x_min, min(x_base, x_max))
        else:
            x_min, x_max = end.x(), start.x()
            x_base = max(x_min, min(x_base, x_max))

        clearance = max(self.EDGE_ORTH_CLEARANCE, self.EDGE_ORTH_CLEARANCE + 0.25 * abs(route_offset))
        lane_step = max(8.0, self.EDGE_ROUTE_LANE_SPACING)

        candidate_xs = [x_base]
        vertical_probe_start = QtCore.QPointF(x_base, start.y())
        vertical_probe_end = QtCore.QPointF(x_base, end.y())
        for rect in obstacle_rects:
            if self._segment_intersects_rect(vertical_probe_start, vertical_probe_end, rect):
                candidate_xs.append(rect.left() - clearance)
                candidate_xs.append(rect.right() + clearance)

        for i in range(1, 8):
            delta = i * lane_step
            candidate_xs.append(x_base - delta)
            candidate_xs.append(x_base + delta)

        unique_candidate_xs = []
        for x_val in candidate_xs:
            if all(abs(x_val - existing) > 1e-9 for existing in unique_candidate_xs):
                unique_candidate_xs.append(x_val)

        best_path = path_points
        best_score = self._edge_path_obstacle_score(best_path, obstacle_rects, x_base, x_base, x_min, x_max)
        for x_mid in unique_candidate_xs:
            candidate_path = self._simplify_path(
                [
                    start,
                    QtCore.QPointF(x_mid, start.y()),
                    QtCore.QPointF(x_mid, end.y()),
                    end,
                ]
            )
            score = self._edge_path_obstacle_score(candidate_path, obstacle_rects, x_mid, x_base, x_min, x_max)
            if score < best_score:
                best_score = score
                best_path = candidate_path

        return best_path

    def _edge_path_obstacle_score(self, path_points, obstacle_rects, x_mid, x_base, x_min, x_max):
        score = 0.0
        total_length = 0.0
        for seg_start, seg_end in self._path_segments(path_points):
            total_length += QtCore.QLineF(seg_start, seg_end).length()
            for rect in obstacle_rects:
                if self._segment_intersects_rect(seg_start, seg_end, rect):
                    score += self.EDGE_NODE_INTERSECTION_PENALTY

        if x_mid < x_min:
            score += 2000.0 + 50.0 * (x_min - x_mid)
        elif x_mid > x_max:
            score += 2000.0 + 50.0 * (x_mid - x_max)

        score += 0.02 * total_length
        score += 0.25 * abs(x_mid - x_base)
        return score

    def _path_intersects_rects(self, path_points, rects):
        for seg_start, seg_end in self._path_segments(path_points):
            for rect in rects:
                if self._segment_intersects_rect(seg_start, seg_end, rect):
                    return True
        return False

    def _loopback_path(self, src_rect, dst_rect, src_center, dst_center):
        """Route a return edge below the process spine (P&ID loopback convention)."""
        clearance = self.EDGE_ORTH_CLEARANCE
        start = QtCore.QPointF(src_rect.right(), src_center.y())
        end = QtCore.QPointF(dst_rect.left(), dst_center.y())
        x_out = src_rect.right() + clearance
        x_in = dst_rect.left() - clearance
        loopback_y = self._layout_loopback_y
        path = [
            start,
            QtCore.QPointF(x_out, src_center.y()),
            QtCore.QPointF(x_out, loopback_y),
            QtCore.QPointF(x_in, loopback_y),
            QtCore.QPointF(x_in, dst_center.y()),
            end,
        ]
        return self._simplify_path(path)

    @classmethod
    def _preferred_horizontal_anchors(cls, src_rect, dst_rect, src_center, dst_center):
        dx = dst_center.x() - src_center.x()
        dy = dst_center.y() - src_center.y()

        # For process-flow layout, prefer left/right entry/exit for readability.
        if abs(dx) >= 1e-9:
            if dx >= 0.0:
                start = QtCore.QPointF(src_rect.right(), src_center.y())
                end = QtCore.QPointF(dst_rect.left(), dst_center.y())
            else:
                start = QtCore.QPointF(src_rect.left(), src_center.y())
                end = QtCore.QPointF(dst_rect.right(), dst_center.y())
            return start, end

        if dy >= 0.0:
            return QtCore.QPointF(src_center.x(), src_rect.bottom()), QtCore.QPointF(dst_center.x(), dst_rect.top())
        return QtCore.QPointF(src_center.x(), src_rect.top()), QtCore.QPointF(dst_center.x(), dst_rect.bottom())

    @classmethod
    def _orthogonal_path(cls, start, end, route_offset=0.0):
        """Build an H-V-H elbow (2 horizontal + 1 vertical segment) between
        left/right anchor points.

        The single vertical jog is placed at the horizontal midpoint of the two
        anchors, shifted by *route_offset* to separate parallel edges.  The path
        always exits and arrives horizontally so the arrowhead direction matches
        the side of the node being entered.
        """
        dx = end.x() - start.x()
        dy = end.y() - start.y()

        # Degenerate: already a straight horizontal or vertical line.
        if abs(dy) < 1e-9 or abs(dx) < 1e-9:
            return [start, end]

        # Vertical jog column: midpoint + lane offset, clamped inside the
        # start→end x-span to prevent U-bends from large route_offset values.
        x_mid = (start.x() + end.x()) * 0.5 + route_offset
        if dx > 0:
            x_mid = max(start.x(), min(x_mid, end.x()))
        else:
            x_mid = min(start.x(), max(x_mid, end.x()))

        return cls._simplify_path([
            start,
            QtCore.QPointF(x_mid, start.y()),
            QtCore.QPointF(x_mid, end.y()),
            end,
        ])

    @staticmethod
    def _simplify_path(path_points):
        if not path_points:
            return []

        simplified = [path_points[0]]
        for point in path_points[1:]:
            prev = simplified[-1]
            if abs(point.x() - prev.x()) < 1e-9 and abs(point.y() - prev.y()) < 1e-9:
                continue
            simplified.append(point)

        if len(simplified) <= 2:
            return simplified

        pruned = [simplified[0]]
        for idx in range(1, len(simplified) - 1):
            a = pruned[-1]
            b = simplified[idx]
            c = simplified[idx + 1]
            collinear_x = abs(a.x() - b.x()) < 1e-9 and abs(b.x() - c.x()) < 1e-9
            collinear_y = abs(a.y() - b.y()) < 1e-9 and abs(b.y() - c.y()) < 1e-9
            if collinear_x or collinear_y:
                continue
            pruned.append(b)
        pruned.append(simplified[-1])
        return pruned

    @staticmethod
    def _path_segments(path_points):
        segments = []
        for idx in range(len(path_points) - 1):
            a = path_points[idx]
            b = path_points[idx + 1]
            if QtCore.QLineF(a, b).length() >= 1e-9:
                segments.append((a, b))
        return segments

    @staticmethod
    def _label_segment_for_path(path_points):
        segments = NetworkTopologyView._path_segments(path_points)
        if not segments:
            return None
        return max(segments, key=lambda seg: QtCore.QLineF(seg[0], seg[1]).length())

    def _draw_edge(self, edge_key, src_component, dst_component, path_points, label_segment, own_segments, edge_connections, src_object_label, dst_object_label, is_feedback=False):
        if is_feedback:
            line_pen = QtGui.QPen(QtGui.QColor(140, 140, 140))
            line_pen.setWidth(2)
            line_pen.setStyle(QtCore.Qt.DashLine)
        else:
            line_pen = QtGui.QPen(QtGui.QColor(70, 70, 70))
            line_pen.setWidth(2)

        if len(path_points) < 2:
            return

        segments = self._path_segments(path_points)
        n_segments = len(segments)
        for segment_index, (seg_start, seg_end) in enumerate(segments):
            self.scene.addLine(QtCore.QLineF(seg_start, seg_end), line_pen)

            # Allow route editing by dragging internal orthogonal segments.
            if segment_index == 0 or segment_index == n_segments - 1:
                continue
            dx = abs(seg_end.x() - seg_start.x())
            dy = abs(seg_end.y() - seg_start.y())
            if dx < 1e-9 and dy < 1e-9:
                continue
            orientation = "vertical" if dx < dy else "horizontal"
            drag_item = _DraggableEdgeSegmentItem(
                QtCore.QLineF(seg_start, seg_end),
                edge_key,
                segment_index,
                orientation,
                self,
            )
            self.scene.addItem(drag_item)

        arrow_start = path_points[-2]
        arrow_end = path_points[-1]
        self._draw_arrowhead(arrow_start, arrow_end)

        # Edge labels are drawn in a separate pass after fitInView.

    def _draw_edge_label(self, start, end, edge_connections, src_object_label, dst_object_label, own_segments=None):
        if not self.show_connection_labels or not edge_connections:
            return

        display_labels = [
            (
                f"{conn['source_short']} "
                f"{self.CONNECTION_LABEL_SEPARATOR} "
                f"{conn['dest_short']}"
            )
            for conn in edge_connections
        ]

        display_text = self._build_smart_label_text(display_labels)
        tooltip_builder = lambda src=src_object_label, dst=dst_object_label, conns=edge_connections: self._build_tooltip_html(src, dst, conns)
        text_item = _EdgeLabelItem(display_text, tooltip_builder=tooltip_builder)
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
        sx = abs(self.view.transform().m11())
        sy = abs(self.view.transform().m22())
        if sx < 1e-9:
            sx = 1.0
        if sy < 1e-9:
            sy = 1.0
        scale_mean = 0.5 * (sx + sy)
        if scale_mean < 1e-9:
            scale_mean = 1.0

        # Keep initial gap roughly constant in device pixels.
        normal_offset_px = 12.0
        normal_offset_scene = normal_offset_px / scale_mean
        center = QtCore.QPointF(mid.x() + nx * normal_offset_scene, mid.y() + ny * normal_offset_scene)

        text_rect = text_item.boundingRect()
        scene_text_rect = QtCore.QRectF(0.0, 0.0, text_rect.width() / sx, text_rect.height() / sy)
        placed_rect = self._place_label_rect(center, scene_text_rect, nx, ny, own_segments=own_segments)
        text_item.setPos(placed_rect.left(), placed_rect.top())

        placed_center = placed_rect.center()
        tx = ny
        ty = -nx
        signed_scene_offset = (placed_center.x() - mid.x()) * nx + (placed_center.y() - mid.y()) * ny
        signed_scene_tangent = (placed_center.x() - mid.x()) * tx + (placed_center.y() - mid.y()) * ty
        signed_offset_px = signed_scene_offset * scale_mean
        signed_tangent_px = signed_scene_tangent * scale_mean
        self._register_edge_leader(
            text_item,
            start,
            end,
            signed_offset_px=signed_offset_px,
            signed_tangent_px=signed_tangent_px,
        )
        self._label_occupied_rects.append(self._expand_rect(placed_rect, 4.0))

    def _register_edge_leader(self, text_item, line_start, line_end, signed_offset_px=12.0, signed_tangent_px=0.0):
        pen = QtGui.QPen(QtGui.QColor(150, 150, 150, 210))
        pen.setWidth(1)
        leader_item = self.scene.addLine(QtCore.QLineF(line_start, line_start), pen)
        self._edge_label_links.append(
            {
                "text_item": text_item,
                "leader_item": leader_item,
                "line_start": QtCore.QPointF(line_start.x(), line_start.y()),
                "line_end": QtCore.QPointF(line_end.x(), line_end.y()),
                "signed_offset_px": float(signed_offset_px),
                "signed_tangent_px": float(signed_tangent_px),
            }
        )

    def _update_edge_leaders(self):
        for link in self._edge_label_links:
            self._update_one_edge_leader(link)

    def _update_one_edge_leader(self, link):
        self._reposition_edge_label_for_zoom(link)

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

    def _reposition_edge_label_for_zoom(self, link):
        text_item = link["text_item"]
        line_start = link["line_start"]
        line_end = link["line_end"]

        line = QtCore.QLineF(line_start, line_end)
        length = line.length()
        if length < 1e-9:
            return

        sx = abs(self.view.transform().m11())
        sy = abs(self.view.transform().m22())
        if sx < 1e-9:
            sx = 1.0
        if sy < 1e-9:
            sy = 1.0

        scale_mean = 0.5 * (sx + sy)
        if scale_mean < 1e-9:
            scale_mean = 1.0

        vx = line_end.x() - line_start.x()
        vy = line_end.y() - line_start.y()
        nx = -vy / length
        ny = vx / length
        tx = ny
        ty = -nx
        signed_scene_offset = float(link.get("signed_offset_px", 12.0)) / scale_mean
        signed_scene_tangent = float(link.get("signed_tangent_px", 0.0)) / scale_mean

        mid = QtCore.QPointF((line_start.x() + line_end.x()) * 0.5, (line_start.y() + line_end.y()) * 0.5)
        label_center = QtCore.QPointF(
            mid.x() + nx * signed_scene_offset + tx * signed_scene_tangent,
            mid.y() + ny * signed_scene_offset + ty * signed_scene_tangent,
        )

        text_rect = text_item.boundingRect()
        scene_w = text_rect.width() / sx
        scene_h = text_rect.height() / sy
        text_item.setPos(label_center.x() - 0.5 * scene_w, label_center.y() - 0.5 * scene_h)

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
        if len(edge_labels) <= max_lines:
            return "\n".join(edge_labels)

        kept = list(edge_labels[:max_lines])
        omitted = len(edge_labels) - max_lines
        if omitted >= 2:
            kept.append(f"... (+{omitted} more)")
        else:
            kept.append(edge_labels[max_lines])
        return "\n".join(kept)

    @staticmethod
    def _scalar_from_value(value):
        """Return a scalar from scalar/array-like input; for arrays use the last element."""
        if value is None:
            return None

        arr = np.asarray(value)
        if arr.size == 0:
            return None
        last = arr.reshape(-1)[-1]
        return last.item() if hasattr(last, "item") else last

    @classmethod
    def _format_rel_err_text(cls, rel_err):
        rel_err_scalar = cls._scalar_from_value(rel_err)
        if rel_err_scalar is None:
            return "-"
        try:
            rel_err_float = float(rel_err_scalar)
        except Exception:
            return "-"
        if not np.isfinite(rel_err_float):
            return "-"
        return f"{rel_err_float:.2e}"

    @staticmethod
    def _build_tooltip_html(src_object_label, dst_object_label, full_connections):
        src_header = escape(str(src_object_label))
        dst_header = escape(str(dst_object_label))

        rows = []
        for conn in full_connections:
            source_short = conn["source_short"]
            dest_short = conn["dest_short"]

            connection = conn["connection"]
            if connection.has_step_history:
                status = connection.last_step_is_converged
                n_iter = connection.last_step_n_iter
                rel_err = connection.last_step_err_rel
            else:
                status = None
                current_iter = getattr(connection, "n_iter", 0)
                n_iter = current_iter if current_iter > 0 else None
                rel_err = connection.err_rel if current_iter > 0 else None

            status = NetworkTopologyView._scalar_from_value(status)
            n_iter = NetworkTopologyView._scalar_from_value(n_iter)

            if status is None:
                status_text = "n/a"
            else:
                status_text = "yes" if bool(status) else "no"

            n_iter_text = "-" if n_iter is None else str(n_iter)
            rel_err_text = NetworkTopologyView._format_rel_err_text(rel_err)

            solve_group = connection.solve_group
            group_text = "-" if solve_group is None else str(solve_group)

            rows.append(
                "<tr>"
                f"<td style='padding:2px 8px 2px 0;'>{escape(source_short)}</td>"
                f"<td style='padding:2px 0 2px 8px;'>{escape(dest_short)}</td>"
                f"<td style='padding:2px 0 2px 8px; white-space:nowrap;'>{escape(status_text)}</td>"
                f"<td style='padding:2px 0 2px 8px; white-space:nowrap; text-align:right;'>{escape(n_iter_text)}</td>"
                f"<td style='padding:2px 0 2px 8px; white-space:nowrap; text-align:right; font-family:Consolas, \"Courier New\", monospace;'>{escape(rel_err_text)}</td>"
                f"<td style='padding:2px 0 2px 8px;'>{escape(group_text)}</td>"
                "</tr>"
            )
        rows_html = "".join(rows)

        return (
            "<html>"
            "<b><span style='color:#cc0000;'>Connections</span></b>"
            "<table style='margin-top:4px; border-collapse:collapse;'>"
            "<thead><tr>"
            f"<th style='text-align:left; padding:0 8px 2px 0;'>{src_header}</th>"
            f"<th style='text-align:left; padding:0 0 2px 8px;'>{dst_header}</th>"
            "<th style='text-align:left; padding:0 0 2px 8px; white-space:nowrap;'>Converged</th>"
            "<th style='text-align:right; padding:0 0 2px 8px; white-space:nowrap;'>n_iter</th>"
            "<th style='text-align:right; padding:0 0 2px 8px; white-space:nowrap;'>rel_err</th>"
            "<th style='text-align:left; padding:0 0 2px 8px;'>group</th>"
            "</tr></thead>"
            f"<tbody>{rows_html}</tbody>"
            "</table>"
            "</html>"
        )

    def _place_label_rect(self, center, text_rect, nx, ny, own_segments=None):
        sx = abs(self.view.transform().m11())
        sy = abs(self.view.transform().m22())
        if sx < 1e-9:
            sx = 1.0
        if sy < 1e-9:
            sy = 1.0
        scale_mean = 0.5 * (sx + sy)
        if scale_mean < 1e-9:
            scale_mean = 1.0

        # Search offsets in pixel-sized increments for zoom-stable spacing.
        base_offset = 12.0 / scale_mean
        step = 12.0 / scale_mean
        tangent_step = 16.0 / scale_mean
        max_tries = 12
        tangent_tries = 6

        best_rect = None
        best_penalty = float("inf")

        normal_candidates = [0.0]
        for i in range(1, max_tries + 1):
            delta = i * step
            normal_candidates.append(delta)
            normal_candidates.append(-delta)

        tangent_candidates = [0.0]
        for i in range(1, tangent_tries + 1):
            shift = i * tangent_step
            tangent_candidates.append(shift)
            tangent_candidates.append(-shift)

        tx = ny
        ty = -nx
        candidate_pairs = [(delta_n, delta_t) for delta_t in tangent_candidates for delta_n in normal_candidates]
        candidate_pairs.sort(key=lambda pair: (abs(pair[0]) + 0.75 * abs(pair[1]), abs(pair[1])))

        min_node_clearance = 10.0 / scale_mean
        for delta_n, delta_t in candidate_pairs:
            offset = base_offset + delta_n
            cx = center.x() + nx * offset + tx * delta_t
            cy = center.y() + ny * offset + ty * delta_t
            rect = QtCore.QRectF(
                cx - text_rect.width() / 2,
                cy - text_rect.height() / 2,
                text_rect.width(),
                text_rect.height(),
            )
            penalty = self._placement_penalty(
                rect,
                own_segments=own_segments,
                min_node_clearance=min_node_clearance,
            )
            if penalty < best_penalty:
                best_penalty = penalty
                best_rect = rect
                if penalty == 0.0:
                    break

        return best_rect

    def _placement_penalty(self, rect, own_segments=None, min_node_clearance=0.0):
        padded_rect = self._expand_rect(rect, 2.0)
        penalty = 0.0

        for node_rect in self._node_rects.values():
            if padded_rect.intersects(node_rect):
                # Treat node overlap as expensive so labels prefer nearby whitespace.
                penalty += 20000.0 + 20.0 * self._intersection_area(padded_rect, node_rect)
            elif min_node_clearance > 0.0:
                gap = self._rect_separation(padded_rect, node_rect)
                if gap < min_node_clearance:
                    penalty += 150.0 * (min_node_clearance - gap)

        for used_rect in self._label_occupied_rects:
            if padded_rect.intersects(used_rect):
                penalty += self._intersection_area(padded_rect, used_rect) + 500.0

        for seg_start, seg_end in self._connector_segments:
            if own_segments is not None and any(
                self._segments_match(seg_start, seg_end, own_seg[0], own_seg[1])
                for own_seg in own_segments
            ):
                continue
            if self._segment_intersects_rect(seg_start, seg_end, padded_rect):
                penalty += self.LABEL_LINE_OVERLAP_PENALTY

        return penalty

    @staticmethod
    def _rect_separation(a, b):
        dx = max(0.0, max(a.left() - b.right(), b.left() - a.right()))
        dy = max(0.0, max(a.top() - b.bottom(), b.top() - a.bottom()))
        return math.hypot(dx, dy)

    @staticmethod
    def _segments_match(a1, a2, b1, b2, tol=1e-9):
        def _close(p, q):
            return abs(p.x() - q.x()) <= tol and abs(p.y() - q.y()) <= tol

        return (_close(a1, b1) and _close(a2, b2)) or (_close(a1, b2) and _close(a2, b1))

    @staticmethod
    def _segment_intersects_rect(seg_start, seg_end, rect):
        if rect.contains(seg_start) or rect.contains(seg_end):
            return True

        p1 = QtCore.QPointF(rect.left(), rect.top())
        p2 = QtCore.QPointF(rect.right(), rect.top())
        p3 = QtCore.QPointF(rect.right(), rect.bottom())
        p4 = QtCore.QPointF(rect.left(), rect.bottom())
        rect_edges = ((p1, p2), (p2, p3), (p3, p4), (p4, p1))

        for e1, e2 in rect_edges:
            if NetworkTopologyView._segments_intersect(seg_start, seg_end, e1, e2):
                return True
        return False

    @staticmethod
    def _segments_intersect(a, b, c, d):
        def orientation(p, q, r):
            val = (q.y() - p.y()) * (r.x() - q.x()) - (q.x() - p.x()) * (r.y() - q.y())
            if abs(val) < 1e-12:
                return 0
            return 1 if val > 0 else 2

        def on_segment(p, q, r):
            return (
                min(p.x(), r.x()) - 1e-12 <= q.x() <= max(p.x(), r.x()) + 1e-12
                and min(p.y(), r.y()) - 1e-12 <= q.y() <= max(p.y(), r.y()) + 1e-12
            )

        o1 = orientation(a, b, c)
        o2 = orientation(a, b, d)
        o3 = orientation(c, d, a)
        o4 = orientation(c, d, b)

        if o1 != o2 and o3 != o4:
            return True

        if o1 == 0 and on_segment(a, c, b):
            return True
        if o2 == 0 and on_segment(a, d, b):
            return True
        if o3 == 0 and on_segment(c, a, d):
            return True
        if o4 == 0 and on_segment(c, b, d):
            return True

        return False

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

    def _on_node_drag_finished(self, comp, new_center):
        """Record the dragged position and schedule a full redraw."""
        if self._snap_to_grid_enabled:
            new_center = self._snap_point_to_grid(new_center)
        layout_key = self._component_layout_keys.get(comp)
        if layout_key is None:
            layout_key = str(getattr(comp, "name", "") or "").strip() or type(comp).__name__
        self._override_positions[layout_key] = (new_center.x(), new_center.y())
        QtCore.QTimer.singleShot(0, self._draw_topology)

    def _on_edge_segment_drag_finished(self, edge_key, segment_index, delta):
        path_points = self._rendered_edge_paths.get(edge_key)
        if path_points is None or len(path_points) < 4:
            QtCore.QTimer.singleShot(0, self._draw_topology)
            return
        if segment_index <= 0 or segment_index >= len(path_points) - 2:
            QtCore.QTimer.singleShot(0, self._draw_topology)
            return

        path_points = self._copy_path_points(path_points)
        p0 = path_points[segment_index]
        p1 = path_points[segment_index + 1]
        dx = abs(p1.x() - p0.x())
        dy = abs(p1.y() - p0.y())

        if dx < dy:
            new_x = p0.x() + float(delta)
            if self._snap_to_grid_enabled:
                new_x = self._snap_point_to_grid(QtCore.QPointF(new_x, 0.0)).x()
            path_points[segment_index].setX(new_x)
            path_points[segment_index + 1].setX(new_x)
        else:
            new_y = p0.y() + float(delta)
            if self._snap_to_grid_enabled:
                new_y = self._snap_point_to_grid(QtCore.QPointF(0.0, new_y)).y()
            path_points[segment_index].setY(new_y)
            path_points[segment_index + 1].setY(new_y)

        path_points = self._simplify_path(path_points)
        if len(path_points) >= 2:
            self._edge_path_overrides[edge_key] = path_points

        QtCore.QTimer.singleShot(0, self._draw_topology)

    def _toggle_snap_to_grid(self, checked):
        self._snap_to_grid_enabled = bool(checked)
        if self._snap_grid_spinbox is not None:
            self._snap_grid_spinbox.setEnabled(self._snap_to_grid_enabled)

    def _set_snap_grid_size(self, value):
        self._snap_grid_size_px = max(1.0, float(value))

    def _snap_point_to_grid(self, point):
        step = float(self._snap_grid_size_px)
        if step <= 0.0:
            return QtCore.QPointF(point.x(), point.y())
        snapped_x = round(point.x() / step) * step
        snapped_y = round(point.y() / step) * step
        return QtCore.QPointF(snapped_x, snapped_y)

    def save_layout(self, path):
        """Save current node positions to a JSON layout file."""
        positions_data = {}
        for comp, rect in self._node_rects.items():
            layout_key = self._component_layout_keys.get(comp)
            if layout_key is None:
                continue
            center = rect.center()
            positions_data[layout_key] = {"x": center.x(), "y": center.y()}

        edge_paths_data = {
            edge_key: self._path_points_to_json(path_points)
            for edge_key, path_points in self._edge_path_overrides.items()
            if path_points and len(path_points) >= 2
        }

        data = {
            "version": 1,
            "positions": positions_data,
            "edge_paths": edge_paths_data,
        }
        with open(path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=2)

    def load_layout(self, path):
        """Load node positions from a JSON layout file and redraw."""
        self._load_layout_file(path)
        self._draw_topology()

    def _load_layout_file(self, path):
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        version = data.get("version")
        if version != 1:
            raise ValueError(f"Unsupported layout file version: {version!r}")
        self._override_positions = {
            name: (float(d["x"]), float(d["y"]))
            for name, d in data.get("positions", {}).items()
        }

        edge_path_overrides = {}
        for edge_key, points_json in data.get("edge_paths", {}).items():
            points = self._path_points_from_json(points_json)
            points = self._simplify_path(points)
            if len(points) >= 2:
                edge_path_overrides[edge_key] = points
        self._edge_path_overrides = edge_path_overrides

    def _toolbar_save(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self.view, "Save Layout", "", "JSON Files (*.json)"
        )
        if path:
            self.save_layout(path)

    def _toolbar_load(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self.view, "Load Layout", "", "JSON Files (*.json)"
        )
        if path:
            self.load_layout(path)

    def _toolbar_reset(self):
        self._override_positions.clear()
        self._edge_path_overrides.clear()
        self._draw_topology()

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
