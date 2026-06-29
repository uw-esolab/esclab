import numpy as np
import pyqtgraph as qtg
from pyqtgraph.Qt import QtCore, QtWidgets, QtGui
import time


class _AxisLockedViewBox(qtg.ViewBox):
    """ViewBox with modifier-key-constrained scroll zoom.
    Ctrl+scroll  -> x-axis only
    Shift+scroll -> y-axis only
    plain scroll -> both axes (default)
    """

    def wheelEvent(self, ev, axis=None):
        mods = ev.modifiers()
        if mods & QtCore.Qt.ControlModifier:
            axis = 0
        elif mods & QtCore.Qt.ShiftModifier:
            axis = 1
        super().wheelEvent(ev, axis=axis)


class _PlotterKeyEventFilter(QtCore.QObject):
    """Global key filter for plotter-window keyboard navigation."""

    def eventFilter(self, obj, event):
        if event.type() == QtCore.QEvent.KeyPress and OnlinePlotter.handle_keyboard_event(event):
            return True
        return super().eventFilter(obj, event)


class OnlinePlotter:
    app = None
    main_window = None
    tab_widget = None
    follow_button = None
    tooltip_button = None
    instances = []
    font_targets = []
    tab_shortcuts = []
    key_event_filter = None
    KEYBOARD_ZOOM_IN_FACTOR = 1.15
    KEYBOARD_ZOOM_OUT_FACTOR = 1.0 / KEYBOARD_ZOOM_IN_FACTOR
    KEYBOARD_PAN_PIXELS = 80
    n_plotters = 0
    font_size_pt = 10
    min_font_size_pt = 6
    max_font_size_pt = 30
    UI_EVENT_PROCESS_INTERVAL_SEC = 0.02
    UI_EVENT_PROCESS_BUDGET_MS = 3
    _last_ui_event_process = 0.0

    @classmethod
    def ensure_window(cls, plotter_size=(0.9, 0.9)):
        if cls.app is None:
            cls.app = qtg.mkQApp()

        if cls.main_window is not None:
            return cls.main_window

        main_window = QtWidgets.QMainWindow()
        main_window.setWindowTitle("ᴇꜱᴄʟᴀʙ ☀")
        tab_widget = QtWidgets.QTabWidget()

        controls_widget = QtWidgets.QWidget()
        controls_layout = QtWidgets.QHBoxLayout(controls_widget)
        controls_layout.setContentsMargins(6, 6, 6, 0)
        controls_layout.setSpacing(6)

        font_up_button = QtWidgets.QPushButton("🗚")
        font_down_button = QtWidgets.QPushButton("🗛")
        font_up_button.setFixedSize(30, 30)
        font_down_button.setFixedSize(30, 30)
        font_up_button.setToolTip("Increase font size")
        font_down_button.setToolTip("Decrease font size")
        controls_layout.addWidget(font_up_button)
        controls_layout.addWidget(font_down_button)
        follow_button = QtWidgets.QPushButton("Follow")
        follow_button.setCheckable(True)
        follow_button.setChecked(True)
        follow_button.setFixedSize(80, 30)
        follow_button.setToolTip("Auto-follow the latest data in the plot (data plotter only)")
        controls_layout.addWidget(follow_button)
        tooltip_button = QtWidgets.QPushButton("Tooltip")
        tooltip_button.setCheckable(True)
        tooltip_button.setChecked(False)
        tooltip_button.setFixedSize(70, 30)
        tooltip_button.setToolTip("Show tooltip with values at the cursor position (data plotter only)")
        controls_layout.addWidget(tooltip_button)
        
        save_png_button = QtWidgets.QPushButton("Save Image")
        save_png_button.setFixedSize(80, 30)
        save_png_button.setToolTip("Save the current plot as a PNG image")
        controls_layout.addWidget(save_png_button)
        
        controls_layout.addStretch()

        container = QtWidgets.QWidget()
        container_layout = QtWidgets.QVBoxLayout(container)
        container_layout.setContentsMargins(0, 0, 0, 0)
        container_layout.setSpacing(0)
        container_layout.addWidget(controls_widget)
        container_layout.addWidget(tab_widget)
        main_window.setCentralWidget(container)

        font_down_button.clicked.connect(lambda: cls.adjust_font_size(-1))
        font_up_button.clicked.connect(lambda: cls.adjust_font_size(1))
        follow_button.clicked.connect(lambda checked: cls.set_auto_follow(checked))
        tooltip_button.clicked.connect(lambda checked: cls.set_plot_tooltip(checked))
        save_png_button.clicked.connect(lambda: cls.save_current_tab_png())
        next_tab_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("Ctrl+Tab"), main_window)
        prev_tab_shortcut = QtWidgets.QShortcut(QtGui.QKeySequence("Ctrl+Shift+Tab"), main_window)
        next_tab_shortcut.setContext(QtCore.Qt.WindowShortcut)
        prev_tab_shortcut.setContext(QtCore.Qt.WindowShortcut)
        next_tab_shortcut.activated.connect(lambda: cls.cycle_tabs(forward=True))
        prev_tab_shortcut.activated.connect(lambda: cls.cycle_tabs(forward=False))

        all_shortcuts = [next_tab_shortcut, prev_tab_shortcut]
        for shortcut in all_shortcuts:
            shortcut.setContext(QtCore.Qt.WindowShortcut)

        cls.tab_shortcuts = all_shortcuts
        cls.follow_button = follow_button
        cls.tooltip_button = tooltip_button

        if cls.key_event_filter is None:
            cls.key_event_filter = _PlotterKeyEventFilter()
            cls.app.installEventFilter(cls.key_event_filter)

        screen = QtWidgets.QApplication.primaryScreen()
        if screen is not None:
            screen_rect = screen.availableGeometry()
            if plotter_size is not None:
                if plotter_size[0] <= 1.0 and plotter_size[1] <= 1.0:
                    main_window.resize(int(screen_rect.width() * plotter_size[0]), int(screen_rect.height() * plotter_size[1]))
                else:
                    main_window.resize(plotter_size[0], plotter_size[1])
            frame = main_window.frameGeometry()
            frame.moveCenter(screen_rect.center())
            main_window.move(frame.topLeft())
        elif plotter_size is not None:
            if plotter_size[0] <= 1.0 and plotter_size[1] <= 1.0:
                main_window.resize(1000, 600)
            else:
                main_window.resize(plotter_size[0], plotter_size[1])

        main_window.show()
        cls.main_window = main_window
        cls.tab_widget = tab_widget
        return main_window

    @classmethod
    def cycle_tabs(cls, forward=True):
        if cls.tab_widget is None:
            return

        count = cls.tab_widget.count()
        if count <= 1:
            return

        current = cls.tab_widget.currentIndex()
        if current < 0:
            current = 0

        step = 1 if forward else -1
        cls.tab_widget.setCurrentIndex((current + step) % count)

    @classmethod
    def save_current_tab_png(cls, dpi=300):
        widget = cls._current_tab_widget()
        if widget is None:
            return
        
        path, _ = QtWidgets.QFileDialog.getSaveFileName(cls.main_window, "Save Image", "", "PNG Files (*.png)")
        if not path:
            return
            
        if not path.lower().endswith('.png'):
            path += '.png'
            
        screen_dpi = getattr(cls.main_window, 'logicalDpiX', lambda: 96.0)()
        scale_factor = dpi / screen_dpi

        orig_size = widget.size()
        scaled_width = int(orig_size.width() * scale_factor)
        scaled_height = int(orig_size.height() * scale_factor)
        scaled_size = QtCore.QSize(scaled_width, scaled_height)
        
        image = QtGui.QImage(scaled_size, QtGui.QImage.Format_ARGB32)
        image.fill(QtGui.QColor("white"))
        
        painter = QtGui.QPainter(image)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        
        if hasattr(widget, 'scene') and widget.scene() is not None:
            target_rect = QtCore.QRectF(0, 0, scaled_size.width(), scaled_size.height())
            source_rect = QtCore.QRectF(0, 0, orig_size.width(), orig_size.height())
            
            # Temporarily disable ItemIgnoresTransformations for all items
            # (such as pyqtgraph's LegendItem) so they scale up correctly
            flags_to_restore = []
            export_items = []
            
            for item in widget.scene().items():
                if int(item.flags()) & QtWidgets.QGraphicsItem.ItemIgnoresTransformations:
                    item.setFlag(QtWidgets.QGraphicsItem.ItemIgnoresTransformations, False)
                    flags_to_restore.append(item)
                
                # Notify pyqtgraph items so cosmetic pens scale up natively
                if hasattr(item, 'setExportMode'):
                    item.setExportMode(True, {'resolutionScale': scale_factor, 'antialias': True})
                    export_items.append(item)
                    
            widget.scene().render(painter, target_rect, source_rect)
            
            # Restore the flags and export mode
            for item in flags_to_restore:
                item.setFlag(QtWidgets.QGraphicsItem.ItemIgnoresTransformations, True)
            for item in export_items:
                item.setExportMode(False)
        else:
            painter.scale(scale_factor, scale_factor)
            widget.render(painter)
            
        painter.end()
        
        # Set DPI metadata for the image file
        dpm = int(dpi * 100 / 2.54)
        image.setDotsPerMeterX(dpm)
        image.setDotsPerMeterY(dpm)
        
        image.save(path, "PNG")

    @classmethod
    def handle_keyboard_event(cls, event):
        """Handle Ctrl+zoom/pan keys for network tabs only."""
        if cls.main_window is None or cls.tab_widget is None:
            return False
        if not cls.main_window.isActiveWindow():
            return False

        mods = event.modifiers()
        if not (mods & QtCore.Qt.ControlModifier):
            return False
        if mods & QtCore.Qt.AltModifier or mods & QtCore.Qt.MetaModifier:
            return False

        key = event.key()
        if key in (QtCore.Qt.Key_Tab, QtCore.Qt.Key_Backtab):
            return False

        if key in (QtCore.Qt.Key_Minus, QtCore.Qt.Key_Underscore):
            return cls.zoom_current_network_tab(zoom_in=False)
        if key in (QtCore.Qt.Key_Plus, QtCore.Qt.Key_Equal):
            return cls.zoom_current_network_tab(zoom_in=True)
        if key == QtCore.Qt.Key_Up:
            return cls.pan_current_network_tab("up")
        if key == QtCore.Qt.Key_Down:
            return cls.pan_current_network_tab("down")
        if key == QtCore.Qt.Key_Left:
            return cls.pan_current_network_tab("left")
        if key == QtCore.Qt.Key_Right:
            return cls.pan_current_network_tab("right")

        return False

    @classmethod
    def _current_tab_widget(cls):
        if cls.tab_widget is None:
            return None
        return cls.tab_widget.currentWidget()

    @classmethod
    def _current_network_view(cls):
        widget = cls._current_tab_widget()
        if isinstance(widget, QtWidgets.QGraphicsView):
            return widget
        return None

    @classmethod
    def zoom_current_network_tab(cls, zoom_in=True):
        view = cls._current_network_view()
        if view is None:
            return False

        factor = cls.KEYBOARD_ZOOM_IN_FACTOR if zoom_in else cls.KEYBOARD_ZOOM_OUT_FACTOR
        current_scale = view.transform().m11()
        if current_scale <= 0.0:
            current_scale = 1.0
        min_scale = getattr(view, "MIN_SCALE", 0.02)
        max_scale = getattr(view, "MAX_SCALE", 50.0)
        target_scale = max(min_scale, min(max_scale, current_scale * factor))
        scale_factor = target_scale / current_scale
        view.scale(scale_factor, scale_factor)
        callback = getattr(view, "transform_changed_callback", None)
        if callback is not None:
            callback()
        return True

    @classmethod
    def pan_current_network_tab(cls, direction):
        view = cls._current_network_view()
        if view is None:
            return False

        step = int(cls.KEYBOARD_PAN_PIXELS)
        hbar = view.horizontalScrollBar()
        vbar = view.verticalScrollBar()
        if direction == "left":
            hbar.setValue(hbar.value() - step)
        elif direction == "right":
            hbar.setValue(hbar.value() + step)
        elif direction == "up":
            vbar.setValue(vbar.value() - step)
        elif direction == "down":
            vbar.setValue(vbar.value() + step)
        else:
            return False

        return True

    """
    y1 | [a, b, c, ...] component input or output values to be plotted on axis 1
    y2 | [d, e, f, ...] component input or output values to be plotted on axis 2
    """

    def __init__(self, y1, y2, y1lim, y2lim, y1label, y2label, nmax_points, update_every, plotter_size=(0.9, 0.9), tab_title=None, show_live=True):
        assert isinstance(y1, type([]))
        assert isinstance(y2, type([])) or y2 is None

        self.show_live = show_live
        self.current_step = -1
        self._fill_idx = 0
        self._auto_follow = True
        self.y1label = y1label
        self.y2label = y2label
        self.nmax_points = nmax_points
        self.update_every = update_every
        self.y1_items = y1
        self.y2_items = y2
        self.x_data = None
        self.y1_data = None
        self.y2_data = None
        self._y1lim = y1lim
        self._y2lim = y2lim
        self._plotter_size = plotter_size
        self._tab_title = tab_title

        if show_live:
            self._build_widgets()
        else:
            self.y1_lines = []
            self.y2_lines = []

    def _build_widgets(self):
        """Create all Qt plot widgets and register this plotter in the shared window."""
        colors = [
            "#377eb8",
            "#ff7f00",
            "#4daf4a",
            "#f781bf",
            "#a65628",
            "#984ea3",
            "#999999",
            "#e41a1c",
            "#dede00",
            "#00ffe5",
            "#5d0ab6",
        ]
        OnlinePlotter.ensure_window(plotter_size=self._plotter_size)
        self.app = OnlinePlotter.app

        self.win = qtg.GraphicsLayoutWidget()
        OnlinePlotter.n_plotters += 1
        tab_label = self._tab_title if self._tab_title not in [None, ""] else self.y1label if self.y1label not in [None, ""] else f"Plot {OnlinePlotter.n_plotters}"
        OnlinePlotter.tab_widget.addTab(self.win, tab_label)

        self.ax1 = self.win.addPlot(viewBox=_AxisLockedViewBox())
        self.ax1.setLabel("bottom", "Time")
        self.ax1.setLabel("left", self.y1label)
        self.legend_y1 = self.ax1.addLegend(offset=(10, 10))
        self.legend_y1.setBrush(QtGui.QBrush(QtGui.QColor(255, 255, 255, 40)))
        self.legend_y2 = None
        self.ax1.showGrid(x=True, y=True, alpha=0.3)
        self._tooltip_vline = qtg.InfiniteLine(
            angle=90,
            movable=False,
            pen=qtg.mkPen(QtGui.QColor(220, 40, 40, 220), width=1.5),
        )
        self._tooltip_vline.setZValue(1000)
        self.ax1.addItem(self._tooltip_vline, ignoreBounds=True)
        self._tooltip_vline.hide()

        if self._y1lim is not None:
            self.ax1.setYRange(*self._y1lim)

        self.y1_lines = []
        c = 0
        for i in range(len(self.y1_items)):
            pen = qtg.mkPen(color=colors[c % len(colors)], width=2)
            line = self.ax1.plot([], [], pen=pen, name=self.y1_items[i].name)
            self.y1_lines.append(line)
            c += 1

        self.y2_lines = []
        if self.y2_items is not None:
            self.ax2 = _AxisLockedViewBox()
            self.ax1.showAxis("right")
            self.ax1.scene().addItem(self.ax2)
            self.ax1.getAxis("right").linkToView(self.ax2)
            self.ax2.setXLink(self.ax1)
            self.ax1.getAxis("right").setLabel(self.y2label, color="#c0c0c0")
            self.legend_y2 = qtg.LegendItem(offset=(-10, 10))
            self.legend_y2.setBrush(QtGui.QBrush(QtGui.QColor(255, 255, 255, 40)))
            self.legend_y2.setParentItem(self.ax1.vb)

            if self._y2lim is not None:
                self.ax2.setYRange(*self._y2lim)

            def updateViews():
                self.ax2.setGeometry(self.ax1.vb.sceneBoundingRect())
                self.ax2.linkedViewChanged(self.ax1.vb, self.ax2.XAxis)

            updateViews()
            self.ax1.vb.sigResized.connect(updateViews)

            for i in range(len(self.y2_items)):
                pen = qtg.mkPen(color=colors[c % len(colors)], width=2, style=QtCore.Qt.DashLine)
                line = qtg.PlotDataItem([], [], pen=pen, name=self.y2_items[i].name)
                self.ax2.addItem(line)
                self.legend_y2.addItem(line, self.y2_items[i].name)
                self.y2_lines.append(line)
                c += 1

        self.win.nextRow()
        self.ax_conv = self.win.addPlot()
        self.ax_conv.setMaximumHeight(18)
        self.ax_conv.setMinimumHeight(18)
        for _axis in ("left", "bottom", "right", "top"):
            self.ax_conv.hideAxis(_axis)
        self.ax_conv.hideButtons()
        self.ax_conv.setMouseEnabled(x=False, y=False)
        self.ax_conv.setYRange(0, 1, padding=0)
        self.ax_conv.setXLink(self.ax1)
        self._conv_bar_item = None

        self.win.nextRow()
        self.ax_iter = self.win.addPlot()
        self.ax_iter.setMaximumHeight(18)
        self.ax_iter.setMinimumHeight(18)
        for _axis in ("left", "bottom", "right", "top"):
            self.ax_iter.hideAxis(_axis)
        self.ax_iter.hideButtons()
        self.ax_iter.setMouseEnabled(x=False, y=False)
        self.ax_iter.setYRange(0, 1, padding=0)
        self.ax_iter.setXLink(self.ax1)
        self._iter_bar_item = None

        if self not in OnlinePlotter.instances:
            OnlinePlotter.instances.append(self)
        self.apply_font_size()
        self.ax1.vb.sigRangeChangedManually.connect(self._on_range_changed_manually)
        self._setup_strip_tooltip()

    def _finalize(self):
        """Build Qt widgets and render all buffered data. Called from wait_for_plots() for deferred plotters."""
        self._build_widgets()
        self.__refresh_plot()

    @classmethod
    def adjust_font_size(cls, delta):
        new_font_size = max(cls.min_font_size_pt, min(cls.max_font_size_pt, cls.font_size_pt + delta))
        if new_font_size == cls.font_size_pt:
            return

        cls.font_size_pt = new_font_size
        for plotter in cls.instances:
            plotter.apply_font_size()
        for target in cls.font_targets:
            target.apply_font_size()

    @classmethod
    def register_font_target(cls, target):
        if target not in cls.font_targets:
            cls.font_targets.append(target)

    def apply_font_size(self):
        tick_font = QtGui.QFont()
        tick_font.setPointSize(OnlinePlotter.font_size_pt)
        label_color = "#c0c0c0"
        label_style = {"font-size": f"{OnlinePlotter.font_size_pt + 2}pt"}

        left_axis = self.ax1.getAxis("left")
        bottom_axis = self.ax1.getAxis("bottom")
        left_axis.setStyle(tickFont=tick_font)
        bottom_axis.setStyle(tickFont=tick_font)

        self.ax1.setLabel("left", self.y1label, color=label_color, **label_style)
        self.ax1.setLabel("bottom", "Time", color=label_color, **label_style)

        if self.legend_y1 is not None:
            label_font = self.legend_y1.font()
            label_font.setPointSize(OnlinePlotter.font_size_pt)
            self.legend_y1.setFont(label_font)

            for _, label_item in self.legend_y1.items:
                label_item.setText(label_item.text, size=f"{OnlinePlotter.font_size_pt}pt")

        if self.y2_items is not None:
            right_axis = self.ax1.getAxis("right")
            right_axis.setStyle(tickFont=tick_font)
            right_axis.setLabel(self.y2label, color=label_color, **label_style)

            if self.legend_y2 is not None:
                for _, label_item in self.legend_y2.items:
                    label_item.setText(label_item.text, size=f"{OnlinePlotter.font_size_pt}pt")

    @staticmethod
    def _fraction_to_brush(f):
        if f <= 0.0:
            return QtGui.QBrush(QtGui.QColor(0, 0, 0, 0))
        if f <= 0.5:
            t = f * 2.0
            r, g, b = 255, int(255 - t * (255 - 165)), 0
        else:
            t = (f - 0.5) * 2.0
            r, g, b = int(255 - t * (255 - 200)), int(165 * (1.0 - t)), 0
        return QtGui.QBrush(QtGui.QColor(r, g, b, 220))

    @staticmethod
    def _iter_fraction_to_brush(f):
        if f <= 0.0:
            return QtGui.QBrush(QtGui.QColor(0, 0, 0, 0))
        alpha = int(f * 220)
        r = int(255 * (1.0 - f))
        g = int(255 * (1.0 - f))
        b = 220
        return QtGui.QBrush(QtGui.QColor(r, g, b, alpha))

    def log_step(self, time, conv_fraction=0.0, iter_fraction=0.0):
        self.current_step += 1
        idx = self._fill_idx

        self.x_data[idx] = time
        self.conv_data[idx] = conv_fraction
        self.iter_data[idx] = iter_fraction

        for j, yval in enumerate(self.y1_items):
            self.y1_data[j, idx] = yval.v

        if self.y2_items is not None:
            for j, yval in enumerate(self.y2_items):
                self.y2_data[j, idx] = yval.v

        self._fill_idx += 1

        if self.show_live and self.current_step % self.update_every == 0:
            self.__refresh_plot()

    def __refresh_plot(self):
        n = self._fill_idx
        if n == 0:
            return
        x_view = self.x_data[:n]

        for i in range(len(self.y1_lines)):
            self.y1_lines[i].setData(x_view, self.y1_data[i, :n])

        for i in range(len(self.y2_lines)):
            self.y2_lines[i].setData(x_view, self.y2_data[i, :n])

        if self._auto_follow:
            x_start = self.x_data[max(0, n - self.nmax_points)]
            x_end = self.x_data[n - 1]
            self.ax1.setXRange(x_start, x_end, padding=0.02)

        if self._conv_bar_item is not None:
            self.ax_conv.removeItem(self._conv_bar_item)
            self._conv_bar_item = None
        nc_mask = self.conv_data[:n] > 0
        if np.any(nc_mask):
            nc_x = x_view[nc_mask]
            brushes = [self._fraction_to_brush(f) for f in self.conv_data[:n][nc_mask]]
            self._conv_bar_item = qtg.BarGraphItem(
                x=nc_x, height=1, width=self._timestep, brushes=brushes, pen=qtg.mkPen(None)
            )
            self.ax_conv.addItem(self._conv_bar_item)

        if self._iter_bar_item is not None:
            self.ax_iter.removeItem(self._iter_bar_item)
            self._iter_bar_item = None
        iter_brushes = [self._iter_fraction_to_brush(f) for f in self.iter_data[:n]]
        self._iter_bar_item = qtg.BarGraphItem(
            x=x_view, height=1, width=self._timestep, brushes=iter_brushes, pen=qtg.mkPen(None)
        )
        self.ax_iter.addItem(self._iter_bar_item)

    def preallocate(self, n_steps, timestep):
        self._fill_idx = 0
        self._timestep = timestep
        self.x_data = np.empty(n_steps)
        self.y1_data = np.empty((len(self.y1_items), n_steps))
        n_y2 = len(self.y2_items) if self.y2_items is not None else 1
        self.y2_data = np.empty((n_y2, n_steps))
        self.conv_data = np.zeros(n_steps)
        self.iter_data = np.zeros(n_steps)

    def _setup_strip_tooltip(self):
        self._strip_proxy = qtg.SignalProxy(self.ax1.scene().sigMouseMoved, rateLimit=20, slot=self._on_strip_mouse_move)
        self._tooltip_last_idx = -1
        self._plot_tooltip_enabled = False
        self._tip_label = QtWidgets.QLabel()
        self._tip_label.setWindowFlags(QtCore.Qt.ToolTip | QtCore.Qt.FramelessWindowHint | QtCore.Qt.WindowStaysOnTopHint)
        self._tip_label.setStyleSheet(
            "QLabel { background-color: #ffffdc; color: #000000; border: 1px solid #aaaaaa; padding: 4px 8px; }"
        )
        self._tip_label.setAttribute(QtCore.Qt.WA_TransparentForMouseEvents)
        self._tip_label.hide()

    def _on_strip_mouse_move(self, args):
        scene_pos = args[0]

        if self.x_data is None or self._fill_idx == 0:
            self._tip_label.hide()
            self._tooltip_vline.hide()
            self._tooltip_last_idx = -1
            return

        in_conv = self.ax_conv.vb.sceneBoundingRect().contains(scene_pos)
        in_iter = self.ax_iter.vb.sceneBoundingRect().contains(scene_pos)
        in_plot = self.ax1.vb.sceneBoundingRect().contains(scene_pos)
        in_y2 = self.y2_items is not None and hasattr(self, "ax2") and self.ax2.sceneBoundingRect().contains(scene_pos)

        if not in_conv and not in_iter and not in_plot and not in_y2:
            self._tip_label.hide()
            self._tooltip_vline.hide()
            self._tooltip_last_idx = -1
            return

        if (in_plot or in_y2) and not in_conv and not in_iter and not self._plot_tooltip_enabled:
            self._tip_label.hide()
            self._tooltip_vline.hide()
            self._tooltip_last_idx = -1
            return

        vb = self.ax_conv.vb if in_conv else (self.ax_iter.vb if in_iter else self.ax1.vb)
        x_val = vb.mapSceneToView(scene_pos).x()

        n = self._fill_idx
        x_slice = self.x_data[:n]
        idx = int(np.searchsorted(x_slice, x_val, side="left"))
        idx = max(0, min(idx, n - 1))
        if idx > 0 and abs(x_slice[idx - 1] - x_val) < abs(x_slice[idx] - x_val):
            idx -= 1

        t = x_slice[idx]
        if self._plot_tooltip_enabled:
            self._tooltip_vline.setPos(t)
            self._tooltip_vline.show()
        else:
            self._tooltip_vline.hide()

        cursor_pos = QtGui.QCursor.pos()
        self._tip_label.move(cursor_pos.x() + 16, cursor_pos.y() + 16)

        if idx == self._tooltip_last_idx and self._tip_label.isVisible():
            return
        self._tooltip_last_idx = idx

        if in_conv or in_iter:
            text = f"time | {t:.1f} s\nConv failures: {self.conv_data[idx]:.1%}\nIter fraction: {self.iter_data[idx]:.1%}"
        else:
            lines = [f"time | {t:.1f} s"]
            display_items = self.y1_items + (self.y2_items if self.y2_items is not None else [])
            ncmax = -1
            for item in display_items:
                ncmax = max(ncmax, len(item.name))
            for j, item in enumerate(display_items):
                unit_str = f" {item.units}" if item.units else ""
                # which axis data are we pointing at?
                axis_data = self.y1_data if j < len(self.y1_items) else self.y2_data
                jj = j if j < len(self.y1_items) else j - len(self.y1_items)
                lines.append(f"{item.name:<{ncmax}s} | {axis_data[jj, idx]:.4g}{unit_str}")
            text = "\n".join(lines)

        self._tip_label.setFont(QtGui.QFont("Consolas", OnlinePlotter.font_size_pt-2))
        self._tip_label.setText(text)
        self._tip_label.adjustSize()
        self._tip_label.show()

    def _on_range_changed_manually(self, viewRange):
        self._auto_follow = False
        if OnlinePlotter.follow_button is not None:
            OnlinePlotter.follow_button.setChecked(False)

    @classmethod
    def set_auto_follow(cls, value):
        for plotter in cls.instances:
            plotter._auto_follow = value
        if cls.follow_button is not None:
            cls.follow_button.setChecked(value)

    @classmethod
    def set_plot_tooltip(cls, value):
        for plotter in cls.instances:
            plotter._plot_tooltip_enabled = value
            if not value:
                plotter._tip_label.hide()
                plotter._tooltip_vline.hide()
                plotter._tooltip_last_idx = -1
        if cls.tooltip_button is not None:
            cls.tooltip_button.setChecked(value)

    @classmethod
    def process_ui_events(cls, force=False):
        """Process Qt events with throttling so UI activity cannot starve simulation."""
        if cls.app is None:
            return

        now = time.perf_counter()
        if not force and (now - cls._last_ui_event_process) < cls.UI_EVENT_PROCESS_INTERVAL_SEC:
            return

        try:
            cls.app.processEvents(QtCore.QEventLoop.AllEvents, cls.UI_EVENT_PROCESS_BUDGET_MS)
        except TypeError:
            # Compatibility fallback for bindings without max-time overload.
            cls.app.processEvents()

        cls._last_ui_event_process = now
