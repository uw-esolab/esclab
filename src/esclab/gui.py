
import wx

class Icon(wx.StaticBitmap):
    def __init__(self, parent, label, pos):
        super().__init__(parent, bitmap=wx.Bitmap(50, 50))
        self.label = label
        self.SetPosition(pos)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
        self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
        self.Bind(wx.EVT_MOTION, self.OnMouseMove)
        self.dragging = False

    def OnLeftDown(self, event):
        if not event.ControlDown():
            self.CaptureMouse()
            self.dragging = True
            self.drag_start_pos = event.GetPosition()
        event.Skip()

    def OnLeftUp(self, event):
        if self.HasCapture():
            self.ReleaseMouse()
        self.dragging = False

    def OnMouseMove(self, event):
        if self.dragging:
            pos = event.GetPosition()
            dx = pos.x - self.drag_start_pos.x
            dy = pos.y - self.drag_start_pos.y
            self.SetPosition(self.GetPosition() + wx.Point(dx, dy))

class Canvas(wx.Panel):
    def __init__(self, parent):
        super().__init__(parent)
        self.SetBackgroundColour(wx.WHITE)
        self.icons = []
        self.lines = []
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_LEFT_DOWN, self.OnLeftDown)
        self.Bind(wx.EVT_LEFT_UP, self.OnLeftUp)
        self.Bind(wx.EVT_MOTION, self.OnMouseMove)
        self.dragging_line = False

    def AddIcon(self, label, pos):
        icon = Icon(self, label, pos)
        self.icons.append(icon)

    def OnPaint(self, event):
        dc = wx.PaintDC(self)
        for line in self.lines:
            dc.DrawLine(line[0], line[1])

    def OnLeftDown(self, event):
        if event.ControlDown():
            pos = event.GetPosition()
            for icon in self.icons:
                if icon.GetRect().Contains(pos):
                    self.dragging_line = True
                    self.line_start_pos = pos
                    self.line_start_icon = icon
                    break

    def OnLeftUp(self, event):
        if self.dragging_line:
            pos = event.GetPosition()
            for icon in self.icons:
                if icon.GetRect().Contains(pos) and icon != self.line_start_icon:
                    self.lines.append((self.line_start_icon.GetPosition() + wx.Point(25, 25), icon.GetPosition() + wx.Point(25, 25)))
                    self.Refresh()
                    break
            self.dragging_line = False

    def OnMouseMove(self, event):
        if self.dragging_line:
            pos = event.GetPosition()
            self.Refresh()
            dc = wx.ClientDC(self)
            dc.DrawLine(self.line_start_pos, pos)

class MainFrame(wx.Frame):
    def __init__(self):
        super().__init__(None, title="Icon Drag and Connect", size=(800, 600))
        self.canvas = Canvas(self)
        toolbar = self.CreateToolBar()
        add_icon_tool = toolbar.AddTool(wx.ID_ANY, "Add Icon", wx.Bitmap(50, 50))
        toolbar.Realize()
        self.Bind(wx.EVT_TOOL, self.OnAddIcon, add_icon_tool)
        self.Show()

    def OnAddIcon(self, event):
        label = wx.GetTextFromUser("Enter icon label:", "Add Icon")
        if label:
            self.canvas.AddIcon(label, wx.Point(100, 100))

if __name__ == "__main__":
    app = wx.App(False)
    frame = MainFrame()
    app.MainLoop()
