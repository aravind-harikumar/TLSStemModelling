# -*- coding: utf-8 -*-

"""
***************************************************************************
    OpenViewerAction.py
    ---------------------
    Date                 : August 2012
    Copyright            : (C) 2012 by Victor Olaya
    Email                : volayaf at gmail dot com
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

__author__ = 'Victor Olaya'
__date__ = 'August 2012'
__copyright__ = '(C) 2012, Victor Olaya'
# This will get replaced with a git SHA1 when you do a git archive
__revision__ = '$Format:%H$'

from sextante.gui.ToolboxAction import ToolboxAction
import os
from PyQt4 import QtGui
import subprocess
from sextante.lidar.fusion.FusionUtils import FusionUtils

class OpenViewerAction(ToolboxAction):

    def __init__(self):
        self.name="Open Fusion LAS viewer"
        self.group="Visualization"

    def getIcon(self):
        return QtGui.QIcon(os.path.dirname(__file__) + "/../images/tool.png")

    def execute(self):
        f = os.path.join(FusionUtils.FusionPath(), "pdq.exe")
        if os.path.exists(f):
            subprocess.Popen(f)
        else:
            QtGui.QMessageBox.critical(None, "Unable to open viewer", "The current Fusion folder does not contain the viewer executable.\nPlease check the configuration in the SEXTANTE settings dialog.")
