# -*- coding: utf-8 -*-

"""
***************************************************************************
    LAStoolsUtils.py
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

from PyQt4.QtCore import *
from PyQt4.QtGui import *
import subprocess
from sextante.core.SextanteLog import SextanteLog
from sextante.core.SextanteConfig import SextanteConfig

class LAStoolsUtils():

    LASTOOLS_FOLDER = "LASTOOLS_FOLDER"

    @staticmethod
    def LAStoolsPath():
        folder = SextanteConfig.getSetting(LAStoolsUtils.LASTOOLS_FOLDER)
        if folder == None:
            folder =""

        return folder

    @staticmethod
    def runLAStools(commands, progress):
        loglines = []
        loglines.append("LAStools console output")
        commandline = " ".join(commands)
        proc = subprocess.Popen(commandline, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE,stderr=subprocess.STDOUT, universal_newlines=False).stdout
        for line in iter(proc.readline, ""):
            loglines.append(line)
        SextanteLog.addToLog(SextanteLog.LOG_INFO, loglines)
