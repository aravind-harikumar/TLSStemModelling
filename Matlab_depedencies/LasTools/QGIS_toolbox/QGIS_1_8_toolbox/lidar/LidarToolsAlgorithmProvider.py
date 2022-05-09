# -*- coding: utf-8 -*-

"""
***************************************************************************
    LidarToolsAlgorithmProvider.py
    ---------------------
    Date                 : August 2012
    Copyright            : (C) 2012 by Victor Olaya
    Email                : volayaf at gmail dot com
    ---------------------
    Date                 : September 2013
    Copyright            : (C) 2013 by Martin Isenburg
    Email                : martin near rapidlasso point com
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

import os
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from sextante.core.AlgorithmProvider import AlgorithmProvider
from sextante.core.SextanteUtils import SextanteUtils
from sextante.core.SextanteConfig import Setting, SextanteConfig
from sextante.lidar.lastools.LAStoolsUtils import LAStoolsUtils
from sextante.lidar.lastools.lasground import lasground
from sextante.lidar.lastools.lasheight import lasheight
from sextante.lidar.lastools.lasclassify import lasclassify
from sextante.lidar.lastools.laszip import laszip
from sextante.lidar.lastools.lasindex import lasindex
from sextante.lidar.lastools.lasclip import lasclip
from sextante.lidar.lastools.lasthin import lasthin
from sextante.lidar.lastools.lasnoise import lasnoise
from sextante.lidar.lastools.lassort import lassort
from sextante.lidar.lastools.lastile import lastile
from sextante.lidar.lastools.lasgrid import lasgrid
from sextante.lidar.lastools.lasview import lasview
from sextante.lidar.lastools.lasboundary import lasboundary
from sextante.lidar.lastools.lasinfo import lasinfo
from sextante.lidar.lastools.las2dem import las2dem
from sextante.lidar.lastools.blast2dem import blast2dem
from sextante.lidar.lastools.las2iso import las2iso
from sextante.lidar.lastools.las2las_filter import las2las_filter
from sextante.lidar.lastools.las2las_project import las2las_project
from sextante.lidar.lastools.las2las_transform import las2las_transform
from sextante.lidar.lastools.blast2iso import blast2iso
from sextante.lidar.lastools.lasprecision import lasprecision
from sextante.lidar.lastools.lasvalidate import lasvalidate
from sextante.lidar.lastools.lasduplicate import lasduplicate
from sextante.lidar.lastools.las2txt import las2txt
from sextante.lidar.lastools.txt2las import txt2las
from sextante.lidar.lastools.las2shp import las2shp
from sextante.lidar.lastools.shp2las import shp2las
from sextante.lidar.lastools.lasmerge import lasmerge
from sextante.lidar.lastools.lassplit import lassplit
from sextante.lidar.lastools.lascanopy import lascanopy
from sextante.lidar.lastools.lascontrol import lascontrol
from sextante.lidar.lastools.lasoverage import lasoverage
from sextante.lidar.lastools.lasoverlap import lasoverlap
"""
"""
from sextante.lidar.fusion.OpenViewerAction import OpenViewerAction
from sextante.lidar.fusion.CanopyMaxima import CanopyMaxima
from sextante.lidar.fusion.CanopyModel import CanopyModel
from sextante.lidar.fusion.ClipData import ClipData
from sextante.lidar.fusion.CloudMetrics import CloudMetrics
from sextante.lidar.fusion.Cover import Cover
from sextante.lidar.fusion.GridMetrics import GridMetrics
from sextante.lidar.fusion.GridSurfaceCreate import GridSurfaceCreate
from sextante.lidar.fusion.GroundFilter import GroundFilter
from sextante.lidar.fusion.MergeData import MergeData
from sextante.lidar.fusion.FilterData import FilterData
from sextante.lidar.fusion.FusionUtils import FusionUtils

class LidarToolsAlgorithmProvider(AlgorithmProvider):

    def __init__(self):
        AlgorithmProvider.__init__(self)
        self.activate = False
        self.algsList = []
        if SextanteUtils.isWindows():
            lastools = [lasground(), lasheight(), lasclassify(), lasclip(), lastile(), lasgrid(), las2dem(),  blast2dem(), las2iso(), blast2iso(), lasview(), lasboundary(), lasinfo(), lasprecision(), lasvalidate(), lasduplicate(), las2txt(), txt2las(), laszip(), lasindex(), lasthin(), lassort(), lascanopy(), lasmerge(), las2shp(), shp2las(), lasnoise(), lassplit(),  las2las_filter(), las2las_project(), las2las_transform(), lasoverage(), lasoverlap()]
        else:
            lastools = [lasinfo(), lasprecision(), lasvalidate(), las2txt(), txt2las(), laszip(), lasindex(), lasmerge(), las2las_filter(), las2las_project(), las2las_transform()]
        for alg in lastools:
            alg.group = "LAStools"
        self.algsList.extend(lastools)

        if SextanteUtils.isWindows():
            self.actions.append(OpenViewerAction())
            fusiontools = [CloudMetrics(), CanopyMaxima(), CanopyModel(), ClipData(), Cover(), FilterData(),
                         GridMetrics(), GroundFilter(), GridSurfaceCreate(), MergeData()]
            for alg in fusiontools:
                alg.group = "Fusion"
            self.algsList.extend(fusiontools)

    def initializeSettings(self):
        AlgorithmProvider.initializeSettings(self)
        SextanteConfig.addSetting(Setting(self.getDescription(), LAStoolsUtils.LASTOOLS_FOLDER, "LAStools folder", LAStoolsUtils.LAStoolsPath()))
        SextanteConfig.addSetting(Setting(self.getDescription(), FusionUtils.FUSION_FOLDER, "Fusion folder", FusionUtils.FusionPath()))

    def getName(self):
        return "lidartools"

    def getDescription(self):
        return "Tools for LiDAR data"

    def getIcon(self):
        return QIcon(os.path.dirname(__file__) + "/../images/tool.png")

    def _loadAlgorithms(self):
        self.algs = self.algsList

    def getSupportedOutputTableExtensions(self):
        return ["csv"]
