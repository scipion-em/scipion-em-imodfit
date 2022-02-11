# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
from ..protocols import ImodfitFlexFitting
from ..constants import *
import pyworkflow.protocol.params as params
from pwem.viewers import Chimera, ChimeraViewer, VmdView, EmProtocolViewer
from distutils.spawn import find_executable

VOLUME_CHIMERA, VOLUME_VMD = 0, 1
FITTED_PDB, MOVIE_PDB = 0, 1

class viewerImodfitVMD(EmProtocolViewer):
  _label = 'Viewer pdb movie'
  _targets = [ImodfitFlexFitting]

  def _defineParams(self, form):
    form.addSection(label='Visualization of output PDB and PDB movie')
    form.addParam('displayOutput', params.EnumParam,
                  choices=['Fitted PDB', 'Movie PDB'],
                  default=MOVIE_PDB,
                  display=params.EnumParam.DISPLAY_HLIST,
                  label='Output PDB to display',
                  help='*FittedPDB*: display final fitted structure\n'
                       '*VMD*: display fitting PDB movie'
                  )
    form.addParam('displayPDB', params.EnumParam,
                  choices=['Chimerax', 'VMD'],
                  default=VOLUME_VMD,
                  display=params.EnumParam.DISPLAY_HLIST,
                  label='Display output PDB with',
                  help='*Chimerax*: display AtomStruct as cartoons with '
                       'ChimeraX.\n *VMD*: display AtomStruct and movies.'
                  )

  def _getVisualizeDict(self):
    return {
      'displayPDB': self._showPDB,
    }

  def _validate(self):
    if (self.displayPDB == VOLUME_CHIMERA
            and find_executable(Chimera.getProgram()) is None):
      return ["chimera is not available. "
              "Either install it or choose option 'slices'. "]
    return []

  # =========================================================================
  # ShowPDBs
  # =========================================================================

  def _showPDB(self, paramName=None):
    if self.displayPDB == VOLUME_CHIMERA:
      return self._showPDBChimera()

    elif self.displayPDB == VOLUME_VMD:
      return self._showPDBVMD()

  def _showPDBChimera(self):
    """ Create a chimera script to visualize selected PDB. """
    if self.displayOutput.get() == FITTED_PDB:
      outputAS = self.protocol.fittedAtomStruct
    elif self.displayOutput.get() == MOVIE_PDB:
      outputAS = self.protocol.movieAtomStruct

    CViewer = ChimeraViewer(project=self.getProject(), protocol=self.protocol)
    return CViewer._visualize(outputAS)

  def _showPDBVMD(self):
    if self.displayOutput.get() == FITTED_PDB:
      tclName = 'fitted'
      outputAS = self.protocol.fittedAtomStruct
    elif self.displayOutput.get() == MOVIE_PDB:
      tclName = 'movie'
      outputAS = self.protocol.movieAtomStruct

    tclFile = self.protocol._getExtraPath('{}.tcl'.format(tclName))
    with open(tclFile, 'w') as f:
        f.write(TCL_STR % (os.path.abspath(outputAS.getFileName())))
                           #os.path.abspath(outputAS.getVolume().getFileName())))
    args = '-e {}'.format(os.path.abspath(tclFile))
    return [VmdView(args, cwd=self.protocol._getExtraPath())]

