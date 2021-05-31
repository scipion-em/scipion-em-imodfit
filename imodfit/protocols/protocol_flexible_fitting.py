# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     you (you@yourinstitution.email)
# *
# * your institution
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
Describe your python module here:
This module will provide the traditional Hello world example
"""
from pyworkflow.protocol import Protocol, params, Integer
from pwem.objects.data import AtomStruct, Volume
from pyworkflow.utils import Message
import pyworkflow.utils as pwutils
import os
from imodfit import Plugin
import shutil

class imodfitFlexFitting(Protocol):
    """
    Performs flexible fitting of a protein structure to a map.
    Tutorial: https://chaconlab.org/hybrid4em/imodfit/imodfit-intro
    """
    _label = 'Flexible fitting'
    FROM_FILE = 0
    FROM_SCIPION = 1

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        importChoices = self._getImportChoices()

        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Volume')
        group.addParam('importFrom', params.EnumParam,
                      choices=importChoices, default=self.FROM_SCIPION,
                      label='Volume from',
                      help='Select the volume source')
        group.addParam('inputVolumeFile', params.FileParam,
                      allowsNull=False,
                      condition='importFrom==FROM_FILE',
                      label="Input volume file",
                      help='Target EM map file (mrc or ccp4)')
        group.addParam('inputVolume', params.PointerParam,
                      pointerClass='Volume', allowsNull=False,
                      condition='importFrom==FROM_SCIPION',
                      label="Input volume",
                      help='Target EM map')

        group = form.addGroup('Atom structure')
        group.addParam('importFromAtom', params.EnumParam,
                       choices=importChoices, default=self.FROM_SCIPION,
                       label='Atom structure from',
                       help='Select the atom structure source')
        group.addParam('inputAtomStruct', params.PointerParam,
                       pointerClass='AtomStruct', allowsNull=False,
                       condition='importFromAtom==FROM_SCIPION',
                       label="Input atom structure",
                       help='Select the atom structure to be fitted in the volume')
        group.addParam('inputAtomStructFile', params.FileParam,
                       allowsNull=False,
                       condition='importFromAtom==FROM_FILE',
                       label="Input pdb file",
                       help='Select the atom structure file to be fitted in the volume')

        form.addSection(label='Parameters')
        form.addParam('resolution', params.IntParam,
                      default=10,
                      label='Resolution in Angstroms',
                      help='Resolution in Angstroms. The resolution criterion follows EMAN package procedures')
        form.addParam('cutoff', params.FloatParam,
                      default=0,
                      label='EM density map threshold',
                      help='EM density map threshold. All density levels below this value will not be considered.')
        form.addParam('extraParams', params.StringParam,
                       default='-t ',
                       label='Extra parameters',
                       help='Extra parameters as expected from the command line.'
                            'https://chaconlab.org/hybrid4em/imodfit/imodfit-intro')

    def _getImportChoices(self):
      return ['File', 'Scipion Object']

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('imodfitStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
      if self.importFrom.get() == self.FROM_SCIPION:
        mrcFile = self.inputVolume.get().getFileName()
      else:
        mrcFile = self.inputVolumeFile.get()
      self.ccp4File = pwutils.replaceBaseExt(mrcFile, 'ccp4')
      shutil.copy(mrcFile, self.ccp4File)

      if self.importFromAtom.get() == self.FROM_SCIPION:
        self.pdbFile = os.path.abspath(self.inputAtomStruct.get().getFileName())
      else:
        self.pdbFile = os.path.abspath(self.inputAtomStructFile.get())

    def imodfitStep(self):
      ccp4File = os.path.abspath(self.ccp4File)
      Plugin.runIMODfit(self, 'imodfit_mkl',
                        args=[self.pdbFile, ccp4File, self.resolution.get(), self.cutoff.get(), self.extraParams.get()],
                        cwd=self._getExtraPath())

    def createOutputStep(self):
        fittedPDB = AtomStruct(self._getExtraPath('imodfit_fitted.pdb'))
        if self.importFrom.get() == self.FROM_SCIPION:
          fittedPDB.setVolume(self.inputVolume.get())
        elif self.importFrom.get() == self.FROM_FILE:
          pass
        self._defineOutputs(fittedPDB=fittedPDB)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods
