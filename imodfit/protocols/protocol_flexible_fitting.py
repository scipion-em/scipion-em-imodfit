# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     you (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from pyworkflow.protocol import Protocol, params
from pwem.objects.data import AtomStruct
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
        cgChoices = self._get_cgChoices()

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
        group = form.addGroup('Parameters')
        group.addParam('resolution', params.IntParam,
                      default=10,
                      label='Resolution in Angstroms',
                      help='Resolution in Angstroms. The resolution criterion follows EMAN package procedures')
        group.addParam('cutoff', params.FloatParam,
                      default=0,
                      label='EM density map threshold',
                      help='EM density map threshold. All density levels below this value will not be considered.')
        group.addParam('maxIter', params.IntParam,
                      default=10000,
                      label='Maximum iterations',
                      help='Maximum number of iterations')
        group.addParam('cgModel', params.EnumParam,
                      choices=cgChoices, default=2,
                      label='Coarse-Grained model',
                      help='Coarse-Grained model: 0=CA, 1=3BB2R, 2=Full-Atom, 3=NCAC(experimental) (default=2)')
        group.addParam('chiAngle', params.BooleanParam,
                      default=False, expertLevel=params.LEVEL_ADVANCED,
                      label='CHI dihedral angle',
                      help='Considers first CHI dihedral angle (default=disabled)')
        group.addParam('modesRange', params.FloatParam,
                      default=0.2, expertLevel=params.LEVEL_ADVANCED,
                      label='Used modes range',
                      help='Used modes range, either number [1,N] <integer>, or ratio [0,1) <float> (default=0.2).')
        group.addParam('excitedModesRange', params.FloatParam,
                      default=0.02, expertLevel=params.LEVEL_ADVANCED,
                      label='Excited modes range',
                      help='Excited modes range, either number [1,nevs] <integer>, or ratio [0,1) <float> (default=0.02)')

        group = form.addGroup('Output')
        group.addParam('outputBasename', params.StringParam,
                      default='imodfit', expertLevel=params.LEVEL_ADVANCED,
                      label='Output basename',
                      help='Output files basename (default=imodfit)')
        group.addParam('fullAtom', params.BooleanParam,
                      default=False, expertLevel=params.LEVEL_ADVANCED,
                      label='Full-atom output',
                      help='Enables full-atom output models')
        group.addParam('outputMovie', params.BooleanParam,
                      default=True,
                      label='Outputs a Multi-PDB',
                      help=' 	Outputs a Multi-PDB trajectory movie (<basename_movie>.pdb)')

        group = form.addGroup('Extra')
        group.addParam('extraParams', params.StringParam,
                       default='',
                       label='Extra parameters',
                       help='Extra parameters as expected from the command line.'
                            'https://chaconlab.org/hybrid4em/imodfit/imodfit-intro')

    def _getImportChoices(self):
      return ['File', 'Scipion Object']
    def _get_cgChoices(self):
      return ['CA', '3BB2R', 'Full-Atom', 'NCAC']

    def _getImodfitArgs(self):
      ccp4File = os.path.abspath(self.ccp4File)
      #Standard arguments
      args = [self.pdbFile, ccp4File, self.resolution.get(), self.cutoff.get()]
      args += ['-i {}'.format(self.maxIter.get()), '-m {}'.format(self.cgModel.get())]
      if self.chiAngle.get():
        args += ['-x']
      args += ['-n {}'.format(self.modesRange.get()), '-e {}'.format(self.excitedModesRange.get())]

      #Output arguments
      if self.outputBasename.get() != 'imodfit':
        args += ['-o {}'.format(self.outputBasename.get())]
      if self.fullAtom.get():
        args += ['-F']
      if self.outputMovie.get():
        args += ['-t']

      #Extra parameters
      if self.extraParams.get() != '':
        args += ['{}'.format(self.extraParams.get())]
      return args

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
      Plugin.runIMODfit(self, 'imodfit_mkl',
                        args=self._getImodfitArgs(),
                        cwd=self._getExtraPath())

    def createOutputStep(self):
        fittedPDB = AtomStruct(self._getExtraPath('{}_fitted.pdb'.format(self.outputBasename.get())))
        moviePDB = AtomStruct(self._getExtraPath('{}_movie.pdb'.format(self.outputBasename.get())))
        if self.importFrom.get() == self.FROM_SCIPION:
          fittedPDB.setVolume(self.inputVolume.get())
          moviePDB.setVolume(self.inputVolume.get())

        elif self.importFrom.get() == self.FROM_FILE:
          pass
        self._defineOutputs(fittedPDB=fittedPDB)
        self._defineOutputs(moviePDB=moviePDB)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods
