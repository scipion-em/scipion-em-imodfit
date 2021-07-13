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
This protocol is used to perform a flexible fitting of a protein structure over a electron map.
A rigid fitting for ensuring their prior best positions is needed before performing this flexible fitting.
"""
from pyworkflow.protocol import Protocol, params
from pwem.objects.data import AtomStruct
from pwem.emlib.image import ImageHandler
from pwem.convert.atom_struct import toPdb
from pyworkflow.utils import Message
import pyworkflow.utils as pwutils
import os, shutil
from imodfit import Plugin

from pwem.convert import Ccp4Header


class imodfitFlexFitting(Protocol):
    """
    Performs flexible fitting of a protein structure to a map.
    A rigid fitting for ensuring their prior best positions is needed before performing this flexible fitting.
    Tutorial: https://chaconlab.org/hybrid4em/imodfit/imodfit-intro
    """
    _label = 'Flexible fitting'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        cgChoices = self._get_cgChoices()

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputVolume', params.PointerParam,
                      pointerClass='Volume', allowsNull=False,
                      label="Input volume",
                      help='Target EM map')

        form.addParam('inputAtomStruct', params.PointerParam,
                       pointerClass='AtomStruct', allowsNull=False,
                       label="Input atom structure",
                       help='Select the atom structure to be fitted in the volume')

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
                      help='Outputs a Multi-PDB trajectory movie (<basename_movie>.pdb)')

        group = form.addGroup('Extra')
        group.addParam('extraParams', params.StringParam,
                       default='',
                       label='Extra parameters',
                       help='Extra parameters as expected from the command line.'
                            'https://chaconlab.org/hybrid4em/imodfit/imodfit-intro')

    def _get_cgChoices(self):
      return ['CA', '3BB2R', 'Full-Atom', 'NCAC']

    def _getImodfitArgs(self):
      ccp4AbsPath = os.path.abspath(self.ccp4File)
      #ccp4AbsPath = '/home/danieldh/i2pc/software/em/iMODfit-1.51/imodfit_test/1oel.ccp4'
      #Standard arguments
      args = [self.pdbFile, ccp4AbsPath, self.resolution.get(), self.cutoff.get()]
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

    def _getMrcInputVol(self):
      inpVol = self.inputVolume.get()
      # Input volume as a ccp4 file
      name, ext = os.path.splitext(inpVol.getFileName())
      if ext != '.mrc':
        mrcFile = self._getExtraPath(pwutils.replaceBaseExt(inpVol.getFileName(), 'mrc'))
        _ih = ImageHandler()
        _ih.convert(inpVol, mrcFile)
      else:
        mrcFile = inpVol.getFileName()
      return mrcFile

    def _getPdbInputStruct(self):
      inpStruct = self.inputAtomStruct.get()
      name, ext = os.path.splitext(inpStruct.getFileName())
      if ext == '.cif':
          cifFile = inpStruct.getFileName()
          pdbFile = self._getExtraPath(pwutils.replaceBaseExt(cifFile, 'pdb'))
          toPdb(cifFile, pdbFile)

      else:
        pdbFile = inpStruct.getFileName()
      return os.path.abspath(pdbFile)
    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('imodfitStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
      inpVol = self.inputVolume.get()
      self.ccp4File = self._getExtraPath(pwutils.replaceBaseExt(inpVol.getFileName(), 'ccp4'))
      shutil.copy(inpVol.getFileName(), self.ccp4File)

      #Ensuring a proper ccp4 file header
      sampling = inpVol.getSamplingRate()
      origin = inpVol.getOrigin(force=True).getShifts()
      ccp4H = Ccp4Header(self.ccp4File)
      ccp4H.copyCCP4Header(origin, sampling, Ccp4Header.ORIGIN)

      self.pdbFile = self._getPdbInputStruct()

    def imodfitStep(self):
      Plugin.runIMODfit(self, 'imodfit_mkl',
                        args=self._getImodfitArgs(),
                        cwd=self._getExtraPath())

    def createOutputStep(self):
        fittedPDB = AtomStruct(self._getExtraPath('{}_fitted.pdb'.format(self.outputBasename.get())))
        moviePDB = AtomStruct(self._getExtraPath('{}_movie.pdb'.format(self.outputBasename.get())))
        fittedPDB.setVolume(self.inputVolume.get())
        moviePDB.setVolume(self.inputVolume.get())

        self._defineOutputs(fittedAtomStruct=fittedPDB)
        self._defineOutputs(movieAtomStruct=moviePDB)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods
