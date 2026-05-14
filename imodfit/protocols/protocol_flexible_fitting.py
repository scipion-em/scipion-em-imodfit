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
    Performs flexible fitting of a protein structure to an electron
    microscopy density map in order to adapt an atomic model to the
    conformational landscape observed experimentally. The protocol is
    intended for situations in which rigid-body fitting is insufficient
    because the structure must undergo internal deformations to match
    the target density. A prior rigid fitting step is recommended to
    ensure that the initial position of the structure is already
    approximately aligned with the map. Tutorial:
    https://chaconlab.org/hybrid4em/imodfit/imodfit-intro

    AI Generated:

    Flexible Fitting Protocol (imodfitFlexFitting) - User Manual
        Overview

        The Flexible Fitting protocol refines the conformation of an
        atomic structure so that it better matches an experimental
        cryo-EM density map. Unlike rigid fitting approaches that only
        translate and rotate the structure as a whole, flexible fitting
        allows internal structural rearrangements, making it especially
        useful for proteins and macromolecular assemblies that undergo
        biologically relevant conformational changes.

        In cryo-EM studies, it is common for experimentally determined
        maps to capture states that differ from the available atomic
        models. Flexible fitting bridges this gap by deforming the model
        in a physically meaningful way while preserving overall
        structural integrity. This approach is particularly valuable
        when studying molecular motions, domain rearrangements, ligand-
        induced conformational changes, or assemblies with flexible
        regions.

        Inputs and Biological Context

        The protocol requires two principal inputs: an electron density
        map and an atomic structure. The map represents the experimental
        target obtained from cryo-EM reconstruction, while the atomic
        structure provides the starting conformation that will be
        adapted to the map.

        A biologically meaningful fitting requires that the initial
        structure already occupies approximately the correct position
        and orientation inside the density. For this reason, a rigid
        fitting stage is strongly recommended before performing flexible
        refinement. If the initial placement is poor, the resulting
        deformations may become unrealistic or fail to converge toward a
        meaningful solution.

        The protocol supports structures originating from both PDB and
        mmCIF formats, making it compatible with modern structural
        biology repositories and modeling workflows.

        Resolution and Density Threshold

        The target map resolution plays an important role in determining
        the level of structural detail that can be reliably interpreted.
        Lower-resolution maps generally support large-scale domain
        motions and overall conformational adaptation, whereas higher-
        resolution maps can guide more detailed local refinements.

        The density threshold determines which regions of the map
        contribute to the fitting process. Proper threshold selection is
        biologically important because it helps distinguish meaningful
        molecular density from background noise. Excessively low
        thresholds may introduce noise-driven distortions, while overly
        restrictive thresholds may exclude flexible or weakly resolved
        regions that are biologically relevant.

        Coarse-Grained Models and Structural Representation

        The protocol offers several coarse-grained representations that
        control the balance between computational efficiency and
        structural detail. Simplified representations are useful for
        exploring large conformational motions efficiently, particularly
        in very large macromolecular assemblies. Full-atom
        representations provide more detailed refinements but generally
        require greater computational effort.

        From a biological perspective, coarse-grained approaches are
        often sufficient for capturing collective domain movements or
        major conformational transitions. Full-atom refinement becomes
        more important when detailed structural interpretation or
        downstream atomic analysis is required.

        Conformational Flexibility and Normal Modes

        The fitting strategy relies on collective motions that describe
        biologically plausible deformations of the structure. The number
        and range of modes included during refinement determine how much
        flexibility is allowed during the adaptation process.

        Restricting the fitting to lower-frequency collective motions
        usually produces smoother and more biologically realistic
        deformations. Allowing too many highly localized motions may
        increase the risk of overfitting the experimental density,
        especially in noisy or low-resolution maps.

        The protocol also allows the inclusion of side-chain flexibility
        through dihedral-angle refinement. This option may improve local
        agreement with the density in well-resolved regions, although it
        also increases conformational complexity.

        Iterative Refinement Strategy

        Flexible fitting is performed iteratively, progressively adapting
        the structure to improve agreement with the density map. The
        maximum number of iterations controls how extensively the
        conformational search is explored.

        In practice, convergence behavior depends strongly on the
        quality of the starting model, the resolution of the map, and
        the intrinsic flexibility of the biological system. Large and
        heterogeneous complexes may require more extensive refinement
        than compact and rigid proteins.

        Outputs and Structural Interpretation

        The protocol produces a fitted atomic structure representing the
        refined conformation adapted to the target map. Optionally, it
        can also generate trajectory-like outputs that capture the
        progressive conformational changes during the fitting process.

        These trajectory outputs are especially valuable for biological
        interpretation because they help visualize the transition between
        the starting conformation and the final fitted state. Such
        information can provide insight into collective domain motions,
        hinge rearrangements, or flexible transitions associated with
        biological activity.

        Full-atom outputs may also be generated for downstream atomic
        analysis, visualization, or structural validation workflows.

        Practical Recommendations

        In routine cryo-EM workflows, it is advisable to begin with a
        well-positioned rigid fit and conservative flexibility settings.
        Excessive flexibility may produce unrealistic deformations,
        particularly when the map resolution is limited or when the
        density contains substantial noise.

        For large conformational transitions, coarse-grained models
        usually provide more stable behavior and better convergence.
        Higher-detail refinements should generally be reserved for maps
        with sufficient resolution to support detailed structural
        interpretation.

        Users should always visually inspect the fitted model relative
        to the density map and verify that biologically important
        regions remain structurally plausible after refinement.

        Final Perspective

        Flexible fitting is one of the most powerful strategies for
        integrating atomic structures with cryo-EM density maps when
        conformational variability is present. By allowing controlled
        structural adaptation guided by experimental density, the
        protocol enables biologically meaningful interpretation of
        dynamic molecular states while preserving realistic structural
        behavior.
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
