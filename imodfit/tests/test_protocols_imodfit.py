# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * MRC Laboratory of Molecular Biology (MRC-LMB)
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

from pyworkflow.tests import BaseTest, setupTestProject
from pwem.protocols import ProtImportVolumes, ProtImportPdb
import pwem
import shutil, os
from ..constants import IMODFIT, IMODFIT_DEFAULT_VERSION
from ..protocols import imodfitFlexFitting
from pwem import Domain


FROM_FILE, FROM_SCIPION = 0, 1

class TestImodfit(BaseTest):
    @classmethod
    def setUpClass(cls):
        testFolder = pwem.Config.EM_ROOT + '/{}-{}/imodfit_test/'.format(IMODFIT, IMODFIT_DEFAULT_VERSION)
        cls.pdbFile = testFolder + '1oel.pdb'
        cls.volFile = testFolder + '1sx4A.ccp4'
        cls.mrcFile = testFolder + '1sx4A.mrc'
        if not os.path.exists(cls.mrcFile):
            shutil.copy(cls.volFile, cls.mrcFile)

        setupTestProject(cls)
        cls._runImportPDB()
        cls._runImportVolumes()

    @classmethod
    def _runImportVolumes(cls):
        """ Run an Import volumes protocol. """
        protImportVol = cls.newProtocol(
            ProtImportVolumes,
            importFrom=0,
            filesPath=TestImodfit.mrcFile,
            samplingRate=2)
        cls.launchProtocol(protImportVol)
        cls.protImportVol = protImportVol

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=TestImodfit.pdbFile)
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    def _runIMODFIT(self):
        protImodfit = self.newProtocol(
            imodfitFlexFitting,
            inputVolume=self.protImportVol.outputVolume,
            inputAtomStruct=self.protImportPDB.outputPdb)

        self.launchProtocol(protImodfit)
        pdbOut = getattr(protImodfit, 'fittedPDB', None)
        self.assertIsNotNone(pdbOut)

    def test_IMODFIT_fromScipion(self):
        self._runIMODFIT()




