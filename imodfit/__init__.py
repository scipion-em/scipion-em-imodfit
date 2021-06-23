# **************************************************************************
# *
# * Authors:     Scipion Team
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from os.path import join
import pwem
import os

from pyworkflow.utils import Environ
from .constants import *

_references = ['Lopez-Blanco JR and Chacon P.2013']
__version__ = '1.51'


class Plugin(pwem.Plugin):
    _homeVar = IMODFIT_HOME
    _pathVars = [IMODFIT_HOME]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(IMODFIT_HOME, IMODFIT + '-' + IMODFIT_DEFAULT_VERSION)

    @classmethod
    def getImodfitEnviron(cls):
        """ Setup the environment variables needed to launch imodfit. """
        environ = Environ(os.environ)
        runtimePath = join(pwem.Config.EM_ROOT, IMODFIT + '-' + IMODFIT_DEFAULT_VERSION)

        # Add required disperse path to PATH and pyto path to PYTHONPATH
        environ.update({'LD_LIBRARY_PATH': os.pathsep.join(['~/intel/oneapi/mkl/latest/lib/intel64',
                                                            '~/intel/oneapi/mkl/latest/lib/ia32'])
                        })
        return environ

    @classmethod
    def defineBinaries(cls, env):
        # At this point of the installation execution cls.getHome() is None, so the em path should be provided
        pluginHome = join(pwem.Config.EM_ROOT, IMODFIT + '-' + IMODFIT_DEFAULT_VERSION)
        installationCmd = 'wget %s -O %s && ' % (cls._getImodfitDownloadUrl(), cls._getImodfitTxz())
        installationCmd += 'tar -xf %s --strip-components 1 && ' % cls._getImodfitTxz()
        installationCmd += 'rm %s && ' % cls._getImodfitTxz()

        #Installing required libraries
        installationCmd += 'wget %s -O %s && ' % (cls._getLibrariesDownloadUrl(), cls._getMKLsh())
        installationCmd += 'chmod +x %s && ' % cls._getMKLsh()
        installationCmd += 'sh %s -a -s --eula accept && ' % \
                           (cls._getMKLsh())

        #Creating validation file
        IMODFIT_INSTALLED = '%s_installed' % IMODFIT
        installationCmd += 'touch %s' % IMODFIT_INSTALLED  # Flag installation finished

        env.addPackage(IMODFIT,
                       version=IMODFIT_DEFAULT_VERSION,
                       tar='void.tgz',
                       commands=[(installationCmd, IMODFIT_INSTALLED)],
                       neededProgs=["wget", "tar"],
                       default=True)

    @classmethod
    def runIMODfit(cls, protocol, program, args, cwd=None):
        """ Run IMODFIT command from a given protocol. """
        protocol.runJob(join(cls.getHome('bin', program)), args, cwd=cwd, env=cls.getImodfitEnviron())

    @classmethod
    def getMCRPath(cls):
        return cls.getHome(IMODFIT)

    @staticmethod
    def _getImodfitTxz():
        return 'iMODFIT_v1.51_Linux_20190228.txz'

    @staticmethod
    def _getMKLsh():
      return 'intel_mkl_lib.sh'

    @classmethod
    def _getImodfitDownloadUrl(cls):
        return '\"https://chaconlab.org/hybrid4em/imodfit/imodfit-donwload?task=callelement&format=raw&item_id=23&element=f85c494b-2b32-4109-b8c1-083cca2b7db6&method=download&args[0]=d3da3168d547801010043ca84f4b3d2f\"'

    @classmethod
    def _getLibrariesDownloadUrl(cls):
        return 'https://registrationcenter-download.intel.com/akdlm/irc_nas/17757/l_onemkl_p_2021.2.0.296_offline.sh'