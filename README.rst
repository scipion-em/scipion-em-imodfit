==============================
Scipion - iMODfit plugin
==============================

Scipion plugin that offers access to the **IMODFIT software** for flexible alignment of protein structures to EM maps
https://chaconlab.org/hybrid4em/imodfit


=====
Setup
=====

- **Install this plugin in devel mode:**

Using the command line:

.. code-block::

    scipion3 installp -p local/path/to/scipion-em-imodfit --devel

iMODfit precompiled binaries and necessary libraries are downloaded during the installation

**- VMD viewer**

**VMD** binaries for the visualization of the results will **NOT** be downloaded automatically with the plugin.

The independent download of VMD by the user is required
for using this visualization software. You can download it in:
https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD

Then, *VMD_HOME* has to be set in *scipion.conf* file. (e.g: VMD_HOME=/usr/local/bin)

- **Contact information:**

If you experiment any problem, please contact us here: scipion-users@lists.sourceforge.net or open an issue
--> https://github.com/scipion-em/scipion-em-imodfit/issues

We'll be pleased to help.

*Scipion Team*


