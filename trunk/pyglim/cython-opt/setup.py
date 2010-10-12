from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("pyglimCGrid", ["pyglimCGrid.pyx"]),
               Extension("pyglimCThck", ["pyglimCThck.pyx"])]

setup(
  name = 'PyGlim',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
