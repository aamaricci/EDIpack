from setuptools import find_packages
import pkgconfig
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
import sys

ED_libraries=pkgconfig.parse('scifor')['libraries']
ED_library_dirs=pkgconfig.parse('scifor')['library_dirs']
ED_include_dirs=pkgconfig.parse('scifor')['include_dirs']

ED_libraries.extend(pkgconfig.parse('mpi')['libraries'])
ED_library_dirs.extend(pkgconfig.parse('mpi')['library_dirs'])
ED_include_dirs.extend(pkgconfig.parse('mpi')['include_dirs'])

ED_libraries.extend(pkgconfig.parse('lanc_ed')['libraries'])
ED_library_dirs.extend(pkgconfig.parse('lanc_ed')['library_dirs'])
ED_include_dirs.extend(pkgconfig.parse('lanc_ed')['include_dirs'])

ext1 = Extension(
    name='ed2py',
    sources=['src/ED_INPUT_VARS.f90','src/ed2py/ed2py.f90'],
    f2py_options=["--quiet"],
    libraries=ED_libraries,
    library_dirs=ED_library_dirs,
    include_dirs=ED_include_dirs,
    extra_f90_compile_args=["-O2", "-ffree-line-length-none","-cpp","-D_MPI"])


setup(
    name = "lancpy",
    version = "0.0.6",
    description = "LANC ED python API",
    author = "Adriano Amaricci",
    author_email = "amaricci@sissa.it",
    url='https://github.com/QcmPlab/LANC_ED',
    package_dir={"": "python"},
    packages=find_packages(where="python"),
    ext_modules=[ext1])


