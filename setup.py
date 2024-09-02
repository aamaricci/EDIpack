import os
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

class MesonBuildExt(build_ext):
    def build_extension(self, ext):
        # Path to the directory where Meson will output build files
        build_dir = os.path.abspath(self.build_temp)

        # Path to the directory where Meson will install the output
        install_dir = os.path.abspath(self.build_lib)

        # Create the build directory if it does not exist
        os.makedirs(build_dir, exist_ok=True)

        # Run Meson
        meson_command = [
            'meson',  # Call meson executable
            build_dir,  # Where to put build artifacts
            '--prefix', install_dir,  # Where to install artifacts
        ]

        # Run Ninja after Meson
        ninja_command = [
            'ninja', '-C', build_dir  # Run ninja in the build directory
        ]

        # Call Meson command
        subprocess.check_call(meson_command)

        # Call Ninja build process
        subprocess.check_call(ninja_command)

        # Optionally install the output if Meson is configured for it
        subprocess.check_call(['ninja', '-C', build_dir, 'install'])
        
# You can define your Python extensions here if needed
ext_modules = [
    Extension('edipy', sources=[])
]

setup(
    name='edipy',
    version='0.1',
    ext_modules=ext_modules,  # Include extensions, even if they are dummy here
    cmdclass={
        'build_ext': MesonBuildExt,  # Use the custom Meson build command
    },
)