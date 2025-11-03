import os
import subprocess
import shutil
from setuptools import setup, find_packages
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.develop import develop as _develop
from setuptools.command.install import install as _install


def read_version():
    # Try to read __version__ from PYXAID/__init__.py if present
    here = os.path.abspath(os.path.dirname(__file__))
    init_py = os.path.join(here, 'PYXAID', '__init__.py')
    version = '0.0.0'
    if os.path.exists(init_py):
        try:
            with open(init_py, 'r') as fh:
                for line in fh:
                    if line.strip().startswith('__version__'):
                        # expected: __version__ = 'x.y.z'
                        parts = line.split('=', 1)
                        if len(parts) > 1:
                            version = parts[1].strip().strip("\"'\n ")
                            break
        except Exception:
            pass
    return version


def run_make():
    here = os.path.abspath(os.path.dirname(__file__))
    src_dir = os.path.join(here, 'PYXAID', 'src_cpp')
    makefile = os.path.join(src_dir, 'Makefile')
    if os.path.exists(makefile):
        # clean first
        subprocess.check_call(['make', '-C', src_dir, 'clean'])

        print('Running make in:', src_dir)
        # Try to detect a C++ compiler from the environment (CXX/CPP) or PATH
        cxx = None
        # Prefer explicit environment variables if set
        for var in ('CXX', 'CPP'):
            val = os.environ.get(var)
            if val:
                # if it's an absolute path or found in PATH, accept it
                if os.path.isabs(val) or shutil.which(val):
                    cxx = val
                    break
        # Fallback: look for common compiler wrapper names (conda-prefixed or system)
        if not cxx:
            for name in ('x86_64-conda-linux-gnu-c++', 'x86_64-conda_cos6-linux-gnu-c++', 'g++', 'c++', 'clang++'):
                path = shutil.which(name)
                if path:
                    cxx = path
                    break

        # Build environment for make by copying current env and setting compiler vars
        env = os.environ.copy()
        if cxx:
            env['CXX'] = cxx
            # Some Makefiles expect CPP to refer to the C++ compiler wrapper; set it too
            env['CPP'] = cxx
            print('Using C++ compiler:', cxx)
        else:
            raise RuntimeError('No suitable C++ compiler found in environment or PATH.')

        subprocess.check_call(['make', '-C', src_dir, 'all', '-j8'], env=env)

        # After successful build, try to find the compiled shared library (*.so)
        # and copy it into the PYXAID package directory so it is included at install
        so_file = os.path.join(here, 'PYXAID', 'pyxaid_core.so')
        if os.path.exists(so_file):
            print('Built shared library found at:', so_file)
        else:
            raise RuntimeError('Build completed but no shared library found at expected location: ' + so_file)

    else:
        raise RuntimeError('Makefile not found in source directory: ' + src_dir)


class build_py(_build_py):
    def run(self):
        run_make()
        super().run()
        
        #move built .so to package directory
        here = os.path.abspath(os.path.dirname(__file__))
        so_file = os.path.join(here, 'PYXAID', 'pyxaid_core.so')
        if os.path.exists(so_file):
            self.mkpath(os.path.join(self.build_lib, 'PYXAID'))
            shutil.copy2(so_file, os.path.join(self.build_lib, 'PYXAID'))
            print('Copied built library to build directory:', os.path.join(self.build_lib, 'PYXAID'))
        else:
            raise RuntimeError('Build failed: no built shared library found')


class develop(_develop):
    def run(self):
        run_make()
        super().run()


class install(_install):
    def run(self):
        super().run()

here = os.path.abspath(os.path.dirname(__file__))
readme = os.path.join(here, 'README.md')
long_description = ''
if os.path.exists(readme):
    with open(readme, 'r', encoding='utf-8') as fh:
        long_description = fh.read()

setup(
    name='pyxaid',
    version=read_version(),
    description='PYXAID: tools for electronic structure / nonadiabatic dynamics',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='(see repository)',
    packages=find_packages(exclude=('dev_modules', 'tutorials')),
    include_package_data=True,
    zip_safe=False,
    cmdclass={
        'build_py': build_py,
        'develop': develop,
        'install': install
    }
)
