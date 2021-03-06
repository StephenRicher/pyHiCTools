#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" 'Set of tools to performing HiC analysis.' """

# Note: To use the 'upload' functionality of this file, you must:
#   $ pipenv install twine --dev

from setuptools import setup, find_packages, Command
from shutil import rmtree
import sys
import os


def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as f:
        return f.read()


def get_info():
    info = {}
    with open('version.py') as fp:
        exec(fp.read(), info)
    return info


class UploadCommand(Command):
    """Support setup.py upload for twine."""

    description = 'Build and publish the package.'
    user_options = []
    version = get_info()['__version__']

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print(f'\033[1m{s}\033[0m')

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self, version=version):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(os.path.dirname(__file__), 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        os.system(f'{sys.executable} setup.py sdist bdist_wheel --universal')

        self.status('Uploading the package to PyPI via Twine…')
        os.system('twine upload dist/*')

        self.status('Pushing git tags…')
        os.system(f'git tag -a v{version}')
        os.system('git push --tags')

        sys.exit()


setup(
    name='pyHiCTools',
    author='Stephen Richer',
    author_email='sr467@bath.ac.uk',
    url='https://github.com/StephenRicher/pyHiCTools',
    scripts=['bin/pyHiCTools'],
    python_requires='>=3.6.0',
    install_requires=['pyCommonTools>=2.0'],
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.6',
        'Natural Language :: English',
    ],
    version=get_info()['__version__'],
    description=__doc__,
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    zip_safe=False,
)
