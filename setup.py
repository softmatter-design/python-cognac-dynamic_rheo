"""Minimal setup file for chain_evaluate."""

from setuptools import setup, find_packages

setup(
    name='dynamic-rheo',
    version='0.0.1',
    license='proprietary',
    description='Module Experiment',

    author='hsasaki',
    author_email='hsasaki@softmatters.net',
    url='None.com',

    packages=find_packages(where='src'),
    package_dir={'': 'src'},

    entry_points={
        "console_scripts": [
          'dr_setup = evaluate_all:main'
        ]
    }
)