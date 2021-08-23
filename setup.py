import os
import glob
from setuptools import setup, find_packages

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

extra_files = package_files('easyvasp/file')+package_files('easyvasp/Pseudopotentials')

setup(
    name = "easyvasp",
    version = "1.0.1",
    keywords = 'high-throughput vasp pymatgen',
    author = "jialiang hou",
    author_email="675118082@qq.com",
    description="easyvasp is a python package to facilitate "
                "High throughput modeling and high throughput analysis for DFT calculation",
    install_requires=["numpy>=1.18.5", "pymatgen>=2020.8.13", "matplotlib>=3.2.2",'scipy>=1.5.0','scikit-learn>=0.23.1','pandas>=1.0.5'],
    license = "MIT License",
    classifiers=[
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.6",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Software Development :: Libraries :: Python Modules"
            ],
    packages = find_packages(),
    include_package_data = True,
    platforms = "any",
    package_data={'': extra_files},
    zip_safe = False,
)
