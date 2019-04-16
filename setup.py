from setuptools import setup, find_packages
from codecs import open
import os

with open("README.rst", "r") as f:
    long_description = f.read()

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()


def get_version():
    here = os.path.abspath(os.path.dirname(__file__))
    version_file = os.path.join(here, 'fruitbat', '__version__.py')

    with open(version_file) as vf:
        lines = vf.read()
        version = lines.strip("__version__ = ").strip("").strip()
        return version


setup(
    name='fruitbat',
    version='1.0.0',
    author='Adam Batten',
    author_email='adamjbatten@gmail.com',
    url='https://github.com/abatten/fruitbat',
    download_url="https://pypi.org/project/fruitbat",
    project_urls={
        'Documentation': "https://fruitbat.readthedocs.io",
        'Source Code': "https://github.com/abatten/fruitbat"
        },
    description='Calculate the redshift of a FRB from its dispersion measure',
    long_description=long_description,
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: BSD License",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        ],
    package_dir={"fruitbat": "fruitbat"},
    packages=find_packages(),
    package_data={'fruitbat': ['*.npz', '*.npy', '*.csv']},
    include_package_data=True,
    keywords=("FRB redshift astronomy astrophysics fast radio burst"),
)
