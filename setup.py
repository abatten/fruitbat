from setuptools import setup, find_packages
from codecs import open

with open("README.rst", "r") as f:
    long_description = f.read()

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()


setup(
    name='fruitbat',
    version='0.0.2',
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
    keywords=("FRB redshift astronomy astrophysics fast radio burst"),
)
