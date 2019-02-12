from setuptools import setup, find_packages
from codecs import open

with open("README.rst", "r") as f:
    long_description = f.read()

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()


setup(
    name='fruitbat',
    author='Adam Batten',
    author_email='adamjbatten@gmail.com',
    url='https://github.com/abatten/fruitbat',
    version='0.0.1rc1',
    description='Calculate the redshift of a FRB from its dispersion measure',
    package_dir={"fruitbat": "fruitbat"},
    long_description=long_description,
    install_requires=requirements,

    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ],
    packages=find_packages()
)
