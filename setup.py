from setuptools import find_packages, setup

def parse_requirements(filename):
    with open(filename, 'r') as f:
        return f.read().splitlines()

setup(
    name='mymd',
    packages=find_packages(),
    install_requires=parse_requirements('requirements.txt'),
)