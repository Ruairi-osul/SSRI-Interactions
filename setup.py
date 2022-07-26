from setuptools import setup, find_packages

with open("requirements.txt") as f:
    reqs = f.read().splitlines()

setup(
    name="ssri_interactions",
    version="0.0.1",
    author="Ruairi O'Sullivan",
    author_email="ruairi.osullivan.work@gmail.com",
    packages=find_packages(),
    python_requires=">=3.6",
    install_requires=reqs,
)
