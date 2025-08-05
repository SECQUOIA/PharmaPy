"""
Legacy setup.py for backward compatibility.
The preferred build method is now using pyproject.toml with PEP 517/518.
"""
from setuptools import setup, find_packages
import os

# Read the content of requirements.txt into a list
with open("requirements.txt", "r", encoding="utf-8") as f:
    requirements = f.read().splitlines()

# Read the README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='PharmaPy',
    version='0.0.1',
    packages=find_packages(),
    author='Daniel Casas-Orozco',
    author_email='dcasasor@purdue.edu',
    description='PharmaPy: Dynamic simulation of pharmaceutical manufacturing systems',
    long_description=long_description,
    long_description_content_type="text/markdown",
    license='MIT',
    url='https://github.com/SECQUOIA/PharmaPy',
    project_urls={
        'Documentation': 'https://pharmapy.readthedocs.io/',
        'Bug Tracker': 'https://github.com/SECQUOIA/PharmaPy/issues',
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    install_requires=requirements,
    include_package_data=True,
    package_data={
        "PharmaPy": ["data/**/*"],
    },
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "sphinx",
            "sphinx-rtd-theme",
            "black",
            "flake8",
            "isort",
        ],
        "test": [
            "pytest>=6.0",
            "pytest-cov",
        ],
        "docs": [
            "sphinx",
            "sphinx-rtd-theme",
        ],
    },
)