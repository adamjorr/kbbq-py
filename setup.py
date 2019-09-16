import setuptools
import kbbq

#see https://packaging.python.org/tutorials/packaging-projects/
#see https://packaging.python.org/guides/distributing-packages-using-setuptools/

with open("README.md", 'r') as fh:
    long_description = fh.read()

with open("requirements.txt", 'r') as fh:
    install_requires = fh.read().split() 

setuptools.setup(
    name="kbbq",
    version=kbbq.__version__,
    author="Adam Orr",
    author_email="ajorr1@asu.edu",
    description="A library and command line tool for reference-free base quality score recalibration.",
    long_description = long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/adamjorr/kbbq",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 2 - Pre-Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires='>=3',
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
            'kbbq=kbbq.main:main'
        ]
    }
)

