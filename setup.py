# Random Walk Segregation Index Multiplex Networks
# See full license in LICENSE.txt

from setuptools import setup, find_packages

classifiers = ['Development Status :: 4- Beta',
               'License :: OSI Approved :: MIT License',
               'Operating System :: OS Independent',
               'Intended Audience :: Science/Research',
               'Natural Language :: English',
               'Topic :: Scientific/Engineering :: GIS',
               'Topic :: Scientific/Engineering :: Information Analysis',
               'Topic :: Scientific/Engineering :: Mathematics',
               'Programming Language :: Python :: 3.6'
               ]
with open('requirements.txt') as f:
    requirements_lines=f.readlines()
install_requires=[r.strip() for r in requirements_lines]

setup(
    name="multiseg",
    version="0.1.0",
    description="Create a multiplex network from various transport networks, and calculate Random Walk Segregation Index.",
    classifiers=classifiers,
    author='Mateo Neira',
    author_email='mateo.neira.16@ucl.ac.uk',
    license='MIT',
    platforms='any',
    packages='multiseg',
    packages=find_packages(exclude=['*test']),
    install_requires=install_requires
)

