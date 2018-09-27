from setuptools import setup
from distutils.util import convert_path

versionholder = {}
ver_path = convert_path('originator/__init__.py')
with open(ver_path) as init_file:
    exec(init_file.read(), versionholder)

setup(name='originator',
    version=versionholder['__version__'],
    description='Rotate circular replicon to start with dnaA gene',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.'
    ],
    keywords=['dnaA','ori','origin','replication','circular'],
    author='Ola Brynildsrud',
    author_email='ola.brynildsrud@fhi.no',
    packages=['originator'],
    install_requires=[
        'biopython'],
    entry_points={
        'console_scripts': ['originator=originator.originator:main']
    },
    package_data={'originator': ['setup.py','originator/db/*']},
    include_package_data=True
)