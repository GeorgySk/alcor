from setuptools import (setup,
                        find_packages)

import alcor
from alcor.config import PROJECT_NAME

project_base_url = 'https://github.com/lycantropos/alcor/'
setup(name=PROJECT_NAME,
      version='0.0.0',
      description=alcor.__doc__,
      long_description=open('README.rst').read(),
      license='MIT',
      author='Georgy Skorobogatov, '
             'Azat Ibrakov',
      author_email='skorobogatov@phystech.edu, '
                   'azatibrakov@gmail.com',
      url=project_base_url,
      download_url=project_base_url + 'archive/master.tar.gz',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Web Environment',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: Implementation :: CPython',
          'Operating System :: POSIX',
          'Operating System :: MacOS :: MacOS X',
          'Operating System :: Microsoft :: Windows',
          'Topic :: Scientific/Engineering :: Physics',
      ],
      keywords=['astrophysics'],
      packages=find_packages(exclude=('tests',)),
      install_requires=[
          'psycopg2>=2.7.1',  # PostgreSQL driver
          'sqlalchemy>=1.1.14',  # ORM
          'PyYAML>=3.12.0',  # settings loading
          'pandas>=0.20.3',  # data analysis
          'numpy>=1.11.3',  # multidimensional arrays computations
          'matplotlib>=2.0.2',  # plotting
          'h5py>=2.7.1',  # hdf5 files support
      ],
      setup_requires=['pytest-runner>=2.11.1'],
      tests_require=[
          'pydevd>=1.1.1',  # debugging
          'pytest>=3.2.1',
          'pytest-cov>=2.4.0',
          'hypothesis>=3.28.0',
          'hypothesis_sqlalchemy>=0.0.2',
      ])
