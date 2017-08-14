from setuptools import (setup,
                        find_packages)

from alcor.config import PROJECT_NAME

setup(name=PROJECT_NAME,
      packages=find_packages(),
      version='0.0.0',
      description='Library for astrophysical researches.',
      author='Georgy Skorobogatov, '
             'Azat Ibrakov',
      author_email='skorobogatov@phystech.edu, '
                   'azatibrakov@gmail.com',
      url='https://github.com/lycantropos/alcor',
      download_url='https://github.com/lycantropos/'
                   'alcor/archive/master.tar.gz',
      keywords=['astronomy'],
      install_requires=[
          'psycopg2>=2.7.1',
          'PyYAML>=3.12.0',  # settings loading
          'matplotlib>=2.0.2'  # plotting
      ],
      setup_requires=['pytest-runner>=2.11.1'],
      tests_require=[
          'pydevd>=1.0.0',  # debugging
          'pytest>=3.2.1',
          'pytest-cov>=2.4.0',
          'hypothesis>=3.13.0',
      ])
