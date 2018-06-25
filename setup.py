from setuptools import (setup,
                        find_packages)

from alcor.config import PROJECT_NAME

project_base_url = 'https://github.com/lycantropos/alcor/'
setup(name=PROJECT_NAME,
      packages=find_packages(),
      version='0.0.0',
      description='Library for astrophysical researches.',
      author='Georgy Skorobogatov, '
             'Azat Ibrakov',
      author_email='skorobogatov@phystech.edu, '
                   'azatibrakov@gmail.com',
      url=project_base_url,
      download_url=project_base_url + 'archive/master.tar.gz',
      keywords=['astrophysics'],
      install_requires=[
          'PyYAML>=3.12.0',  # settings loading
      ])
