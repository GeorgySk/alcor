============================================
``hypothesis`` strategies for ``SQLAlchemy``
============================================

.. image:: https://travis-ci.org/wolvespack/alcor.svg?branch=master
  :target:  https://travis-ci.org/wolvespack/alcor

.. image:: https://codecov.io/gh/wolvespack/alcor/branch/master/graph/badge.svg
  :target: https://codecov.io/gh/wolvespack/alcor

.. image:: https://badge.fury.io/py/alcor.svg
  :target: https://badge.fury.io/py/alcor


In what follows ``python3`` is an alias for ``python3.5``
or any later version (``python3.6`` and so on).

Installation
------------
Install the latest ``pip`` & ``setuptools`` packages versions

.. code-block:: bash

  python3 -m pip install --upgrade pip setuptools

Release
~~~~~~~
Download and install the latest stable version from ``PyPI`` repository

.. code-block:: bash

  python3 -m pip install --upgrade alcor

Developer
~~~~~~~~~
Download and install the latest version from ``GitHub`` repository

.. code-block:: bash

  git clone https://github.com/wolvespack/alcor.git
  cd alcor
  python3 setup.py install

Running tests
-------------
Plain

.. code-block:: bash

    python3 setup.py test

Inside ``Docker`` container

.. code-block:: bash

    docker-compose up

Inside ``Docker`` container with remote debugger

.. code-block:: bash

    ./set-dockerhost.sh docker-compose -f docker-compose.debug.yml up
