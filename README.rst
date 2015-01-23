pymacula
========
Python wrapper module for ``macula`` starspot code: https://www.cfa.harvard.edu/~dkipping/macula.html.

Installation
------------

Install from PyPI::

    % pip install pymacula
    
or from the current repository::

    % git clone https://github.com/timothydmorton/pymacula.git
    % cd pymacula
    % python setup.py install

If you come across an error in the installation regarding something like ``undefined reference to main``, then you may have to install with

.. code-block::

    % FFLAGS=-fPIC LDFLAGS=-shared python setup.py install
    
For usage examples, see `example notebook <http://nbviewer.ipython.org/github/timothydmorton/pymacula/blob/master/notebooks/examples.ipynb>`_. 
