pymacula
--------
Python wrapper module for ``macula`` starspot code: https://www.cfa.harvard.edu/~dkipping/macula.html.

.. code-block::

    % git clone https://github.com/timothydmorton/pymacula.git
    % cd pymacula
    % python setup.py install

If you come across errors in the installation regarding something like "undefined reference to main," then you may have to install with

.. code-block::

    % FFLAGS=-fPIC LDFLAGS=-shared python setup.py install
    
    
