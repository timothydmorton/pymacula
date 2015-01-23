pymacula
========
Python wrapper module for ``macula`` starspot code: https://www.cfa.harvard.edu/~dkipping/macula.html.

.. code-block::

    % git clone https://github.com/timothydmorton/pymacula.git
    % cd pymacula
    % python setup.py install

If you come across an error in the installation regarding something like ``undefined reference to main``, then you may have to install with

.. code-block::

    % FFLAGS=-fPIC LDFLAGS=-shared python setup.py install
    
Usage
-----
You can generate a random starspot model::

    >>> from pymacula import MaculaModel
    >>> model = MaculaModel() #default is 3 random spots

You can call the model at random times, e.g., for plotting::

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> ts = np.arange(0,100,0.1)
    >>> plt.plot(ts, model(ts))
    
You can also make a model with non-default star parameters::

    >>> model = MaculaModel(star={'Peq':15}) #rotation period of 15 days
    
    
