# -*- coding: utf-8 -*-

__version__ = "0.0.1"

try:
    __PYMACULA_SETUP__
except NameError:
    __PYMACULA_SETUP__ = False

if not __PYMACULA_SETUP__:
    from .macula import Star, Spot, MaculaModel, macula
        
