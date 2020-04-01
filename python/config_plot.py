# coding: utf8
import matplotlib


def set_backend():
    from sys import platform as _platform
    if _platform == "darwin":  # MAC OS X
        # https://markhneedham.com/blog/2018/05/04/python-runtime-error-osx-matplotlib-not-installed-as-framework-mac/
        matplotlib.use('TkAgg')
