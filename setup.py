#!/usr/bin/env python3
# -*- coding:utf-8 -*-

from setuptools import setup

with open("./LongDescription", 'r') as README_FILE:
    long_description = README_FILE.read()

setup(
    name='PathPandem',
    version='1.0.2.3',
    description='Simulate Pandemic Pathogen Outbreak',
    license="GPLv3",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Pradyumna Paranjape',
    author_email='pradyparanjpe@rediffmail.com',
    url="",
    packages=['PathPandem'],
    install_requires=['numpy', 'gooey', 'matplotlib'],
    scripts=['bin/PathPandem',],
    package_data={'PathPandem': ['reverse_cfr_database.pkl']},
)
