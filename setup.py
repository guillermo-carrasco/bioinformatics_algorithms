#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='bioinformatics algorithms',
      description="Bioinformatic tools and routines used to solve the problems \
              proposed on coursera course Bioinformatics Algorithms 1",
      author='Guillermo Carrasco',
      author_email='guille.ch.88@gmail.com',
      url='http://mussolblog.wordpress.com',
      packages=find_packages(exclude=['tests']),
)
