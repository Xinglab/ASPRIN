import sys, os

# should be able to safely do this now.
from setuptools import setup, find_packages

setup(name='asprin',
      version='1.0.0',
      packages = find_packages('src'),  # include all packages under src
                        package_dir = {'':'src'},   # all distutils packages are under src
      entry_points={'console_scripts': ['asprin=asprin.scripts.asprin:main']},
                        description = 'asprin',
                        author = 'Emad Bahrami-Samani',
                        author_email = 'ebs@ucla.edu',
                        url = 'https://github.com/xinglab/asprin',
                        download_url = 'https://github.com/xinglab/asprin/tarball/asprin_1.0.0',
                        license='GPL3',
                        keywords = [],
                        classifiers = [],
)

