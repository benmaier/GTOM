from setuptools import setup

setup(name='gtom',
      version='0.3',
      description='computes the general topological overlap matrix of an undirected network',
      url='https://bitbucket.org/bfmaier/gtom',
      author='Benjamin Maier',
      author_email='bfmaier@physik.hu-berlin.de',
      license='MIT',
      packages=['gtom'],
      install_requires=[
          'numpy',
          'scipy',
      ],
      zip_safe=False)
