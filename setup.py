from setuptools import setup

setup(name='gtom',
      version='0.4',
      description='computes the general topological overlap matrix of an undirected network',
      url='https://github.com/benmaier/gtom',
      author='Benjamin Maier',
      author_email='bfmaier@physik.hu-berlin.de',
      license='MIT',
      packages=['gtom'],
      install_requires=[
          'numpy',
          'scipy',
      ],
      zip_safe=False)
