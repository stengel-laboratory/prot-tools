import versioneer
import setuptools


setuptools.setup(name='prot_tools',
      description='Miscellaneous selection of scripts extracting information from pdb files',
      author='Kai-Michael Kammer',
      author_email='kai-michael.kammer@uni-konstanz.de',
      url='https://github.com/stengel-laboratory/prot-tools',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      packages=setuptools.find_packages(),
      install_requires=[
          'numpy',
          'pandas',
          'biopython',
          'scipy',
          'matplotlib',
          'python-libsbml'
          'jinja2',
          'plotly',
      ],
      license='MIT',
      python_requires='>=3.10'
      )
