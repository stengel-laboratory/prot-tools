import versioneer
import setuptools


setuptools.setup(name='prot_tools',
      description='Miscellaneous selection of scripts extracting information from pdb files',
      author='Kai Kammer',
      author_email='kai-michael.kammer@uni-konstanz.de',
      url='https://git.uni-konstanz.de/kai-michael-kammer/prot_tools',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      packages=setuptools.find_packages(),
      install_requires=[
          'numpy',
          'pandas',
          'biopython',
      ],
      license='MIT',
      python_requires='>=3.6'
      )
