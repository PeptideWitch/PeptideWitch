from distutils.core import setup

setup(
    name='PepWitch',
    version='1.2',
    packages=['build.PepWitch1.1.setuptools-27.2.0-py3.6.egg.setuptools',
              'build.PepWitch1.1.setuptools-27.2.0-py3.6.egg.setuptools.extern',
              'build.PepWitch1.1.setuptools-27.2.0-py3.6.egg.setuptools.command',
              'build.PepWitch1.1.setuptools-27.2.0-py3.6.egg.pkg_resources',
              'build.PepWitch1.1.setuptools-27.2.0-py3.6.egg.pkg_resources.extern',
              'build.PepWitch1.1.setuptools-27.2.0-py3.6.egg.pkg_resources._vendor',
              'build.PepWitch1.1.setuptools-27.2.0-py3.6.egg.pkg_resources._vendor.packaging'],
    url='',
    license='Free to use, edit and distribute on condition of a citation for any research that employs the use of this code or any of its derivatives (in whole or in part)',
    author='David Handler',
    author_email='david.handler@students.mq.edu.au',
    description='NSAF and high stringency protein ID calculations based on adjusted Q value analysis.'
)

install_requires=['numpy 1.13.3',
                  'matplotlib 2.1.0',
                  'scipy 1.0.0',
                  'statsmodels 0.8.0']
