from setuptools import setup

setup(name='Lace',
      version='1.14.1',
      packages=['Lace'],
      license='GPL3',
      install_requires=['pandas',
                        'networkx',
                        'numpy',
                        'matplotlib'],
      entry_points={
          'console_scripts': [
              'BuildSuperTranscript=Lace.BuildSuperTranscript:main',
              'Lace_Checker=Lace.Checker:main',
              'Lace=Lace.Lace:main',
              'Mobius-as=Lace.Mobius_as:main',
              'Mobius=Lace.Mobius:main',
              'STViewer=Lace.STViewer:main'
          ]
      },
)
