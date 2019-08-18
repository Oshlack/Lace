from setuptools import setup

setup(name='Lace',
      version='1.14',
      packages=['Lace'],
      license='GPL3',
      install_requires=['pandas',
                        'networkx',
                        'numpy',
                        'matplotlib'],
      entry_points={
          'console_scripts': [
              'BuildSuperScript=Lace.BuildSuperScript:main',
              'Lace_Checker=Lace.Checker:main',
              'Lace=Lace.Lace:main',
              'Mobius-as=Lace.Mobius-as:main',
              'Mobius=Lace.Mobius:main',
              'STViewer=Lace.STViewer:main',
          ]
      },
)
