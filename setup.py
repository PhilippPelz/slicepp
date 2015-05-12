from distutils.core import setup

setup(name='qstem',
      version='1.0',
      description='Python Distribution Utilities',
      author='philipp',
      author_email='gward@python.net',
      url='https://www.python.org/sigs/distutils-sig/',
      packages=['qstem', 'qstem.fileio','qstem.models','qstem.dialogs','qstem.models.raw_models'],
      package_dir={'qstem':'python'}
     ) 
