from distutils.core import setup, Extension

contactlooper = Extension('latticeproteins.contactlooper', sources = ['src/contactlooper.c'])

setup (name = 'latticeproteins', 
       fullname = 'Lattice Protein Simulation Package',
       version = '0.1',
       author = 'Jesse D. Bloom',
       author_email = 'jbloom@fhcrc.rog',
       description = 'Code for lattice protein simulations.',
       platforms = 'Tested on Mac OS X and Linux',
       packages = ['latticeproteins'],
       package_dir = {'latticeproteins' : 'src'},
       ext_modules = [contactlooper]
)
