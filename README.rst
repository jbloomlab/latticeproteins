=================================
Lattice protein simulator
=================================

This is a lattice protein simulator written by `Jesse Bloom`_ and extended by `Zachary Sailer`_ in the `Harms lab`_.

This simulator has been used in the following publications; if you use this program please cite these publications:

    * `Protein stability promotes evolvability`_

    * `Stability and the evolvability of function in a model protein`_

This software is distributed under the `GNU Public License`_, meaning you are free to use this software for pretty much whatever you want provided that you retain the license.

Unfortunately, I wrote this software package back in my PhD days before I was aware of the appropriate procedures for package-level documentation. So although the source code is well documented, there isn't any higher level documentation.

The package is written in Python.

`Version 0.1`_ should work with recent versions of Python 2 (such as 2.7 and probably 2.6).

`Version 0.2`_ has been ported to Python 3 by `Zachary Sailer`_.

, and should work with recent versions 2.*. It also includes a few C extensions, and so compilation requires the ``gcc`` compiler. To install the package from source::

    python setup.py build
    sudo python setup.py install

This package calculates exact thermodynamic stabilities within the model by summing over the entire ensemble of conformations to compute the partition function. When run, it first creates a database that stores all of these conformations, which can take a substantial amount of time. It will run quickly for short proteins <20 residues, but get increasingly slow and memory-intensive after that.


.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`Protein stability promotes evolvability`: http://www.ncbi.nlm.nih.gov/pubmed/16581913
.. _`Stability and the evolvability of function in a model protein`: http://www.ncbi.nlm.nih.gov/pubmed/15111394
.. _`GNU Public License`: http://www.gnu.org/licenses/gpl.html
.. _`Zachary Sailer`: https://github.com/Zsailer
.. _`Harms lab`: http://harmslab.uoregon.edu/
.. _`Version 0.1`: https://github.com/jbloom/latticeproteins/tree/v0.1
.. _`Version 0.2`: https://github.com/jbloom/latticeproteins/tree/v0.2
