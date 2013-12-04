## MQAP
# Copyright (C) 2013 Gokcen Eraslan
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import prody

import tempfile
import textwrap
import os

def convert(pdbStructure, dir=None):
    if os.path.isfile(pdbStructure):
        pdbStructure = prody.parsePDB(pdbStructure)

    fastahandle, fastaname = tempfile.mkstemp(suffix='.fasta', dir=dir,
                                              text=True)

    for chain in pdbStructure.iterChains():
        sequence = chain.getSequence()
        chid = chain.getChid()

        os.write(fastahandle, '>chain_%s\n%s\n' % (chid,
                                      '\n'.join(textwrap.wrap(sequence, 60))))
    os.close(fastahandle)
    return fastaname

