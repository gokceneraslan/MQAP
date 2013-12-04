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

import os
import tempfile
import shutil

import numpy as np

from prody.atomic import AtomGroup
from prody.utilities import gunzip, which

from prody import parsePDB
from prody import fetchPDB

import fasta

__all__ = ['execPSIPRED', 'parsePSIPRED', 'performPSIPRED']


def execPSIPRED(pdb, outputname=None, outputdir=None):
    runpsipred = which('runpsipred')
    if runpsipred is None:
        raise EnvironmentError('command not found: runpsipred executable is not '
                               'found in the system path')

    assert outputname is None or isinstance(outputname, str),\
        'outputname must be a string'
    assert outputdir is None or isinstance(outputdir, str),\
        'outputdir must be a string'

    if not os.path.isfile(pdb):
        pdb = fetchPDB(pdb, compressed=False)
    if pdb is None:
        raise ValueError('pdb is not a valid PDB identifier or filename')
    if os.path.splitext(pdb)[1] == '.gz':
        if outputdir is None:
            pdb = gunzip(pdb, os.path.splitext(pdb)[0])
        else:
            pdb = gunzip(pdb, os.path.join(outputdir,
                         os.path.split(os.path.splitext(pdb)[0])[1]))

    tempdir = tempfile.mkdtemp()
    fastaname = fasta.convert(pdb, dir=tempdir)

    if outputdir is None:
        outputdir = '.'
    if outputname is None:
        out = os.path.join(outputdir,
                           os.path.splitext(os.path.split(pdb)[1])[0] +
                           '.psipred')
    else:
        out = os.path.join(outputdir, outputname + '.psipred')

    oldpath = os.getcwd()
    os.chdir(tempdir)

    status = os.system('{0} {1}'.format(runpsipred, fastaname))

    os.chdir(oldpath)
    ssfile = os.path.splitext(fastaname)[0] + ".ss2"
    shutil.copy(ssfile, out)
    shutil.rmtree(tempdir, ignore_errors=True)

    if status == 0:
        return out


def parsePSIPRED(psipred, ag):
    if not os.path.isfile(psipred):
        raise IOError('{0} is not a valid file path'.format(psipred))
    if not isinstance(ag, AtomGroup):
        raise TypeError('ag argument must be an AtomGroup instance')

    psipred = open(psipred)

    n_atoms = ag.numAtoms()
    SS = np.array(['-'] * n_atoms)
    COILSCORE = np.zeros(n_atoms, float)
    HELIXSCORE = np.zeros(n_atoms, float)
    STRANDSCORE = np.zeros(n_atoms, float)

    for line in psipred:
        if line.startswith('#'):
            continue

        psipred = line.split()
        if len(psipred) != 6:
            continue

        resindex = int(psipred[0])
        resname = psipred[1]
        ss = psipred[2]
        coilscore = float(psipred[3])
        helixscore = float(psipred[4])
        strandscore = float(psipred[5])

        # get the last chain
        # TODO: generate multiple fasta files for oligomeric PDB files
        chainid = list(ag.iterChains())[-1].getChid()

        if ag.ca.copy()[chainid].getSequence()[resindex - 1].upper() != resname.upper():
            continue

        resnum = list(ag.ca.copy()[chainid].iterResidues())[resindex-1].getResnum()
        indices = ag[chainid][resnum].getIndices()

        SS[indices] = ss
        COILSCORE[indices] = coilscore
        HELIXSCORE[indices] = helixscore
        STRANDSCORE[indices] = strandscore

    ag.setData('psipred_ss', SS)
    ag.setData('psipred_coilscore', COILSCORE)
    ag.setData('psipred_helixscore', HELIXSCORE)
    ag.setData('psipred_strandscore', STRANDSCORE)
    return ag


def performPSIPRED(pdb):
    pdb = fetchPDB(pdb, compressed=False)
    return parsePSIPRED(execPSIPRED(pdb), parsePDB(pdb))
