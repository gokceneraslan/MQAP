# -*- coding: utf-8 -*-

import os

import numpy as np

from prody.atomic import AtomGroup
from prody.utilities import gunzip, which

from prody import parsePDB
from prody import fetchPDB

import fasta

__all__ = ['execNetSurfP', 'parseNetSurfP', 'performNetSurfP']


def execNetSurfP(pdb, outputname=None, outputdir=None):
    netsurfp = which('netsurfp')
    if netsurfp is None:
        raise EnvironmentError('command not found: netsurfp executable is not '
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

    fastaname = fasta.convert(pdb)

    if outputdir is None:
        outputdir = '.'
    if outputname is None:
        out = os.path.join(outputdir,
                           os.path.splitext(os.path.split(pdb)[1])[0] +
                           '.netsurfp')
    else:
        out = os.path.join(outputdir, outputname + '.netsurfp')

    status = os.system('{0} -a -i {1} > {2}'.format(netsurfp, fastaname, out))
    os.remove(fastaname)
    if status == 0:
        return out


def parseNetSurfP(netsurfp, ag):
    if not os.path.isfile(netsurfp):
        raise IOError('{0} is not a valid file path'.format(netsurfp))
    if not isinstance(ag, AtomGroup):
        raise TypeError('ag argument must be an AtomGroup instance')

    netsurfp = open(netsurfp)

    n_atoms = ag.numAtoms()
    EXPOSURE = np.array(['-'] * n_atoms)
    RSA = np.zeros(n_atoms, float)
    ASA = np.zeros(n_atoms, float)
    ALPHASCORE = np.zeros(n_atoms, float)
    BETASCORE = np.zeros(n_atoms, float)
    COILSCORE = np.zeros(n_atoms, float)

    for line in netsurfp:
        if line.startswith('#'):
            continue

        netsurfp = line.split()
        if len(netsurfp) != 10:
            continue

        exposure = netsurfp[0]
        resname = netsurfp[1]
        seqname = netsurfp[2]
        resindex = int(netsurfp[3])
        rsa = float(netsurfp[4])
        asa = float(netsurfp[5])
        alphascore = float(netsurfp[7])
        betascore = float(netsurfp[8])
        coilscore = float(netsurfp[9])

        #sequence names in fasta must be in 'chain_X' format, otherwise skip
        if not 'chain_' in seqname:
            continue

        chainid = seqname.split('_')
        chainid = chainid[1] if len(chainid) == 2 else ' '

        if ag.ca.copy()[chainid].getSequence()[resindex - 1].upper() != resname.upper():
            continue

        resnum = list(ag.ca.copy()[chainid].iterResidues())[resindex-1].getResnum()
        indices = ag[chainid][resnum].getIndices()

        EXPOSURE[indices] = exposure
        RSA[indices] = rsa
        ASA[indices] = asa
        ALPHASCORE[indices] = alphascore
        BETASCORE[indices] = betascore
        COILSCORE[indices] = coilscore

    ag.setData('netsurfp_exposure', EXPOSURE)
    ag.setData('netsurfp_rsa', RSA)
    ag.setData('netsurfp_asa', ASA)
    ag.setData('netsurfp_alphascore', ALPHASCORE)
    ag.setData('netsurfp_betascore', BETASCORE)
    ag.setData('netsurfp_coilscore', COILSCORE)
    return ag


def performNetSurfP(pdb):
    pdb = fetchPDB(pdb, compressed=False)
    return parseNetSurfP(execNetSurfP(pdb), parsePDB(pdb))
