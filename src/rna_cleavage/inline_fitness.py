#! /usr/bin/env python

'''inline_fitness.py

Calculate inline fitness (i.e., cleavage propensity) for each residue in
each RNA chain in a PDB file.  '''

import sys
import pdb
import math

import numpy
from prody import parsePDB

__version__ = '0.1'

def inline_fitness(pdb_id, verbose):

    structure = parsePDB(pdb_id)

    for chain in structure.iterChains():

        if not is_rna(chain): continue

        for residue in chain.iterResidues():

            try:
                fitness = calc_inline_fitness(residue, verbose)
            except AttributeError:
                # end of the chain
                continue

            # pos = the position *downstream* of the examined
            # internucleotide.

            chain_id = chain.getChid()
            res_num = residue.getResnum()
            res_id = residue.getResname()

            fields = (chain_id, res_id, res_num, fitness)

            print '\t'.join(map(str, fields))

def is_rna(chain):
    ''' determines whether chain is an RNA or not.

    Args:
        chain (chain): XXX

    Returns:
        bool: True if RNA, False otherwise
    '''
    nucs = set(['A','G','C','U'])
    chainres = set(chain.getResnames())

    if nucs.intersection(chainres): return True

    return False

def calc_inline_fitness(residue, verbose):
    ''' Calculate in-line fitness (F) as described in Breaker & Soukup
    (1999) RNA. 

    Quote from the paper:

    The tau angle can be derived from a triangle whose sides are defined
    by the O2'-P and O2'-O5' interatomic distances and the P-O' bond
    length. A tau angle of 180 degrees is optimal for in-line attack
    Because of steric constraints, a tau angle of 45 degrees approximates
    the most unfavorable orientation that can be achieved by RNA

    Definitons:

        - dist_O2P = distance from O2' to P
        - dist_PO5 = distance from P to O5'
        - dist_O2O5 = distance from O2' to O5'
        
    Eqauation 2 from the paper:

        cos(tau) = (dist_O2P**2 + dist_PO5**2 - dist_O2O5**2) / 
                   (dist_O2P * dist_PO5 * 2)
   
    Equation 3 from the paper:

        F = ((tau - 45) / (180 - 45)) * (3**3 / dist_O2P**3)

    The first term of this equation quantifies the contribution of tau
    angle, and the second term quantifies contribution of attack distance.

    Args:
        residue (residue): pymol chain residue
        verbose (bool): verbose reporting

    Returns:
        F (float): in-line fitness
    '''

    pass

    dist_O2P = calc_distance(residue["O2'"], residue['P'])
    dist_PO5 = calc_distance(residue["P"], residue["O5'"])
    dist_O2O5 = calc_distance(residue["O2'"], residue["O5'"])

    cos_tau = (dist_O2P**2 + dist_PO5**2 - dist_O2O5**2) / \
              (2 * dist_O2P * dist_PO5)

    tau = math.degrees(math.acos(cos_tau))
   
    F = ((tau - 45.0) / (180.0 - 45.0)) * (3.0**3 / dist_O2P**3.0)

    return F

def calc_distance(res1, res2):
    ''' calculate distance between two xyz coordinates
    
    Solution from http://stackoverflow.com/a/1401828/4726866
    
    '''

    dist = numpy.linalg.norm(res1.getCoords() - res2.getCoords())

    return dist

def parse_options(args):
    from optparse import OptionParser

    usage = "%prog [OPTION]... PDB_ID"
    version = "%%prog %s" % __version__
    description = ("calculate inline fitness values for RNA from PDB files")

    parser = OptionParser(usage=usage, version=version,
                          description=description)

    parser.add_option("-v", "--verbose", action="store_true",
        default=False, help="verbose output (default: %default)")

    options, args = parser.parse_args(args)

    if len(args) != 1:
        parser.error("specify PDB_ID")

    return options, args

def main(args=sys.argv[1:]):

    options, args = parse_options(args)

    kwargs = {'verbose':options.verbose}

    pdb_id = args[0]

    return inline_fitness(pdb_id, **kwargs)

if __name__ == '__main__':
    sys.exit(main())
