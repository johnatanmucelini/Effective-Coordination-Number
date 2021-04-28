import sys
import numpy as np
import ase.io
import argparse
from scipy.spatial.distance import cdist

def pbcspreadrowdata(pbc, data, shell=1, seconddimshepe=-1):
    refdata = []
    ix, jx, kx = np.array(pbc, dtype=int)*shell
    for i in range(-ix, ix+1):
        for j in range(-jx, jx+1):
            for k in range(-kx, kx+1):
                refdata.append(data)
    refdata = np.array(refdata).reshape(seconddimshepe)
    return refdata

# Variables ref<atomsdata> are equal to <atomsdata> for molecules,
# in contrast, for atomic structures with any PBC (periodic
# boundery conditions) they are equal to <atomsdata> for the atoms
# in a hypercell of the atomic structure input, larger in the axis
# with PBC


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='The script mensure the Effective Coordination Number (ECN) and the average bound distance (dav), for one or more atomic structures. It work with any atomic structuresfile readble by ase.io.read(), including structures with periodic boundary conditions.')
    parser.add_argument('--files', nargs='+', help='the molecules (xyz, geometry.in, etc) to analyze.')
    parser.add_argument('--json_file', metavar='file.json', action='store', default=None, help='the name of a json file to save the data.')
    args = parser.parse_args()

    # creating variables to save a json file
    if args.json_file:
        import pandas as pd
        meus_dados = pd.DataFrame()
        all_positions = []
        all_chemical_symbols = []
        all_ecn = []
        all_dav = []

    for mfile in args.files:
        print("Reading file: " + mfile)
        structure = ase.io.read(mfile)

        # Getting atomi information
        pbc = structure.pbc
        positions = np.array(structure.positions)
        positions = positions - positions.mean(axis=0)
        cheme = np.array(structure.get_chemical_symbols())
        refcheme = pbcspreadrowdata(pbc, cheme)
        pbc = structure.pbc
        if np.any(pbc):
            refpositions = []
            ix, jx, kx = np.array(structure.pbc, dtype=int)
            ai, aj, ak = structure.cell.array
            for i in range(-ix, ix+1):
                for j in range(-jx, jx+1):
                    for k in range(-kx, kx+1):
                        refpositions.append(positions + ai*i + aj*j + ak*k)
            refpositions = np.array(refpositions).reshape(-1,3)
        else:
            refpositions = positions

        dij = cdist(positions, refpositions) 
        dij = dij + 100*np.array(dij < 0.001, dtype=int)
        refdij = pbcspreadrowdata(pbc, dij, seconddimshepe=(-1,len(dij[0])))

        dav = np.max(cdist(positions, refpositions), axis=1)
        refdav = pbcspreadrowdata(pbc, dav)

        ecn = np.zeros_like(dav)
        refecn = pbcspreadrowdata(pbc, ecn)

        refdav_pre = np.zeros_like(dav)
        refecn_pre = np.zeros_like(dav)

        step = 0
        conv_val = 1E-8
        print('Convergence: ' + str(conv_val) + '\n   dav         ECN')

        # Initializing self consistency
        while step < 2 or delta_ecn > conv_val or delta_dav > conv_val:
            refdav_pre = refdav * 1.
            refecn_pre = refecn * 1.
            refpij = np.exp(1 - (2*refdij/(refdav.reshape(-1, 1) + refdav.reshape(1, -1)))**6)

            refecn = np.sum(refpij, axis=1)
            refdav = np.sum(refpij * refdij, axis=1) / np.sum(refpij, axis=1)
            refecn = np.sum(refpij, axis=1)

            delta_ecn = np.average(np.abs(refecn_pre - refecn))
            delta_dav = np.average(np.abs(refdav_pre - refdav)) 
            print("%8.3E   %8.3E" %(delta_ecn, delta_dav))
            step += 1
    
        print("The convergence criteria has reached.")

        print('index element ecn   dav')
        for index, element in enumerate(cheme):
            print("%3d     %2s    %4.2f  %4.2f" % (index, element, refecn[index], refdav[index]))
        print('')

        if args.json_file:
            all_positions.append(positions)
            all_chemical_symbols.append(cheme)
            all_ecn.append(refecn[:len(ecn)])
            all_dav.append(refdav[:len(ecn)])

    if args.json_file:
        meus_dados['files'] = args.files
        meus_dados['positions'] = all_positions
        meus_dados['cheme'] = all_chemical_symbols
        meus_dados['ecn'] = all_ecn
        meus_dados['dav'] =  all_dav
        
        print('Writing json_file: {}'.format(args.json_file))
        meus_dados.to_json(args.json_file)

        