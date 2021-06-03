from Bio.PDB.PDBParser import PDBParser

def tension_betarelativa(distance, chain_length, AC_MD):
    """ This function calculate the relative tension of a AA chain
        given that the mean distance of a alpha carbon between to  
        continues AA residues is 3.5 A as it can be found in Stryer. 
    """

    t = distance/(AC_MD*chain_length)  # AC_MD: Alpha carbon mean distance.
    return t


def epitope_tension(structure_id, filename):

    parser = PDBParser()
    structure = parser.get_structure(structure_id, filename)

    model_1 = structure[0]

    chain = model_1['A']

    residue_0 = chain[1]
    chain_length = len(chain)
    residue_f = chain[chain_length]

    d_ext = residue_0['CA'] - residue_f['CA']

    tension = tension_betarelativa(d_ext, chain_length)
    return tension, d_ext

# An exaple of how to use the program

if __name__ == "__main__":
    structure_id = "ejemplo_1"
    filename_1 = "lll.pdb"
    tension_1, d_ext_1 = epitope_tension(structure_id, filename_1)

    filename_2 = "lmf.pdb"
    tension_2, d_ext_2 = epitope_tension(structure_id, filename_2)

    print(f'La distancia entre los carbonos alfa del primer y último AA es: {d_ext_1} Amstrongs')
    print(f'La distancia entre los carbonos alfa del primer y último AA es: {d_ext_2} Amstrongs')

    print(f'La tensión beta relativa del epitope en la HLA es de {tension_1}')
    print(f'La tensión beta relativa del epitope en la HLA es de {tension_2}')




    