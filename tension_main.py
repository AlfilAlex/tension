from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import is_aa
import numpy as np

def tension_betarelativa(distance, chain_length, AC_MD = 3.5):
    """ This function calculate the relative tension of a AA chain given that 
    the mean distance of a alpha carbon between to continues AA residues is 3.5 A 
    as it can be found in Stryer.
    """

    t = distance/(AC_MD*chain_length)  # AC_MD: Alpha carbon mean distance.
    return t


def epitope_tension(structure_id, filename):

    parser = PDBParser()
    structure = parser.get_structure(structure_id, filename)

    model_1 = structure[0]


    if len(model_1) == 1:
        # Es necesario determinar el nombre.
        chain_list = model_1.get_list()  # <--- Esto es una lista de objetos chain "<Chain id=C>"
        chainItem = model_1[chain_list[0].get_id()]
        chain = solvent_delete(chainItem)

        if len(chain_list) >= 16:
            print('CUIDADO, la cadena es muy extensa')
   
    else:
        chain_list = model_1.get_list()
        min_len = 100_000_000
        max_len = 1
        for chainItem in chain_list:

            SF_chainItem = solvent_delete(chainItem)
            chain_length = len(SF_chainItem)

            if chain_length < min_len:
                chain = SF_chainItem
                min_len = chain_length
            if chain_length > max_len:
                chain_max = SF_chainItem
                max_len = chain_length
    

    # Los resiudos iniciales y finales (residue_0, residue_f)
    residue_0 = chain[0]
    chain_length = len(chain)
    residue_f = chain[chain_length - 1] # Ahora trabajamos con una lista que tiene indice 0.
    
    matriz_angules = angule_change(chain)
    
    # Por la documentación se sabe que esta operación entre elementos de
    # la clase residuos es posible.
    d_ext = residue_0['CA'] - residue_f['CA']
    tension = tension_betarelativa(d_ext, chain_length)

    return round(tension, 5), round(d_ext, 3), matriz_angules


def solvent_delete(pdbChain):
    SF_chainItem = [residue for residue in pdbChain if is_aa(residue)] # [residue for residue in pdbChain if "W" not in residue.get_id()]
    return SF_chainItem

# Estaría padre generar una función que pueda calcular el angulo que hay entre
# El vector que va de un extremo del epitope al otro y el plano formado por la 
# superficie relativa de la HLA.
def angule(chain_epitope, chain_hla):
    # Cómo rayos encontramos una superficie? 
    pass  

# También una forma de medir la distribución de los angulos en el epitope
def angule_change(chain):
    """ This function takes the epitope as input and first calculates the vector
         that goes from one extreme to another. Later for each of the carbons
         subsequent alphas, calculates the angles from the first alpha carbon
         from the epitope to the next alpha carbon until the last one, whose
         angle must be 0.
    """
    # Chain es una lista de residuos y ya no un objeto chain de BIO.PDB
    coord_CA0 = chain[0]['CA'].get_coord()
    coord_CAf = chain[-1]['CA'].get_coord()
    base_vector = coord_CAf - coord_CA0
    print('el tipo de objeto es: ' + str(type(base_vector)))

    angule_matrix = [0] # De el primer AA
    for i in range(1, len(chain)):
        base_i = chain[i]['CA'].get_coord() - coord_CA0 
        dot_product = np.dot(base_vector, base_i)
        module_mult = np.linalg.norm(base_i)*np.linalg.norm(base_vector)

        angule = np.arccos(dot_product/module_mult)
        angule_matrix.append(np.degrees(angule))

    print(angule_matrix)
    return angule_matrix


# Esto es un ejemplo de cómo usar las funciones. Básicamente el  
# if __name__ == "__main__"... permite que el programa se ejecute, 
# ya que es una condición que se va a cumplir.

if __name__ == "__main__":
    structure_id = "ejemplo_1"
    filename_1 = "lmf.pdb"
    tension_1, d_ext_1, angule = epitope_tension(structure_id, filename_1)

    print(f'The distance between the alpha carbons of the first and last AA is: {d_ext_1} A')
    print(f'The relative beta stress of the epitope in HALLA is {tension_1*100} %')
    print(f'The distribution of the angles is: {angule} ')




    