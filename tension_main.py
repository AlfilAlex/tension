from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import is_aa

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

    # Si model solo tiene una cadena, entonces se considera que solo
    # incluye el epitope. Sin embargo es necesario determinar el nombre.
    # Como este programa aún no determina el nombre, hace el intento con
    # "A" y "C". Sin embargo, otro nombre no sería reconocido.
        chain_list = model_1.get_list()  # <--- Esto es una lista de objetos chain "<Chain id=C>"
        chainItem = model_1[chain_list[0].get_id()]
        chain = solvent_delete(chainItem)

    # A demás debe avisar si se está utilizando una cadena muy extensa
    # o al menos avisar de su tamaño.
        if len(chain_list) >= 15:
            print('CUIDADO, la cadena es muy extensa')






    
    
    else:
        chain_list = model_1.get_list()
        min_len = 100_000_000
        for chainItem in chain_list:

            SF_chainItem = solvent_delete(chainItem)
            chain_length = len(SF_chainItem)

            if chain_length < min_len:
                chain = SF_chainItem
                min_len= chain_length
    

    # Una vez que se cuenta con el chain correcto, se pueden seleccionar
    # Los resiudos iniciales y finales (residue_0, residue_f)
   # print(chain.get_list())
    residue_0 = chain[0]
    chain_length = len(chain)
    residue_f = chain[chain_length - 1] # Ahora trabajamos con una ista que tiene indice 0.


    # Por la documentación se sabe que esta operación entre elementos de
    # la clase residuos es posible.
    d_ext = residue_0['CA'] - residue_f['CA']

    tension = tension_betarelativa(d_ext, chain_length)

    return round(tension, 5), round(d_ext, 3)

# Estaría padre generar una función que pueda calcular el angulo que hay entre
# El vector que va de un extremo del epitope al otro y el plano formado por la 
# superficie relativa de la HLA.
def angule(chain_epitope, chain_hla):
    # Cómo rayos encontramos una superficie? 
    pass


def solvent_delete(pdbChain):
    SF_chainItem = [residue for residue in pdbChain if is_aa(residue)]
    #SF_chainItem = [residue for residue in pdbChain if "W" not in residue.get_id()]

    return SF_chainItem
    



# Esto es un ejemplo de cómo usar las funciones. Básicamente el  
# if __name__ == "__main__"... permite que el programa se ejecute, 
# ya que es una condición que se va a cumplir.

if __name__ == "__main__":
    structure_id = "ejemplo_1"
    filename_1 = "6uon.pdb"
    tension_1, d_ext_1 = epitope_tension(structure_id, filename_1)

    print(f'La distancia entre los carbonos alfa del primer y último AA es: {d_ext_1} Amstrongs')
    print(f'La tensión beta relativa del epitope en la HLA es de {tension_1*100} %')




    