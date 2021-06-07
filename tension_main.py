from Bio.PDB.PDBParser import PDBParser

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

    # A demás debe avisar si se está utilizando una cadena muy extensa
    # o al menos avisar de su tamaño.

    # Esto esta mal, es  decir, no es la forma correcta de usar el 
    # try y el except, pero es la unica forma que se me ocurre por
    # ahora.
        try:
            chain = model_1['A']
        except:
            chain = model_1['C']
    
    
    else:
        chain_list = model_1.get_list()
        max_len = 1
        for chainItem in chain_list:
            # Probablemente aquí deba agreagar alguna función 
            # "is_solvent" o algo así, que verifique que la
            #  cadena no está contando el solvente. O algo
            # que lo elimine.
            chain_length = len(chainItem)
            if chain_length > max_len:
                chain = model_1[chainItem]
                max_len= chain_length
    # Aquí deberiamos análizar todas las cadenas que hay en model 1 
    # Y entonces escoger la cadena que tenga de 9 a 12 aminoácidos.
        


    # Una vez que se cuenta con el chain correcto, se pueden seleccionar
    # Los resiudos iniciales y finales (residue_0, residue_f)
    residue_0 = chain[1]
    chain_length = len(chain)
    residue_f = chain[chain_length]


    # Por la documentación se sabe que esta operación entre elementos de
    # la clase residuos es posible.
    d_ext = residue_0['CA'] - residue_f['CA']

    tension = tension_betarelativa(d_ext, chain_length)

    return tension, d_ext

# Estaría padre generar una función que pueda calcular el angulo que hay entre
# El vector que va de un extremo del epitope al otro y el plano formado por la 
# superficie relativa de la HLA.
def angule(chain_epitope, chain_hla):
    # Cómo rayos encontramos una superficie? 
    pass


# Esto es un ejemplo de cómo usar las funciones. Básicamente el  
# if __name__ == "__main__"... permite que el programa se ejecute, 
# ya que es una condición que se va a cumplir.

if __name__ == "__main__":
    structure_id = "ejemplo_1"
    filename_1 = "lmf.pdb"
    tension_1, d_ext_1 = epitope_tension(structure_id, filename_1)

    filename_2 = "5c0g.pdb"
    tension_2, d_ext_2 = epitope_tension(structure_id, filename_2)

    print(f'La distancia entre los carbonos alfa del primer y último AA es: {d_ext_1} Amstrongs')
    print(f'La distancia entre los carbonos alfa del primer y último AA es: {d_ext_2} Amstrongs')

    print(f'La tensión beta relativa del epitope en la HLA es de {tension_1*100} %')
    print(f'La tensión beta relativa del epitope en la HLA es de {tension_2*100} %')




    