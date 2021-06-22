


with open('E6-aln.txt', 'r') as f:      # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< 1 de 3
    raw_file = f.read()

raw_file_lines = [line for line in raw_file.split('\n') if line != '']

fasta_dict = {}

for line in raw_file_lines:
    if line[0] == '>':
        actual_key = line
        fasta_dict[actual_key] = ''
        
    else:
        fasta_dict[actual_key] =  fasta_dict[actual_key] + line




# En esta parte vamos a calcular los parámetros que necesitamos para 
# dar el formato necesario a las temás keys

# Suponiendo que sabemos de antemano el nombre del ID del epitope
# En el caso ejemplo es >E1


seq_epitope = fasta_dict['>E6']   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< 2 de 3

# Número de elementos antes la cadena

counter = 0
while seq_epitope[counter] == '-':
    counter += 1
start_repeted = counter

counter = -1
while seq_epitope[counter] == '-':
    counter -= 1
back_repeted = counter + 1


final_file = ''
for key,value in fasta_dict.items():
    final_file = final_file + key + '\n' + value[start_repeted:back_repeted] + '\n'
    # print(value[start_repeted:back_repeted])



with open('E6_aln-cuted.txt', 'w') as f:   # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< 2 de 3
    f.write(final_file)