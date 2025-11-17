from Bio.Affy import CelFile

# Define the path to your .CEL file
cel_file_path = '/Users/jamesorr/Documents/Imperial/Project_1/AlphaGenome_IBD/Gene_expression/Data/Bp040_410329M151.CEL'

# Open and parse the .CEL file
with open(cel_file_path, 'rb') as file:
    cel_data = CelFile.read(file)

# Says that this is not a cell file