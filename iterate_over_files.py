import os

#inititate varible for current working dcirectory
cwd = os.getcwd()

#inititate variable for sample folder
sample_folder = '\samples'

#create final path to samples folder
path_samples = cwd+sample_folder

#initiate list with sample names
samples = []

# iterate over samples in sample folder 
for filename in os.listdir(path_samples):
    # use only .fastq-files
    if filename.endswith(".fastq"):
        samples.append(filename)
    else:
        continue