import sys
import pandas as pd
import chemoUtils as cu

#Get commmand line arguments
drug_file = sys.argv[1]
target_file = sys.argv[2]
protein_nodes_file = sys.argv[3]

#Default values
cutoff = 0.50
n_iterations = 100
random_seed = 214
p_cutoff = .05

#Put data files into pandas dataframe
drug_df = pd.read_csv(drug_file)
target_df = pd.read_csv(target_file)
protein_nodes_df = pd.read_csv(protein_nodes_file)

#Sort protein_nodes_df by uniprot_accession
protein_nodes_df = protein_nodes_df.sort_values(by='uniprot_accession').reset_index(drop=True)

#Open output files SIF and name/indications nodeattr to write to
sif_file = open('network.sif', 'w')
name_file = open('name.nodeAttr', 'w')
indication_file = open('indication.nodeAttr', 'w')
name_file.write('name\n')
indication_file.write('indication\n')

#Loop over proteins in protein_nodes.csv
n_proteins = len(protein_nodes_df.index)
maccs_dict = cu.create_maccs_dict(drug_df)
for a in range(n_proteins):
    protein_A = protein_nodes_df.ix[a, 'uniprot_accession']

    for b in range(a + 1, n_proteins):
        protein_B = protein_nodes_df.ix[b, 'uniprot_accession']
        p_bootstrap = cu.compute_p_bootstrap(maccs_dict, protein_A, protein_B, drug_df, target_df, cutoff, n_iterations, random_seed)

        if p_bootstrap <= p_cutoff:
            sif_string = protein_A + " edge " + protein_B + '\n'
            sif_file.write(sif_string)

    name_string = protein_A + " = " + protein_nodes_df.ix[a, 'uniprot_id'] + '\n'
    indication_string = protein_A + " = " + protein_nodes_df.ix[a, 'indications'] + '\n'
    name_file.write(name_string)
    indication_file.write(indication_string)


#Close output files
sif_file.close()
name_file.close()
indication_file.close()