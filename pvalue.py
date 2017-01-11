import argparse
import pandas as pd
import chemoUtils as cu

#Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("-n", action="store", dest="iterations", type=int, default=100)
parser.add_argument("-r", action="store", dest="random_seed", type=int, default=214)
parser.add_argument("drug_file")
parser.add_argument("target_file")
parser.add_argument("protein_A")
parser.add_argument("protein_B")
args = parser.parse_args()
n_iterations = args.iterations
random_seed = args.random_seed
drug_file = args.drug_file
target_file = args.target_file
protein_A = args.protein_A
protein_B = args.protein_B
cutoff = 0.5

#Put data into pandas dataframe
drug_df = pd.read_csv(drug_file)
target_df = pd.read_csv(target_file)
maccs_dict = cu.create_maccs_dict(drug_df)

#Compute pbootstrap and print
p_bootstrap = cu.compute_p_bootstrap(maccs_dict, protein_A, protein_B, drug_df, target_df, cutoff, n_iterations, random_seed)
print(p_bootstrap)
