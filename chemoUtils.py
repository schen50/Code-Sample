import random

# Inputs: Two maccs fingerprint arrays
# Output: The Tanimoto Coefficient as a float
def tanimoto (f1, f2):
    s1 = set(f1)
    s2 = set(f2)
    intersection = s1.intersection(s2)
    return float(len(intersection))/(len(f1) + len(f2) - len(intersection))

#Input: pandas drugs.csv dataframe
#Output: dict mapping index in df to an array representation of its maccs fingerprint
def create_maccs_dict(drug_df):
    maccs_dict = drug_df['maccs'].to_dict()
    return {key: list(map(int, value.split())) for key, value in maccs_dict.items()}


#Input: A protein swissprot id ("uniprot_accession" and pandas drugs.csv/targets.csv dataframe
#Output: List of indices corresponding to drug_df of the binding drugs
def find_binding_drugs(protein_id, target_df, drug_df):
    drug_ids =  target_df.loc[target_df['uniprot_accession'] == protein_id, 'db_id']
    return drug_df[drug_df['db_id'].isin(drug_ids)].index


#Input: Pandas drug.csv dataframe, two lists of row indices corresponding to drug_df, cutoff for tanimoto score
#Output: T_sum
#Computes pairwise Tanimoto similarities and sums values greater than cutoff
def compute_T_sum(maccs_dict, rows_a, rows_b, cutoff):
    T_sum = 0
    for index_a in rows_a:
        maccs_a = maccs_dict[index_a]
        for index_b in rows_b:
            maccs_b = maccs_dict[index_b]
            tanimoto_score = tanimoto(maccs_a, maccs_b)
            if tanimoto_score > cutoff:
                T_sum += tanimoto_score
    return T_sum

#Input: Maccs dict, two protein swissprot/uniprot_accession ids, pandas drug df and target df, cutoff for computing T_sum, number of iterations, and random seed
#Output: PBootstrap
#Calls: compute_T_sum, find_binding_drugs and tanimoto
def compute_p_bootstrap(maccs_dict, protein_A, protein_B, drug_df, target_df, cutoff, n_iterations, random_seed):
    # Create row indices vector of the binding drugs
    rows_a = find_binding_drugs(protein_A, target_df, drug_df)
    rows_b = find_binding_drugs(protein_B, target_df, drug_df)

    # Compute T_sum using actual ligands
    T_sum = compute_T_sum(maccs_dict, rows_a, rows_b, cutoff)

    # Repeat n times with random ligand sets
    random.seed(random_seed)
    na = rows_a.size
    nb = rows_b.size

    n_greater = 0
    for i in range(n_iterations):
        rows_a = random.sample(drug_df.index, na)
        rows_b = random.sample(drug_df.index, nb)
        rand_T_sum = compute_T_sum(maccs_dict, rows_a, rows_b, cutoff)
        if rand_T_sum >= T_sum:
            n_greater += 1

    #Return p_bootstrap
    return float(n_greater) / n_iterations
