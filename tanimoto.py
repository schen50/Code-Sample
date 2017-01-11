import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import chemoUtils as cu


# Input: a db_id and pandas dataframe for targets.csv
# Output:
def find_targets(db_id):
    return targets_data.loc[targets_data['db_id'] == db_id, 'uniprot_accession'].tolist()


############## MAIN FUNCTION #########################
# Get arguments. Put input data into pandas dataframes
drugs_data = pd.read_csv(sys.argv[1])
targets_data = pd.read_csv(sys.argv[2])
output_file = sys.argv[3]

# Sort drugs_data by db_id
drugs_data = drugs_data.sort_values(by='db_id').reset_index(drop=True)

# Create dict mapping a drug's db_id to its target uniprot_accessions
targets_dict = {}
for index, db_id in drugs_data['db_id'].iteritems():
    targets_dict[db_id] = find_targets(db_id)

# Data structures to store info for output
drug1_vec = []
drug2_vec = []
score_vec = []
share_vec = []

# FOR HISTOGRAMS
"""all_tanimoto_vec = []
shared_tanimoto_vec = []
notshared_tanimoto_vec = []"""

for drug_1 in range(len(drugs_data.index)):
    d1_id = drugs_data.ix[drug_1, 'db_id']
    targets_1_set = set(targets_dict[d1_id])

    for drug_2 in range(drug_1 + 1, len(drugs_data.index)):
        d2_id = drugs_data.ix[drug_2, 'db_id']
        targets_2_set = set(targets_dict[d2_id])

        #Compute tanimoto between drug_1 and drug_2
        tanimoto_score = cu.tanimoto(drugs_data.ix[drug_1, 'maccs'].split(), drugs_data.ix[drug_2, 'maccs'].split())
        formatted_tanimoto_score = format(tanimoto_score, '.6f')

        #See if drug_1 and drug_2 share a target
        share = 0
        intersection = targets_1_set.intersection(targets_2_set)
        if len(intersection) > 0:
            share = 1

        #Append results for output file
        drug1_vec.append(d1_id)
        drug2_vec.append(d2_id)
        score_vec.append(formatted_tanimoto_score)
        share_vec.append(share)

        #Append results for histograms
        """all_tanimoto_vec.append(tanimoto_score)
        if share == 1:
            shared_tanimoto_vec.append(tanimoto_score)
        else:
            notshared_tanimoto_vec.append(tanimoto_score)"""

result = np.column_stack((drug1_vec, drug2_vec, score_vec, share_vec))
column_names = ['drug_1', 'drug_2', 'tanimoto', 'share_target']
result_df = pd.DataFrame(data=result, columns=column_names)
result_df.to_csv(output_file, header=False, index=False)

#Generate histograms
"""plt.hist(all_tanimoto_vec)
plt.xlabel('Tanimoto Values')
plt.ylabel('Frequency')
plt.title('schen50 All')
plt.show()

plt.hist(shared_tanimoto_vec)
plt.xlabel('Tanimoto Values')
plt.ylabel('Frequency')
plt.title('schen50 Shared')
plt.show()

plt.hist(notshared_tanimoto_vec)
plt.xlabel('Tanimoto Values')
plt.ylabel('Frequency')
plt.title('schen50 Not Shared')
plt.show()"""
