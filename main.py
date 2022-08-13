from fun_file import *

input_protein = 'rash'
primary_list = related_family(input_protein)
uniprot_input_protein = uniprot_mapping('ACC+ID', 'ACC+ID', primary_list[0][0])
clear_folder(os.getcwd() + '/pdb_files')
for i in primary_list[2]:
    pdb_structure_alphafold(i)
foldseek(uniprot_input_protein)
alignment_matrix = foldseek_alig_matrix('%s.m8' % uniprot_input_protein, uniprot_input_protein)
percentage_of_alignment = alignment_matrix[2]
# biomuta_profile = biomuta_resorce(primary_list[2], primary_list[3])
# biomuta_profile.to_pickle("./input_biomuta_%s.pkl" % uniprot_input_protein)
biomuta_profile = pd.read_pickle("./input_biomuta_%s.pkl" % uniprot_input_protein)
# hotspot_profile = hotspot_profiler(primary_list[2], primary_list[3])
# hotspot_profile.to_pickle("./input_hotspot_%s.pkl" % uniprot_input_protein)
hotspot_profile = pd.read_pickle("./input_hotspot_%s.pkl" % uniprot_input_protein)
# alignment_matrix = foldseek_alig_matrix('%s.m8' % uniprot_input_protein)
result = result_matrix(uniprot_input_protein, alignment_matrix, biomuta_profile, hotspot_profile, primary_list)

# # uncomment to see results of loocv for pfmi3dsa
"""
# result[3].to_pickle("./score_matrix_%s.pkl" % uniprot_input_protein)
# result[4].to_pickle("./probability_matrix_%s.pkl" % uniprot_input_protein)
loocv_count = loocv_input_protein (primary_list, uniprot_input_protein, result)
"""
# # uncomment to draw the RIN graph
"""
rin_protein(uniprot_input_protein, result[2])
plot_rin(uniprot_input_protein) # using R
"""
print(result)

# qcov (Fraction of query sequence covered by alignment)
# tcov (Fraction of target sequence covered by alignment)
# qseq Query sequence
# tseq Target sequence


# ass = pd.read_pickle("req_cont_mat rhoa 10a.pkl")
# # peines = pd.read_pickle("req_dist_mat rhoa 15a.pkl")
# p01112_mut = [122, 124, 169]
# ass = ass.replace({True: 1, False: 0})
# ass.index += 1
# ass.columns += 1
# data = []
# for column in ass:
#         if column in p01112_mut:
#             pass
#         else:
#             ass[column] = 0
# ass.to_csv(os.getcwd() + '/P61586_PYT.csv', sep='\t')
print('hi')
