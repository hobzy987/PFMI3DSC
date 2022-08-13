"""
Created on Sat May 21 10:05:42 2022

@author: Hoballa
"""

# Library file
import os
import string
import urllib.parse
import urllib.request
import pandas as pd
import copy
import shutil
import pcmap
from collections import Counter
from itertools import chain
from difflib import SequenceMatcher
import Levenshtein

def uniprot_mapping(fromtype, totype, identifier):
    """Takes an identifier, and types of identifier
    (to and from), and calls the UniProt mapping service"""
    base = 'http://www.uniprot.org'
    tool = 'mapping'
    params = {'from': fromtype,
              'to': totype,
              'format': 'tab',
              'query': identifier,
              }
    # urllib turns the dictionary params into an encoded url suffix
    data = urllib.parse.urlencode(params)
    # construct the UniProt URL
    url = base + '/' + tool + '?' + data
    # and grab the mapping
    response = urllib.request.urlopen(url).read().decode()
    id = response.split('\t')
    # response.read() provides tab-delimited output of the mapping
    return id[-1].split()[0]


# get related family of protein and UniProt human members
def related_family(p):
    global temp_family, protein_family, protein_list
    protein_name = ' ' + p.upper() + '_HUMAN'
    my_file = open("similar.txt", 'r')
    my_file.read()
    with open("similar.txt", 'r') as file:
        x = file.readlines()
    for i in range(len(x)):
        if ',' not in x[i]:
            temp_family = x[i]
        if protein_name in x[i]:
            protein_family = temp_family
            r = x.index(protein_family)
            for j in range(r + 1, len(x)):
                if ',' not in x[j]:
                    break
                f = str(x[r + 1:j + 1]).split(',')
                protein_list = [word.strip(string.punctuation) for word in f]
    protein_family = protein_family.split('\n')[0]
    human_protein_list = []
    for protein in protein_list:
        # protein = protein.split()
        if '_HUMAN' in protein:
            human_protein_list.append(protein)
    # Index_Input_Protein = [human_protein_list.index(i) for i in human_protein_list if protein_name in i]
    # Index_Input_Protein = [int(i) for i in Index_Input_Protein]
    uniprot_list = []
    for i in range(len(human_protein_list)):
        uniprot_list.append(human_protein_list[i][human_protein_list[i].find('(') + 1:human_protein_list[i].find(')')])
    gene_name = []
    for i in uniprot_list:
        gene_id = uniprot_mapping('ACC+ID', 'GENENAME', i)
        gene_name.append(gene_id)

    return [protein_name.split(), [protein_family], uniprot_list, gene_name]


# get PDB file from Alphafold protein structure database

def pdb_structure_alphafold(p):
    w = os.getcwd() + '/pdb_files/{}.pdb'.format(p)
    url = "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v2.pdb" % p
    urllib.request.urlretrieve(url, w)


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def union(lst1, lst2):
    final_list = list(set(lst1) | set(lst2))
    return final_list


def remover(l):
    newlist = []
    for i in l:
        if isinstance(i, list) and len(i) != 0:
            newlist.append(remover(i))
        # if not a list
        if not isinstance(i, list):
            newlist.append(i)
    return newlist


def remove_nestings(l):
    output = []
    for i in range(len(l)):
        if len(l[i]) == 1:
            output.append(l[i][0])
        else:
            output.append(l[i])
    return output


def removal_index(string, list_of_index):
    list_of_index = list(sorted(list_of_index))
    sep_str = []

    if max(list_of_index) < len(string):
        sep_str.append(string[0:list_of_index[0]])
        for z in range(len(list_of_index)):
            z = z - 1
            sep_str.append(string[list_of_index[z] + 1:list_of_index[z + 1]])
        sep_str.append(string[list_of_index[-1] + 1:])
        list_of_index.pop(-1)
        result = ''.join(sep_str)
        return result


def biomuta_resorce(uniprot_set, gene_name):
    req_cols = ['uniprot_ac', 'uniprot_pos', 'ref_aa', 'alt_aa', 'do_name']
    df = pd.read_csv("biomuta-master.csv", usecols=req_cols)
    temp = []
    for j in uniprot_set:
        contain_values = df[df['uniprot_ac'].str.contains(j)]
        s = contain_values['uniprot_pos']
        w = sorted(s.unique().tolist())
        temp.append(w)
    mutation_profile = pd.DataFrame(list(zip(uniprot_set, gene_name, temp)),
                                    columns=['uniprot_ID', 'genes', 'biomuta_profile'])
    return mutation_profile


def hotspot_profiler(uniprot_set, gene_name_list):
    hotspot_profile1 = []
    dhotspots = pd.read_csv("3dhotspots.csv")
    cancerhotspots = pd.read_csv("cancerhotspots.csv")
    dhotspots_list = []
    cancerhotspots_list = []
    for gene in gene_name_list:
        value = dhotspots[dhotspots['Gene'].str.contains(gene)]
        s = value['Residue']
        dhotspots_list.append(s.tolist())
        value = cancerhotspots[cancerhotspots['Gene'].str.contains(gene)]
        t = value['Residue']
        cancerhotspots_list.append(t.tolist())
    for i in range(len(gene_name_list)):
        onion = union(dhotspots_list[i], cancerhotspots_list[i])
        hotspot_i = []
        for j in onion:
            s = j.split('-')
            if len(s) == 1:
                n = int(''.join(filter(str.isdigit, j)))
                hotspot_i.append(n)
            else:
                p = int(''.join(filter(str.isdigit, s[0])))
                q = int(''.join(filter(str.isdigit, s[1])))
                l = [i for i in range(p, q + 1)]
                hotspot_i.extend(l)
        hotspot_profile1.append(sorted(list(set(hotspot_i))))
    profile_last = {'uniprot_ID': uniprot_set, 'genes': gene_name_list, 'hotspot_profile': hotspot_profile1}
    hotspot_profile = pd.DataFrame(profile_last)
    return hotspot_profile


def String_character_counter(string, start):
    string = str(string)
    start = int(start)
    counter_list = []
    for char in string:
        if char.isalpha() is True:
            counter_list.append(start)
            start = start + 1
        if not char.isalpha():
            counter_list.append(0)

    return counter_list


# def seq_sim(seq1, seq2):
#     length_of_shorter_seq = min(len(seq1), len(seq2))
#     shared = 0
#     different = 0
#     for i in range(length_of_shorter_seq):
#         if seq1[i] == seq2[i]:
#             shared += 1
#         else:
#             different += 1
#
#     sequence_identity = shared / length_of_shorter_seq
#     x = round(sequence_identity, 2)
#     return x



def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()


def foldseek_alig_matrix(alignment_file, uniprot_input_protein):
    s = pd.read_csv(alignment_file, sep='\t', header=None)
    s.columns = ['query', 'target', 'qaln', 'taln', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'qlen', 'tlen',
                 'tcov', 'qseq', 'tseq']
    get_indexes = lambda x, xs: [i for (y, i) in zip(xs, range(len(xs))) if x == y]
    query_seq_str = []
    target_seq_str = []
    query_seq_list = []
    target_seq_list = []
    start_query = []
    start_target = []
    end_query = []
    end_target = []
    counter_of_char_query = []
    counter_of_char_target = []
    index_to_remove = []
    protein_names = []
    protein_lengths = []
    seq_sim = []
    for i in s.iloc():
        q = [sub for sub in i['qaln']]
        t = [sub for sub in i['taln']]
        query_seq_str.append(i['qaln'])
        target_seq_str.append(i['taln'])
        query_seq_list.append(q)
        target_seq_list.append(t)
        start_query.append(i['qstart'])
        start_target.append(i['tstart'])
        end_query.append((i['qend']))
        end_target.append(i['tend'])
        protein_names.append(i['target'].split('.')[0])
        protein_lengths.append(i['tlen'])
    # percentage_of_alignment = s.set_index('target').to_dict()['tcov']
    s['target'] = protein_names
    percentage_of_alignment = s[['target', 'tcov']]
    for index, row in s.iterrows():
        seq1 = row['qseq']
        seq2 = row['tseq']
        sim = Levenshtein.ratio(seq1, seq2) * 100
        seq_sim.append(sim)
    percentage_of_alignment['seq_sim'] = seq_sim
    percentage_of_alignment.to_csv(os.getcwd() + '/' + f'{uniprot_input_protein} percentage of alignment' + '.csv', sep='\t')
    for i in range(len(start_query)):
        w = String_character_counter(query_seq_str[i], start_query[i])
        x = String_character_counter(target_seq_str[i], start_target[i])
        counter_of_char_query.append(w)
        counter_of_char_target.append(x)
        s = get_indexes(0, w)
        index_to_remove.append(s)
    for j in range(len(counter_of_char_query)):
        if counter_of_char_query[j][0] != 1:
            n = counter_of_char_query[j][0]
            f = list(range(1, n))
            counter_of_char_query[j] = f + counter_of_char_query[j]

            r = [0] * len(f)
            seq_add = ['-'] * len(f)
            counter_of_char_target[j] = r + counter_of_char_target[j]
            target_seq_list[j] = seq_add + target_seq_list[j]

    for j in range(len(index_to_remove)):
        if len(index_to_remove[j]) > 0:
            for h in index_to_remove[j]:
                counter_of_char_query[j] = counter_of_char_query[j][:h] + counter_of_char_query[j][h + 1:]
                counter_of_char_target[j] = counter_of_char_target[j][:h] + counter_of_char_target[j][h + 1:]
                query_seq_list[j] = query_seq_list[j][:h] + query_seq_list[j][h + 1:]
                target_seq_list[j] = target_seq_list[j][:h] + target_seq_list[j][h + 1:]
        else:
            continue
    for r in range(len(target_seq_list)):
        if len(target_seq_list[r]) < len(target_seq_list[0]):
            q = len(target_seq_list[0]) - len(target_seq_list[r])
            s_add = ['-'] * q
            count_add = [0] * q
            target_seq_list[r] = target_seq_list[r] + s_add
            counter_of_char_target[r] = counter_of_char_target[r] + count_add
    alignment_matrix_seq = pd.DataFrame(target_seq_list)
    alignment_matrix_seq = alignment_matrix_seq.transpose()
    alignment_matrix_seq.columns = protein_names
    alignment_matrix_res = pd.DataFrame(counter_of_char_target)
    alignment_matrix_res = alignment_matrix_res.transpose()
    alignment_matrix_seq.loc['protein_length'] = protein_lengths
    alignment_matrix_res.columns = protein_names
    return alignment_matrix_seq, alignment_matrix_res, percentage_of_alignment


def result_matrix(uniprot_input_protein, alinment_matrix, biomuta_profile, hotspot_profile, primary_list):
    residue_matrix = alinment_matrix[1]
    mutation_probabilty = {}
    hotspot_probabilty = {}
    nonmutation_probability = {}
    for column in residue_matrix:
        biomuta_mut = biomuta_profile[biomuta_profile['uniprot_ID'].isin([column])]['biomuta_profile'].tolist()[0]
        hotspot_mut = hotspot_profile[hotspot_profile['uniprot_ID'].isin([column])]['hotspot_profile'].tolist()[0]
        if column == uniprot_input_protein:
            pass
        else:
            for row in residue_matrix[column]:
                if row in biomuta_mut and row in hotspot_mut:
                    residue_matrix[column] = residue_matrix[column].replace([row], 2)
                if row in biomuta_mut and row not in hotspot_mut:
                    residue_matrix[column] = residue_matrix[column].replace([row], 1)
                if row in hotspot_mut and row not in biomuta_mut:
                    residue_matrix[column] = residue_matrix[column].replace([row], 2)
                if row not in hotspot_mut and row not in biomuta_mut:
                    residue_matrix[column] = residue_matrix[column].replace([row], 0)
            mutation_probabilty[column] = ((len(biomuta_mut)-len(hotspot_mut)) / alinment_matrix[0].loc['protein_length', column])#### nonhotspotmutation prob = #nonhot mutation / len protein
            hotspot_probabilty[column] = (len(hotspot_mut) / alinment_matrix[0].loc['protein_length', column])#### hotspot mutation probability = #hotspot/len protein
            ###non_mutation = 1-(prob_hotspot + prob_nonhotspot)
    uniprot_list = copy.deepcopy(primary_list[2])
    columns = [s for s in uniprot_list if s != uniprot_input_protein]
    scores_of_residues = residue_matrix[columns].sum(axis=1)
    alinment_matrix[0]['scores'] = scores_of_residues
    alinment_matrix[0].drop(index=alinment_matrix[0].index[-1], axis=0, inplace=True)
    score_matrix = copy.deepcopy(alinment_matrix[1])
    # alinment_matrix[0].to_pickle("./score_loocv_%s.pkl" %uniprot_input_protein)
    for column in residue_matrix:
        if column == uniprot_input_protein:
            pass
        else:
            c = mutation_probabilty.get(column)
            e = hotspot_probabilty.get(column)
            r = 1 - (c + e)
            # e = d
            if c == 0:
                residue_matrix.loc[residue_matrix[column] == 0, column] = 1
            else:
                residue_matrix.loc[residue_matrix[column] == 0, column] = r
                residue_matrix.loc[residue_matrix[column] == 1, column] = c
                residue_matrix.loc[residue_matrix[column] == 2, column] = e

    prob_matrix = residue_matrix[columns]
    probability_of_residues = prob_matrix[columns].prod(axis=1)
    alinment_matrix[0]['probability'] = probability_of_residues
    residue_matrix['probability'] = probability_of_residues
    matrix = alinment_matrix[0]
    probability_matrix = copy.deepcopy(alinment_matrix[1])
    matrix.insert(loc=0, column='residue_number', value=range(1, len(matrix) + 1))
    # matrix.drop(index=matrix.index[-1], axis=0, inplace=True)
    # residue_matrix.to_pickle("./probability_loocv_%s.pkl" % uniprot_input_protein)
    threshold = 0.01 / (len(primary_list[2]) - 1)
    threshold1 = len(primary_list[2]) / 2
    result = matrix.loc[
        (matrix['probability'] <= threshold) & (matrix['scores'] > threshold1)]
    result = result.reset_index(drop=True)
    result = result.iloc[:, [0, 1, -2, -1]]
    b = result['residue_number'].tolist()
    pd.set_option('display.float_format', '{:.3g}'.format)
    html1 = matrix.to_html()
    html2 = result.to_html()
    headers_html = "Input Protein : {},<br> Family of Input Protein : {},<br> Uniprot ID's of Family : {},<br> Gene Names of Family : {}<br>".format(
        primary_list[0], primary_list[1], primary_list[2], primary_list[3])
    text_file = open("%s.html" % uniprot_input_protein, "w")
    text_file.write(headers_html)
    text_file.write('matrix of alignment<br>')
    text_file.write(html1)
    text_file.write(
        'result residues,<br> thresholds = {} - {},<br> result residue : {}<br>'.format(threshold, threshold1, b))
    text_file.write(html2)
    text_file.close()
    return [matrix, result, b, score_matrix, probability_matrix]


def clear_folder(path):
    folder = path
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


def foldseek(uniprot_input_protein):
    s = 'export PATH=$(pwd)/foldseek/bin/:$PATH'
    x = "foldseek easy-search pdb_files/{}.pdb pdb_files/ {}.m8 tmpFolder --format-output 'query,target,qaln,taln,gapopen,qstart,qend,tstart,tend,qlen,tlen,tcov,qseq,tseq'".format(
        uniprot_input_protein, uniprot_input_protein)
    os.system(s + '\n' + x)


def rin_protein(uniprot_input_protein, list_of_residue):
    c2 = pcmap.contactMap(os.getcwd() + '/' + 'pdb_files/' + uniprot_input_protein + '.pdb')
    df = pd.DataFrame(c2['data'])
    for row in df['root']:
        c = int(row['resID'])
        df['root'] = df['root'].replace([row], c)
    partner1 = []
    interacting_residue = []
    # if you have string and integer
    for row in df['partners']:
        partners = []

        for item in row:
            partners.append(int(item['resID']))
        partner1.append(partners)
        # df.replace(to_replace=row, value=partners)
    df['partner'] = partner1
    del df["partners"]
    matrix = [[0 for i in range(len(df))] for j in range(len(df))]
    columns = [i for i in range(1, len(df) + 1)]
    matrix = pd.DataFrame(matrix)
    matrix.index += 1
    matrix.columns = columns
    for row in df['root']:
        for i in list_of_residue:
            if row == i:
                interacting_residue.append(partner1[row])

    for i in range(len(list_of_residue)):
        for j in interacting_residue[i]:
            matrix.at[list_of_residue[i], j] = 1
    matrix.to_csv(os.getcwd() + '/' + uniprot_input_protein + '.csv', sep='\t')


def loocv_input_protein(primary_list, uniprot_input_protein, result):
    uniprot_list = copy.deepcopy(primary_list[2])
    columns_a = [s for s in uniprot_list if s != uniprot_input_protein]
    # score_matrix = pd.read_pickle("./score_matrix_%s.pkl" % uniprot_input_protein)
    score_matrix = result[3]
    scores_of_residues = score_matrix[columns_a].sum(axis=1)
    score_matrix['scores'] = scores_of_residues
    # probability_matrix = pd.read_pickle("./probability_matrix_%s.pkl" % uniprot_input_protein)
    probability_matrix = result[4]
    loocv_count = []
    for i in columns_a:
        iterative_column = [s for s in columns_a if s != i]
        scores_of_residues = score_matrix[iterative_column].sum(axis=1)
        probability_of_residues = probability_matrix[iterative_column].prod(axis=1)
        loocv_results = pd.DataFrame()
        loocv_results['residues'] = score_matrix[uniprot_input_protein]
        loocv_results['scores'] = scores_of_residues
        loocv_results['probability'] = probability_of_residues
        threshold = 0.01 / (len(iterative_column))
        threshold1 = len(iterative_column) / 2
        result = loocv_results.loc[(loocv_results['probability'] <= threshold) & (loocv_results['scores'] > threshold1)]
        result = result.reset_index(drop=True)
        result = result.iloc[:, [0, 1, -2, -1]]
        b = result['residues'].tolist()
        loocv_count.append(b)
    c = Counter(chain.from_iterable(loocv_count))
    result_counter_loocv = pd.DataFrame.from_dict(c, orient='index').transpose()
    result_counter_loocv.to_csv(r'%s_loocv.csv' % uniprot_input_protein, index=False)
    return loocv_count

# def plot_rin(uniprot_input_protein):
#       from rpy2.robjects.lib.grdevices import jpeg
#       from rpy2.robjects import r
#       draw_circos = r('''draw_circos <- function(matrix, jpeg) {
#       #install and load the package
#       install.packages("circlize", repos = 'https://cran.um.ac.ir/')
#       library(circlize)
#       #read the csv file
#       file_name = paste(getwd(),"/",matrix,".csv", sep ="")
#       data <- read.csv(file_name,sep="")
#       #convert the table to a martix
#       data <- as.matrix(data)
#       #create a chord diagram
#       #chordDiagram(data)#
#       #create a chord diagram but without labeling
#       chordDiagram(data, annotationTrack = "grid", preAllocateTracks = 1)
#
#       #add the labels and axis
#       circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
#         xlim = get.cell.meta.data("xlim")
#         ylim = get.cell.meta.data("ylim")
#         sector.name = get.cell.meta.data("sector.index")
#
#         #print labels
#         circos.text(mean(xlim), ylim[1] + 2.5, sector.name,
#                     facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.6)
#
#         #print axis
#         circos.axis(h = "top", labels.cex = 0.5, major.tick.length = 0.5,
#                     sector.index = sector.name, track.index = 2)
#       }, bg.border = NA)
#       title(main = paste("\n",'The Residue Interaction Network of ',matrix,"\n",'The R indicates the proposed functional mutations the X indicate the interacting residues.'))
#       #saving the plot (high definition)
#       dev.copy(jpeg, paste(getwd(), matrix,'.png',sep=""), width=16, height=16, units="in", res=500)
#       dev.off()
#     }''')
#     return draw_circos(uniprot_input_protein, jpeg)
