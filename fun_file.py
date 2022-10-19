"""
Created on Sat May 21 10:05:42 2022

@author: Hoballa
"""

# Library file
from pandas import *
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
import urllib.request
import re
from Bio import SeqIO
from dash import Dash, html
import dash_bio as dashbio
from Bio import PDB
from sklearn.metrics import *
import alv
import numpy as np
from Bio import AlignIO
import Bio.AlignIO
import rich
from rich_msa import RichAlignment
from tempfile import gettempdir
from urllib.request import urlopen
import json
import matplotlib.pyplot as plt
import ast


def msa_viz(filename):
    app = Dash(__name__)
    data = open(filename).read()
    app.layout = html.Div([
        dashbio.AlignmentChart(
            id='bio-alignmentchart-x-alignment-viewer',
            data=data
        ),
    ])
    if __name__ == '__main__':
        app.run_server(debug=True, use_reloader=False)
    return app.run_server(debug=True, use_reloader=False)


def find_car(s, pat):
    pat = r'(\w*%s\w*)' % pat  # Not thrilled about this line
    return re.findall(pat, s)


# get related family of protein and UniProt human members
def related_family_iphi(p):  ###super family / human/ protein members
    global temp_family, protein_family, protein_list, superfamily, uniprot_input_protein
    super_proteinlist = []
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
    sub_list = []
    if 'superfamily' in protein_family:
        superfamily = protein_family.split('.')[0]
        for l in range(len(x)):
            if superfamily in x[l]:
                sub_list.append(x[l])
        for i in sub_list:
            r1 = x.index(i)
            for j in range(r1 + 1, len(x)):
                if ',' not in x[j]:
                    r2 = j
                    f1 = str(x[r1 + 1:r2 + 1]).split(',')
                    protein_list1 = [word.strip(string.punctuation) for word in f1]
                    super_proteinlist.append(protein_list1)
                    break
    human_protein_list = []
    if len(super_proteinlist) == 0:
        for protein in protein_list:
            # protein = protein.split()
            if '_HUMAN' in protein:
                human_protein_list.append(protein)
        # Index_Input_Protein = [human_protein_list.index(i) for i in human_protein_list if protein_name in i]
        # Index_Input_Protein = [int(i) for i in Index_Input_Protein]
    else:
        for k in super_proteinlist:
            for w in k:
                if '_HUMAN' in w:
                    human_protein_list.append(w)
    uniprot_list = []
    for i in range(len(human_protein_list)):
        uniprot_list.append(human_protein_list[i][human_protein_list[i].find('(') + 1:human_protein_list[i].find(')')])
        if len(find_car(human_protein_list[i], p.upper())) != 0:
            uniprot_input_protein = human_protein_list[i][
                                    human_protein_list[i].find('(') + 1:human_protein_list[i].find(')')]
    gene_name = []
    for i in range(len(human_protein_list)):
        gene_id = find_car(human_protein_list[i], '_')[0]
        gene_id = gene_id.split('_')
        gene_name.append(gene_id[0])
    short_fam_name = "_".join([protein_family][0].split('.')[-1].split(' '))

    return [protein_name.split(), [protein_family], uniprot_list, gene_name, uniprot_input_protein, short_fam_name]


# get PDB file from Alphafold protein structure database

def pdb_structure_alphafold(p):
    w = os.getcwd() + '/pdb_files/{}.pdb'.format(p)
    url = "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v3.pdb" % p
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


def hotspot_profiler(uniprot_set, genename_list):
    hotspot_profile1 = []
    dhotspots = pd.read_csv("3dhotspots.csv")
    cancerhotspots = pd.read_csv("cancerhotspots.csv")
    dhotspots_list = []
    cancerhotspots_list = []
    for gene in uniprot_set:
        value = dhotspots[dhotspots['Gene'].str.contains(gene)]
        s = value['Residue']
        dhotspots_list.append(s.tolist())
        value = cancerhotspots[cancerhotspots['Gene'].str.contains(gene)]
        t = value['Residue']
        cancerhotspots_list.append(t.tolist())
    for i in range(len(uniprot_set)):
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
    profile_last = {'uniprot_ID': uniprot_set, 'genes': genename_list, 'hotspot_profile': hotspot_profile1}
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


def foldseek_alig_matrix(alignment_file):  # , uniprot_input_protein):
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
    # percentage_of_alignment['seq_sim'] = seq_sim
    # percentage_of_alignment.to_csv(os.getcwd() + '/' + f'{uniprot_input_protein} percentage of alignment' + '.csv', sep='\t')
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
    return alignment_matrix_seq, alignment_matrix_res  # , percentage_of_alignment


def result_matrix(uniprot_input_protein, alinment_matrix, biomuta_profile, hotspot_profile, primary_list):
    residue_matrix = alinment_matrix[1]
    residue_matrix = residue_matrix.fillna(0)
    residue_matrix = residue_matrix.astype(int)
    mutation_probabilty = {}
    hotspot_probabilty = {}
    # nonmutation_probability = {}
    for column in residue_matrix:
        biomuta_mut = biomuta_profile[biomuta_profile['uniprot_ID'].isin([column])]['biomuta_profile'].tolist()[0]
        hotspot_mut = hotspot_profile[hotspot_profile['uniprot_ID'].isin([column])]['hotspot_profile'].tolist()[0]
        mutation_probabilty[column] = ((len(biomuta_mut) - len(hotspot_mut)) / alinment_matrix[0].loc[
            'protein_length', column])  #### nonhotspotmutation prob = #nonhot mutation / len protein
        hotspot_probabilty[column] = (len(hotspot_mut) / alinment_matrix[0].loc[
            'protein_length', column])  #### hotspot mutation probability = #hotspot/len protein
        ###non_mutation = 1-(prob_hotspot + prob_nonhotspot)
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
    uniprot_list = list(residue_matrix.columns)  # copy.deepcopy(primary_list[2])
    columns = [s for s in uniprot_list if s != uniprot_input_protein]
    scores_of_residues = residue_matrix[columns].sum(axis=1)
    max_score = max(scores_of_residues)
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
    # maximum = max(probability_of_residues)
    # minimum = min(probability_of_residues)
    # new_p = maximum-minimum
    alinment_matrix[0]['probability'] = probability_of_residues
    residue_matrix['probability'] = probability_of_residues
    matrix = alinment_matrix[0]
    probability_matrix = copy.deepcopy(alinment_matrix[1])
    matrix.insert(loc=0, column='residue_number', value=range(1, len(matrix) + 1))
    # matrix.drop(index=matrix.index[-1], axis=0, inplace=True)
    # residue_matrix.to_pickle("./probability_loocv_%s.pkl" % uniprot_input_protein)
    threshold = 0.01 / (len(primary_list[2]) - 1)
    threshold1 = len(primary_list[2]) / 2  # max_score/2
    # percentage = (threshold*100)/maximum
    result = matrix.loc[
        (matrix['probability'] < threshold) & (matrix['scores'] > threshold1)]
    result = result.reset_index(drop=True)
    result = result.iloc[:, [0, 1, -2, -1]]
    b = result['residue_number'].tolist()
    pd.set_option('display.float_format', '{:.3g}'.format)
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
    x = "foldseek easy-search pdb_files/{}.pdb pdb_files/ m8/{}.m8 tmpFolder --format-output 'query,target,qaln,taln," \
        "gapopen,qstart,qend,tstart,tend,qlen,tlen,tcov,qseq,tseq'".format(
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


"""
def plot_rin(uniprot_input_protein):
      from rpy2.robjects.lib.grdevices import jpeg
      from rpy2.robjects import r
      draw_circos = r('''draw_circos <- function(matrix, jpeg) {
      #install and load the package
      install.packages("circlize", repos = 'https://cran.um.ac.ir/')
      library(circlize)
      #read the csv file
      file_name = paste(getwd(),"/",matrix,".csv", sep ="")
      data <- read.csv(file_name,sep="")
      #convert the table to a martix
      data <- as.matrix(data)
      #create a chord diagram
      #chordDiagram(data)#
      #create a chord diagram but without labeling
      chordDiagram(data, annotationTrack = "grid", preAllocateTracks = 1)

      #add the labels and axis
      circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")

        #print labels
        circos.text(mean(xlim), ylim[1] + 2.5, sector.name,
                    facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.6)

        #print axis circos.axis(h = "top", labels.cex = 0.5, major.tick.length = 0.5, sector.index = sector.name, 
        track.index = 2) }, bg.border = NA) title(main = paste("\n",'The Residue Interaction Network of ',matrix,
        "\n",'The R indicates the proposed functional mutations the X indicate the interacting residues.')) #saving 
        the plot (high definition) dev.copy(jpeg, paste(getwd(), matrix,'.png',sep=""), width=16, height=16, 
        units="in", res=500) dev.off() }''') return draw_circos(uniprot_input_protein, jpeg) """


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
        if len(find_car(human_protein_list[i], p.upper())) != 0:
            uniprot_input_protein = human_protein_list[i][
                                    human_protein_list[i].find('(') + 1:human_protein_list[i].find(')')]
    gene_name = []
    for i in range(len(human_protein_list)):
        gene_id = find_car(human_protein_list[i], '_')[0]
        gene_id = gene_id.split('_')
        gene_name.append(gene_id[0])
    short_fam_name = "_".join([protein_family][0].split('.')[-1].split(' '))

    return [protein_name.split(), [protein_family], uniprot_list, gene_name, uniprot_input_protein, short_fam_name]


def write_fasta(dictionary, filename):
    """
    Takes a dictionary and writes it to a fasta file
    Must specify the filename when caling the function
    """
    import textwrap
    with open(filename, "a") as outfile:
        for key, value in dictionary.items():
            outfile.write('>' + key + "\n")
            outfile.write("\n".join(textwrap.wrap(value, 60)))
            outfile.write("\n")
    return outfile
    # print("Success! File written")


def clustlo(fasta_file, output_filename):
    s = "./clustalo -i {in_seqs} -o {out_seqs} -v ".format(in_seqs=fasta_file, out_seqs=output_filename)
    os.system(s)


# iterate over family and get equlibrium of predictions
def iter_till_equi(primary_list, biomuta_profile, hotspot_profile):
    iter_hot_results = [primary_list[2], hotspot_profile['hotspot_profile'].tolist()]
    hotspot_profile1 = copy.deepcopy(hotspot_profile)
    equi_result = []
    equi_prog = []
    iter_num = 0
    while len(equi_result) < len(primary_list[2]):
        temp_result = []
        equi_result = []
        for i in primary_list[2]:
            alignment_matrix = foldseek_alig_matrix('./m8/%s.m8' % i)
            result = result_matrix(i, alignment_matrix, biomuta_profile, hotspot_profile1, primary_list)
            temp_result.append(result[2])
        pre_result = hotspot_profile1['hotspot_profile'].tolist()
        for k in range(len(primary_list[2])):
            new_set = list(set(temp_result[k]) | set(pre_result[k]))
            hotspot_profile1.at[k, 'hotspot_profile'] = new_set
        new_result = hotspot_profile1['hotspot_profile'].tolist()
        for j in range(len(primary_list[2])):
            if set(new_result[j]) == set(pre_result[j]):
                equi_result.append(1)
        equi_prog.append(sum(equi_result))
        iter_num += 1

    number_of_iterations = iter_num
    iter_hot_results.append(temp_result)
    iter_hot_results = pd.DataFrame(iter_hot_results).transpose()
    # iter_hot_results.to_pickle("./pkl/iter_hot_results%s.pkl" % "_".join(primary_list[1][0].split('.')[-1].split('
    # '))) iter_hot_results = pd.read_pickle("./iter_hot_results_%s.pkl" % primary_list[1])
    iter_hot_results = iter_hot_results  # .transpose()
    # print('there were : ', number_of_iterations, 'iteration')
    return number_of_iterations, iter_hot_results


# one iteration run
def iter_one_time(primary_list, biomuta_profile, hotspot_profile):
    iter_hot_results = [primary_list[2], hotspot_profile['hotspot_profile'].tolist()]
    hotspot_profile1 = copy.deepcopy(hotspot_profile)
    temp_result = []
    # equi_result = []
    # equi_prog = []
    for i in primary_list[2]:
        alignment_matrix = foldseek_alig_matrix('./m8/%s.m8' % i)
        result = result_matrix(i, alignment_matrix, biomuta_profile, hotspot_profile1, primary_list)
        temp_result.append(result[2])
    pre_result = hotspot_profile1['hotspot_profile'].tolist()
    for k in range(len(primary_list[2])):
        new_set = list(set(temp_result[k]) | set(pre_result[k]))
        hotspot_profile1.at[k, 'hotspot_profile'] = sorted(new_set)
    new_result = hotspot_profile1['hotspot_profile'].tolist()
    iter_hot_results.append(new_result)
    iter_hot_results = pd.DataFrame(iter_hot_results).transpose()
    # iter_hot_results.to_pickle("./pkl/iter_hot_results%s.pkl" % "_".join(primary_list[1][0].split('.')[-1].split('
    # '))) iter_hot_results = pd.read_pickle("./iter_hot_results_%s.pkl" % primary_list[1]) iter_hot_results =
    # iter_hot_results  # .transpose()
    return iter_hot_results


# get protein seq from alignment
def get_seq(primary_list):
    qseq_list = []
    for i in primary_list[2]:
        s = pd.read_csv('./m8/%s.m8' % i, sep='\t', header=None)
        s.columns = ['query', 'target', 'qaln', 'taln', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'qlen', 'tlen',
                     'tcov', 'qseq', 'tseq']
        seq_dic = {i: s['qseq'][0]}
        qseq_list.append(seq_dic)
    return qseq_list


def i_phi(primary_list, biomuta_profile, hotspot_profile, primary_list2):
    threshold = len(primary_list2) / 2
    global index_in_msa, msa_aa_list, msa_position_counters, number_of_predictions, scores_of_residues, columns
    short_fam_name = "_".join(primary_list[1][0].split('.')[-1].split(' '))
    number_of_predictions = ''
    hotspot_profile_i = copy.deepcopy(hotspot_profile)
    equi_result = []
    iter_num = 0
    equi_prog = []
    old_result = []
    ###
    qseq_list = get_seq(primary_list)
    if os.path.exists('./fa/fasta_seq%s.fa' % short_fam_name):
        os.remove('./fa/fasta_seq%s.fa' % short_fam_name)
    if os.path.exists('./msa/msa_file%s.fa' % short_fam_name):
        os.remove('./msa/msa_file%s.fa' % short_fam_name)
    for j in qseq_list:
        write_fasta(j, './fa/fasta_seq%s.fa' % short_fam_name)
    clustlo('./fa/fasta_seq%s.fa' % short_fam_name, './msa/msa_file%s.fa' % short_fam_name)
    msa_list = []
    fasta_sequences = SeqIO.parse(open('./msa/msa_file%s.fa' % short_fam_name), 'fasta')
    uniprot_list = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        msa_list.append({name: sequence})
        uniprot_list.append(name)
    msa_aa_list = []
    msa_position_counters = []
    for j in msa_list:
        for s in j.values():
            t = [*s]
            msa_position_counters.append(String_character_counter(s, 1))
            msa_aa_list.append(t)
    msa_aa_list = pd.DataFrame(msa_aa_list)
    msa_position_counters = pd.DataFrame(msa_position_counters)
    msa_aa_list.index = uniprot_list
    msa_position_counters.index = uniprot_list
    msa_aa_list = msa_aa_list.transpose()
    msa_position_counters = msa_position_counters.transpose()
    msa_position_counters.to_pickle("./pkl/msa_position_%s.pkl" % short_fam_name)
    ###
    while len(equi_result) < len(primary_list2):
        iter_hot_results = iter_one_time(primary_list, biomuta_profile, hotspot_profile_i)
        ###
        msa_residue_list = []
        for i in msa_list:
            for v in i.values():
                x = String_character_counter(v, 1)
                msa_residue_list.append(x)
        msa_residue_list = pd.DataFrame(msa_residue_list)
        msa_residue_list.index = uniprot_list
        msa_residue_list = msa_residue_list.transpose()
        for column in msa_residue_list:
            iter_hotspots = iter_hot_results[iter_hot_results[0].isin([column])][2].tolist()[0]
            for row in msa_residue_list[column]:
                if row in iter_hotspots:
                    msa_residue_list[column] = msa_residue_list[column].replace([row], 1)
                else:
                    msa_residue_list[column] = msa_residue_list[column].replace([row], 0)
        columns = [s for s in uniprot_list]
        scores_of_residues = msa_residue_list[columns].sum(axis=1)
        msa_aa_list['scores'] = scores_of_residues
        msa_aa_list.to_pickle("./pkl/msa_aa_list_%s.pkl" % short_fam_name)
        scores_adjusted = []
        for i in scores_of_residues:
            if i > len(columns) / 2:
                scores_adjusted.append(1)
            else:
                scores_adjusted.append(0)
        dic_to_alig = {'score': ''.join(map(str, scores_adjusted))}
        write_fasta(dic_to_alig, './msa/msa_file%s.fa' % short_fam_name)
        index_in_msa = [i for i, j in enumerate(scores_adjusted) if j == 1]
        if index_in_msa == old_result:
            break
        else:
            old_result = copy.deepcopy(index_in_msa)
        new_hot_from_msa = []
        for column in msa_position_counters:
            col_hot = []
            for j in index_in_msa:
                q = msa_position_counters[column][j]
                if q != 0:
                    col_hot.append(q)
            new_hot_from_msa.append(col_hot)
        pre_result = hotspot_profile_i['hotspot_profile'].tolist()
        new_result = new_hot_from_msa  # hotspot_profile_i['hotspot_profile'].tolist()
        equi_result = []
        for j in range(len(columns)):
            if set(new_result[j]) == set(pre_result[j]):
                equi_result.append(1)
            else:
                hotspot_profile_i.at[j, 'hotspot_profile'] = sorted(list(set(new_result[j]) | set(pre_result[j])))
        equi_prog.append(sum(equi_result))
        number_of_predictions = sum(scores_adjusted)
        iter_num += 1
        # print('this iteration in the number %s iteration' % iter_num)
    scores_adjusted = []
    for i in scores_of_residues:
        if i > threshold:
            scores_adjusted.append(1)
        else:
            scores_adjusted.append(0)
    dic_to_alig = {'score': ''.join(map(str, scores_adjusted))}
    write_fasta(dic_to_alig, './msa/msa_file%s.fa' % short_fam_name)
    index_in_msa = [i for i, j in enumerate(scores_adjusted) if j == 1]
    new_hot_from_msa = []
    for column in msa_position_counters:
        col_hot = []
        for j in index_in_msa:
            q = msa_position_counters[column][j]
            if q != 0:
                col_hot.append(q)
        new_hot_from_msa.append(col_hot)
    for k in range(len(columns)):
        new_set = list(sorted(set(new_hot_from_msa[k])))
        hotspot_profile_i.at[k, 'hotspot_profile'] = new_set
    return index_in_msa, iter_num, msa_aa_list, hotspot_profile_i, msa_position_counters, number_of_predictions


def analyzie_mutaligner(url):  # the json url obtained from the site of mutationaligner
    # url = "http://mutationaligner.org/api/domains/Ras"
    response = urlopen(url)
    data_json = json.loads(response.read())
    df_positions = pd.DataFrame(data_json['positions'])
    # df_genes = pd.DataFrame(data_json['genes'])
    df_positions1 = df_positions['mutations']
    listall = []
    for i in df_positions1:
        for j in i:
            listall.append(j)
    df_positions3 = pd.DataFrame(listall)
    for j in range(len(df_positions3['position'])):
        df_positions3['position'][j] = int(re.findall('[0-9]+', df_positions3['position'][j])[0])
    unique_gene_list = df_positions3.geneName.unique()
    muta_aligner_results = []
    for i in unique_gene_list:
        x = sorted(list(set(df_positions3.loc[df_positions3['geneName'] == i, 'position'].tolist())))
        dict_spec = {i: x}
        muta_aligner_results.append(dict_spec)
    return muta_aligner_results


def edrn_muta_profile(primary_list):
    uniprot_set = primary_list[2]
    gene_name = primary_list[3]
    req_cols = ['UniProtAC', 'Gene_symbol', 'Position(A)']
    df = pd.read_csv("biomuta-edrn.csv", usecols=req_cols)
    temp = []
    for j in uniprot_set:
        contain_values = df[df['UniProtAC'].str.contains(j)]
        s = contain_values['Position(A)']
        w = sorted(s.unique().tolist())
        temp.append(w)
    mutation_profile = pd.DataFrame(list(zip(uniprot_set, gene_name, temp)),
                                    columns=['uniprot_ID', 'genes', 'hotspot_profile'])
    return mutation_profile


def similar_iphi_database():
    my_file = open("similar.txt", 'r')
    my_file.read()
    list_of_families = []
    index_of_families = []
    with open("similar.txt", 'r') as file:
        x = file.readlines()
    for i in range(len(x)):
        if ',' not in x[i] and 'family' in x[i]:
            temp_family = x[i]
            list_of_families.append(temp_family.split('\n')[0])
            index_of_families.append(i)
    all_protein_list = []
    for i in range(len(index_of_families)):
        try:
            w = index_of_families[i + 1]
            human_protein_list = []
            f = str(x[index_of_families[i] + 1:w]).split(',')
            protein_list = [word.strip(string.punctuation) for word in f]
            for protein in protein_list:
                # protein = protein.split()
                if '_HUMAN' in protein:
                    human_protein_list.append(protein)
            all_protein_list.append(human_protein_list)
        except:
            w = len(x)
            human_protein_list = []
            f = str(x[index_of_families[i] + 1:w]).split(',')
            protein_list = [word.strip(string.punctuation) for word in f]
            for protein in protein_list:
                # protein = protein.split()
                if '_HUMAN' in protein:
                    human_protein_list.append(protein)
            all_protein_list.append(human_protein_list)
    uniprot_list = []
    gene_name_list = []
    for i in all_protein_list:
        temp_uniprot = []
        temp_genename = []
        if len(i) != 0:
            for j in i:
                temp_uniprot.append(j[j.find('(') + 1:j.find(')')])
                gene_id = find_car(j, '_')[0]
                gene_id = gene_id.split('_')
                temp_genename.append(gene_id[0])
            uniprot_list.append(temp_uniprot)
            gene_name_list.append(temp_genename)
        else:
            uniprot_list.append(None)
            gene_name_list.append(None)
    results = [list_of_families, gene_name_list, uniprot_list]
    columns = ['families', 'gene name', 'uniprot_ACCID']
    family_Dataframe = pd.DataFrame(results).transpose()
    family_Dataframe.columns = columns
    family_Dataframe.to_csv("./csv/family_Dataframe.csv", index=True)
    family_Dataframe.to_pickle("./pkl/family_Dataframe.pkl")


# try:
#     biomuta_profile = pd.read_pickle("./pkl/input_biomuta%s.pkl" % short_fam_name)
#     hotspot_profile = pd.read_pickle("./pkl/input_hotspot%s.pkl" % short_fam_name)
# except:
#     biomuta_profile = biomuta_resorce(primary_list[2], primary_list[3])
#     biomuta_profile.to_pickle("./pkl/input_biomuta%s.pkl" % short_fam_name)
#     hotspot_profile = hotspot_profiler(primary_list[2], primary_list[3])
#     hotspot_profile.to_pickle("./pkl/input_hotspot%s.pkl" % short_fam_name)
# biomuta_profile = biomuta_resorce(primary_list[2], primary_list[3])
# biomuta_profile.to_pickle("./pkl/input_biomuta%s.pkl" % short_fam_name)
# biomuta_profile = pd.read_pickle("./pkl/input_biomuta%s.pkl" % short_fam_name)
# hotspot_profile = hotspot_profiler(primary_list[2], primary_list[3])
# hotspot_profile.to_pickle("./pkl/input_hotspot%s.pkl" % short_fam_name)
# hotspot_profile = pd.read_pickle("./pkl/input_hotspot%s.pkl" % short_fam_name)
'''
hotspot_profile_edrn = edrn_muta_profile(primary_list)
hotspot_profile_edrn.to_pickle("./pkl/edrn_input_hotspot%s.pkl" % short_fam_name)
hotspot_profile_edrn = pd.read_pickle("./pkl/edrn_input_hotspot%s.pkl" % short_fam_name)
'''

# view file in HTML (holds debug in code)
# msa_viz('./msa/msa_file%s.fa' % short_fam_name)
'''
# view file in consol
msa = AlignIO.read('./msa/msa_file%s.fa' % short_fam_name, 'fasta')
alv.view(msa)
# print inside box in consol
msa = Bio.AlignIO.read('./msa/msa_file%s.fa' % short_fam_name, "fasta")
viewer = RichAlignment(
    names=[record.id for record in msa],
    sequences=[str(record.seq) for record in msa],
)
panel = rich.panel.Panel(viewer, title='msa_file%s.fa' % short_fam_name)
rich.print(panel)
'''


# input_protein = 'p53' #XCR1
# primary_list = related_family(input_protein)
# short_fam_name = primary_list[5] #G-protein coupled receptor 1 family
# uniprot_input_protein = primary_list[4]

def IPHI_system():
    database = pd.read_csv("./final_database.csv")
    search_option = input('what do you want to search: "uniprot/family/multiple families"')
    if search_option == 'uniprot':
        gene_input = str(input('please inter gene name to find family :'))
        database["Indexes"] = database['uniprot_ACCID'].str.find(gene_input)
        w = database['Indexes'][database['Indexes'] == 2].index[0]
        contain_values = database.iloc[w]
        protein_family = contain_values['families']
        uniprot_list = ast.literal_eval(contain_values['uniprot_ACCID'])
        gene_name = ast.literal_eval(contain_values['gene name'])
        uniprot_input_protein = gene_input
        short_fam_name = "_".join([protein_family][0].split('.')[-1].split(' '))  # protein_family
        primary_list = [gene_input, protein_family, uniprot_list, gene_name, uniprot_input_protein, short_fam_name]
        biomuta_profile = ast.literal_eval(contain_values['biomuta'])
        hotspot_profile = ast.literal_eval(contain_values['hotspot'])
    if search_option == 'family':
        family_input = str(input('please provide the family name to find members : '))
        # df1[df1['ids'] == "aball"]
        database["Indexes"] = database['families'].str.find(family_input)
        s = database['Indexes']
        w = database['Indexes'][database['Indexes'] == 0].index[0]
        contain_values = database.iloc[w]
        protein_family = contain_values['families']
        uniprot_list = ast.literal_eval(contain_values['uniprot_ACCID'])
        gene_name = ast.literal_eval(contain_values['gene name'])
        uniprot_input_protein = ''
        short_fam_name = "_".join([protein_family][0].split('.')[-1].split(' '))  # protein_family
        primary_list = ['', protein_family, uniprot_list, gene_name, '', short_fam_name]
        biomuta_profile = ast.literal_eval(contain_values['biomuta'])
        hotspot_profile = ast.literal_eval(contain_values['hotspot'])
    if search_option == 'multiple families':
        multiple_input = input('please provide a list of families you want to search :')
        multiple_input = multiple_input.split(',')
        # multiple_input = ast.literal_eval([multiple_input])
        family_names = multiple_input
        uniprot_list_all = []
        gene_list_all = []
        # short_family_names = []
        biomuta_all = []
        hotspot_all = []
        for i in multiple_input:
            database["Indexes"] = database['families'].str.find(i)
            s = database['Indexes']
            w = database['Indexes'][database['Indexes'] == 0].index[0]
            contain_values = database.iloc[w]
            protein_family = contain_values['families']
            uniprot_list = ast.literal_eval(contain_values['uniprot_ACCID'])
            gene_name = ast.literal_eval(contain_values['gene name'])
            uniprot_input_protein = ''
            biomuta_profile = ast.literal_eval(contain_values['biomuta'])
            hotspot_profile = ast.literal_eval(contain_values['hotspot'])
            uniprot_list_all.append(uniprot_list)
            gene_list_all.append(gene_name)
            biomuta_all.append(biomuta_profile)
            hotspot_all.append(hotspot_profile)
        uniprot_list_all = [item for sublist in uniprot_list_all for item in sublist]
        gene_list_all = [item for sublist in gene_list_all for item in sublist]
        biomuta_profile = [item for sublist in biomuta_all for item in sublist]
        hotspot_profile = [item for sublist in hotspot_all for item in sublist]
        primary_list = ['', family_names, uniprot_list_all, gene_list_all, '', family_names]

    clear_folder(os.getcwd() + '/pdb_files')
    removed_uni_structures = []
    removed_genes = []
    for i in primary_list[2]:
        try:
            pdb_structure_alphafold(i)
        except:
            j = primary_list[2].index(i)
            removed_uni_structures.append(i)
            removed_genes.append(primary_list[3][j])
    primary_list[2] = [x for x in primary_list[2] if
                       x not in removed_uni_structures]  # primary_list[2] - removed_uni_structures
    primary_list[3] = [x for x in primary_list[3] if x not in removed_genes]  # primary_list[3] - removed_genes

    for i in primary_list[2]:
        if os.path.exists('./m8/%s.m8' % i):
            continue
        else:
            foldseek(i)

    iterative_system = i_phi(primary_list, biomuta_profile, hotspot_profile, primary_list[2])
    msa_position = iterative_system[4]
    msa_aa_list = iterative_system[2]
    results_in_proteins = iterative_system[3]
    consensue_results = iterative_system[0]
    results_in_proteins.to_csv("./csv/results_in_proteins.csv", index=False)

    set1, set2, set3, unmutated_results, new_hotspots, result = '', '', '', '', '', ''
    results = []
    for i in range(len(primary_list[3])):
        set1 = results_in_proteins.at[i, 'hotspot_profile']
        set2 = biomuta_profile.at[i, 'biomuta_profile']
        set3 = hotspot_profile.at[i, 'hotspot_profile']
        unmutated_results = [x for x in set1 if x not in set2]
        mutated_results = [x for x in set1 if x in set2]
        new_hotspots = [x for x in set1 if x in mutated_results]
        iphi_new_hotspots = [x for x in mutated_results if x not in set3]
        result = [primary_list[2][i], primary_list[3][i], set2, set3, set1, new_hotspots, iphi_new_hotspots,
                  unmutated_results]
        results.append(result)
    columns = ['UniProt', 'GeneName', 'Biomuta', 'HotSpots', 'I-phi', 'all Hotspots', 'New Htspots',
               'Unmutated results']
    iphi_report = pd.DataFrame(results, columns=columns)
    iphi_report.to_csv("./csv/iphi_report.csv", index=False)

    return iphi_report


def flat(lis):
    flatList = []
    # Iterate with outer list
    for element in lis:
        if type(element) is list:
            # Check if type is list than iterate through the sublist
            for item in element:
                flatList.append(item)
        else:
            flatList.append(element)
    return flatList


def get_pdb_chain(uni, pdb_id, chain_id):
    # pdb_id, chain_id = '1ATP', 'I'
    work_dir = './pdb_files'
    PDB.PDBList(verbose=False).retrieve_pdb_file(pdb_id, pdir=work_dir, file_format='pdb')
    # os.rename('{0}/pdb{1}.ent'.format(work_dir, pdb_id.lower()), '{0}/{1}.pdb'.format(work_dir, uni))
    biopdb_name = '{0}/pdb{1}.ent'.format(work_dir, pdb_id.lower())
    # Read the PDB file and extract the chain from structure[0]
    model = PDB.PDBParser(PERMISSIVE=1, QUIET=1).get_structure(pdb_id, biopdb_name)[0]
    io = PDB.PDBIO()
    io.set_structure(model[chain_id])
    io.save('{0}/{1}.pdb'.format(work_dir, uni))
    test = os.listdir(work_dir)
    for item in test:
        if item.endswith(".ent"):
            os.remove(os.path.join(work_dir, item))


def iphi_rcsb_alfa(uni_list, link_database):
    database = pd.read_csv(link_database)
    reu_compelet = []
    reu_compelet1 = []
    for i in uni_list:
        # for j in database.iloc():
        wanted = database[(database
                           .apply(lambda r: r.astype('string').str.contains(i)
                                  .any(), axis=1))]
        primary_list = [i, wanted['families'].to_string(),
                        ast.literal_eval(wanted['uniprot_ACCID'].iloc[:].tolist()[0]),
                        ast.literal_eval(wanted['gene name'].iloc[:].tolist()[0]), i, wanted['families'].to_string()]
        biomuta_profile = pd.DataFrame(
            {'uniprot_ID': ast.literal_eval(wanted['uniprot_ACCID'].iloc[:].tolist()[0]),
             'genes': ast.literal_eval(wanted['gene name'].iloc[:].tolist()[0]),
             'biomuta_profile': ast.literal_eval(wanted['biomuta'].iloc[:].tolist()[0])})
        hotspot_profile = pd.DataFrame(
            {'uniprot_ID': ast.literal_eval(wanted['uniprot_ACCID'].iloc[:].tolist()[0]),
             'genes': ast.literal_eval(wanted['gene name'].iloc[:].tolist()[0]),
             'hotspot_profile': ast.literal_eval(wanted['hotspot'].iloc[:].tolist()[0])})
        try:
            #############################

            # # download structure of alphafold

            clear_folder(os.getcwd() + '/pdb_files')
            removed_uni_structures = []
            removed_genes = []
            for j in primary_list[2]:
                try:
                    pdb_structure_alphafold(j)
                except:
                    k = primary_list[2].index(j)
                    removed_uni_structures.append(k)
                    removed_genes.append(primary_list[3][k])
            primary_list[2] = [x for x in primary_list[2] if
                               x not in removed_uni_structures]  # primary_list[2] - removed_uni_structures
            primary_list[3] = [x for x in primary_list[3] if x not in removed_genes]  # primary_list[3] - removed_genes

            try:
                foldseek(primary_list[4])

                alignment_matrix = foldseek_alig_matrix('./m8/%s.m8' % primary_list[4])
                result = result_matrix(primary_list[4], alignment_matrix, biomuta_profile, hotspot_profile,
                                       primary_list)
                wantedhot = hotspot_profile[(hotspot_profile
                                             .apply(lambda r: r.astype('string').str.contains(i)
                                                    .any(), axis=1))]
                wantedmuta = biomuta_profile[(biomuta_profile
                                              .apply(lambda r: r.astype('string').str.contains(i)
                                                     .any(), axis=1))]
                reu_compelet.append({'input_protein': primary_list[0],
                                     'o_biomuta': wantedmuta['biomuta_profile'].tolist()[0],
                                     'o_hotspot': wantedhot['hotspot_profile'].tolist()[0],
                                     'p_positions': result[2]})

            except:
                wantedhot = hotspot_profile[(hotspot_profile
                                             .apply(lambda r: r.astype('string').str.contains(i)
                                                    .any(), axis=1))]
                wantedmuta = biomuta_profile[(biomuta_profile
                                              .apply(lambda r: r.astype('string').str.contains(i)
                                                     .any(), axis=1))]
                reu_compelet.append({'input_protein': primary_list[0],
                                     'o_biomuta': wantedmuta['biomuta_profile'].tolist()[0],
                                     'o_hotspot': wantedhot['hotspot_profile'].tolist()[0],
                                     'p_positions': None})
        except:
            wantedhot = hotspot_profile[(hotspot_profile
                                         .apply(lambda r: r.astype('string').str.contains(i)
                                                .any(), axis=1))]
            wantedmuta = biomuta_profile[(biomuta_profile
                                          .apply(lambda r: r.astype('string').str.contains(i)
                                                 .any(), axis=1))]
            reu_compelet.append({'input_protein': primary_list[0],
                                 'o_biomuta': wantedmuta['biomuta_profile'].tolist()[0],
                                 'o_hotspot': wantedhot['hotspot_profile'].tolist()[0],
                                 'p_positions': None})
        reu_compelet = pd.DataFrame(reu_compelet)
        reu_compelet.to_pickle('./reu_complete_ALFA.pkl')
        try:
            ############################

            # download from RCSB

            clear_folder(os.getcwd() + '/pdb_files')
            look_up = pd.read_pickle('./look_up.pkl')
            rcsb_found = []
            for j in primary_list[2]:
                for k in look_up.iloc():
                    if j in k['UNI']:
                        rcsb_found.append({'UNI': j, 'PDB': k['PDB'], 'CHAIN': k['CHAIN']})
            rcsb_found = pd.DataFrame(rcsb_found)
            for e in rcsb_found.iloc():
                get_pdb_chain(e['UNI'], e['PDB'], e['CHAIN'])
            primary_list[2] = rcsb_found['UNI'].tolist()  # primary_list[2] - removed
            primary_list[3] = rcsb_found['UNI'].tolist()  # primary_list[3] - removed_genes

            ########################
            try:
                foldseek(primary_list[4])

                alignment_matrix = foldseek_alig_matrix('./m8/%s.m8' % primary_list[4])
                result = result_matrix(primary_list[4], alignment_matrix, biomuta_profile, hotspot_profile,
                                       primary_list)
                wantedhot = hotspot_profile[(hotspot_profile
                                             .apply(lambda r: r.astype('string').str.contains(i)
                                                    .any(), axis=1))]
                wantedmuta = biomuta_profile[(biomuta_profile
                                              .apply(lambda r: r.astype('string').str.contains(i)
                                                     .any(), axis=1))]
                reu_compelet.append({'input_protein': primary_list[0],
                                     'o_biomuta': wantedmuta['biomuta_profile'].tolist()[0],
                                     'o_hotspot': wantedhot['hotspot_profile'].tolist()[0],
                                     'p_positions': result[2]})

            except:
                wantedhot = hotspot_profile[(hotspot_profile
                                             .apply(lambda r: r.astype('string').str.contains(i)
                                                    .any(), axis=1))]
                wantedmuta = biomuta_profile[(biomuta_profile
                                              .apply(lambda r: r.astype('string').str.contains(i)
                                                     .any(), axis=1))]
                reu_compelet1.append({'input_protein': primary_list[0],
                                      'o_biomuta': wantedmuta['biomuta_profile'].tolist()[0],
                                      'o_hotspot': wantedhot['hotspot_profile'].tolist()[0],
                                      'p_positions': None})
        except:
            wantedhot = hotspot_profile[(hotspot_profile
                                         .apply(lambda r: r.astype('string').str.contains(i)
                                                .any(), axis=1))]
            wantedmuta = biomuta_profile[(biomuta_profile
                                          .apply(lambda r: r.astype('string').str.contains(i)
                                                 .any(), axis=1))]
            reu_compelet1.append({'input_protein': primary_list[0],
                                  'o_biomuta': wantedmuta['biomuta_profile'].tolist()[0],
                                  'o_hotspot': wantedhot['hotspot_profile'].tolist()[0],
                                  'p_positions': None})
        reu_compelet1 = pd.DataFrame(reu_compelet1)
        reu_compelet1.to_pickle('./reu_complete_ALFA.pkl')

    return reu_compelet, reu_compelet1


def protein_length(uniprot_id):
    data = urlopen("https://rest.uniprot.org/uniprotkb/" + uniprot_id + ".txt").read()
    data1 = data.decode()
    data = data1.split('\n')
    len_Protein = [int(s) for s in data[0].split() if s.isdigit()]
    len_Protein1 = len_Protein[0]
    return len_Protein1


def iphi_rcsb_alfa_mix(uni_list, link_database):
    database = pd.read_csv(link_database)
    reu_compelet = []
    for i in uni_list:
        # for j in database.iloc():
        wanted = database[(database
                           .apply(lambda r: r.astype('string').str.contains(i)
                                  .any(), axis=1))]
        primary_list = [i, wanted['families'].to_string(),
                        ast.literal_eval(wanted['uniprot_ACCID'].iloc[:].tolist()[0]),
                        ast.literal_eval(wanted['gene name'].iloc[:].tolist()[0]), i, wanted['families'].to_string()]
        biomuta_profile = pd.DataFrame(
            {'uniprot_ID': ast.literal_eval(wanted['uniprot_ACCID'].iloc[:].tolist()[0]),
             'genes': ast.literal_eval(wanted['gene name'].iloc[:].tolist()[0]),
             'biomuta_profile': ast.literal_eval(wanted['biomuta'].iloc[:].tolist()[0])})
        hotspot_profile = pd.DataFrame(
            {'uniprot_ID': ast.literal_eval(wanted['uniprot_ACCID'].iloc[:].tolist()[0]),
             'genes': ast.literal_eval(wanted['gene name'].iloc[:].tolist()[0]),
             'hotspot_profile': ast.literal_eval(wanted['hotspot'].iloc[:].tolist()[0])})

            ############################

            # mixed if no experimental then predicted
        try:
            clear_folder(os.getcwd() + '/pdb_files')
            look_up = pd.read_pickle('./look_up.pkl')
            rcsb_found = []
            alpha_found = []
            listo = look_up['UNI'].tolist()
            for j in primary_list[2]:
                if j in listo:
                    s = look_up.loc[look_up['UNI'] == j]
                    rcsb_found.append(s)
                else:
                    alpha_found.append(j)
            rcsb_found = pd.concat(rcsb_found, ignore_index=True)
            for e in rcsb_found.iloc():
                get_pdb_chain(e['UNI'], e['PDB'], e['CHAIN'])
            removed_uni_structures = []
            removed_genes = []
            for d in alpha_found:
                try:
                    pdb_structure_alphafold(d)
                except:
                    k = primary_list[2].index(d)
                    removed_uni_structures.append(k)
                    removed_genes.append(primary_list[3][k])
            primary_list[2] = [x for x in primary_list[2] if
                               x not in removed_uni_structures]  # primary_list[2] - removed_uni_structures
            primary_list[3] = [x for x in primary_list[3] if x not in removed_genes]

            ########################

            foldseek(primary_list[4])

            alignment_matrix = foldseek_alig_matrix('./m8/%s.m8' % primary_list[4])
            result = result_matrix(primary_list[4], alignment_matrix, biomuta_profile, hotspot_profile,
                                   primary_list)
            wantedhot = hotspot_profile[(hotspot_profile
                                         .apply(lambda r: r.astype('string').str.contains(i)
                                                .any(), axis=1))]
            wantedmuta = biomuta_profile[(biomuta_profile
                                          .apply(lambda r: r.astype('string').str.contains(i)
                                                 .any(), axis=1))]
            reu_compelet.append({'input_protein': primary_list[0],
                                 'o_biomuta': wantedmuta['biomuta_profile'].tolist()[0],
                                 'o_hotspot': wantedhot['hotspot_profile'].tolist()[0],
                                 'p_positions': result[2]})
        except:
            wantedhot = hotspot_profile[(hotspot_profile
                                         .apply(lambda r: r.astype('string').str.contains(i)
                                                .any(), axis=1))]
            wantedmuta = biomuta_profile[(biomuta_profile
                                          .apply(lambda r: r.astype('string').str.contains(i)
                                                 .any(), axis=1))]
            reu_compelet.append({'input_protein': primary_list[0],
                                 'o_biomuta': wantedmuta['biomuta_profile'].tolist()[0],
                                 'o_hotspot': wantedhot['hotspot_profile'].tolist()[0],
                                 'p_positions': None})
    reu_compelet2 = pd.DataFrame(reu_compelet)
    reu_compelet2.to_pickle('./reu_complete_mix.pkl')

    return reu_compelet2

# rcsb = pd.read_pickle('./reu_complete_RCSB.pkl')
# alphe = pd.read_pickle('./reu_complete_ALFA.pkl')
# mix = pd.read_pickle('./reu_complete_mix.pkl')
#
#
# s = mix['p_positions'].tolist()
# s = [[] if i is None else i for i in s]
# mix['m_positions'] = s

# selection2 = []
# for i in alphe.iloc:
#     o_hotspot = i['o_hotspot']
#     o_biomuta = i['o_biomuta']
#     p_hotspot = i['p_positions']
#     e_hotspot = i['e_positions']
#     pro_name = i['input_protein']
#     len_prot = protein_length(pro_name)
#     if len(o_biomuta) != 0:
#         in_hot_e = []
#         in_bio_e = []
#         non_mut_e = []
#         in_hot_p = []
#         in_bio_p = []
#         non_mut_p = []
#         for j in p_hotspot:
#             if j in o_biomuta and j not in o_hotspot:
#                 in_bio_p.append(j)
#             if j in o_hotspot:
#                 in_hot_p.append(j)
#             if j not in o_biomuta and j not in o_hotspot:
#                 non_mut_p.append(j)
#         ptp = len(in_bio_p) + len(in_hot_p)
#         pfp = len(non_mut_p)
#         pfn = len(o_biomuta) - ptp
#         ptn = (len_prot - len(o_biomuta)) - len(non_mut_p)
#         for k in e_hotspot:
#             if k in o_biomuta and k not in o_hotspot:
#                 in_bio_e.append(k)
#             if k in o_hotspot:
#                 in_hot_e.append(k)
#             if k not in o_biomuta and k not in o_hotspot:
#                 non_mut_e.append(k)
#         etp = len(in_bio_e) + len(in_hot_e)
#         efp = len(non_mut_e)
#         efn = len(o_biomuta) - etp
#         etn = (len_prot - len(o_biomuta)) - len(non_mut_e)
#         selection2.append(
#             {'protein': pro_name, 'Protein_len': len_prot, 'predicted_p': p_hotspot, 'predicted_e': e_hotspot,
#              'in_bio_p': in_bio_p, 'in_hot_p': in_hot_p, 'non_mut_p': non_mut_p, 'in_bio_e': in_bio_e,
#              'in_hot_e': in_hot_e, 'non_mut_e': non_mut_e, 'ptp': ptp, 'pfp': pfp, 'pfn': pfn, 'ptn': ptn, 'etp': etp,
#              'efp': efp, 'efn': efn, 'etn': etn})
#     else:
#         selection2.append(
#             {'protein': pro_name, 'Protein_len': len_prot, 'predicted_p': None, 'predicted_e': None, 'in_bio_p': None,
#              'in_hot_p': None, 'non_mut_p': None, 'in_bio_e': None, 'in_hot_e': None, 'non_mut_e': None, 'ptp': 0,
#              'pfp': 0, 'pfn': 0, 'ptn': 0, 'etp': 0, 'efp': 0, 'efn': 0, 'etn': 0})
#
# selection2 = pd.DataFrame(selection2)
# selection2.to_pickle('./pkl/acc_p_e.pkl')
# print(selection2)

# acc_p_e = pd.read_pickle('./pkl/acc_p_e.pkl')
# tp_p = sum(acc_p_e['ptp'])
# fp_p = sum(acc_p_e['pfp'])
# tn_p = sum(acc_p_e['ptn'])
# fn_p = sum(acc_p_e['pfn'])
# sensitivity_p = tp_p/(tp_p+fn_p) #recall
# specificty_p = tn_p/(tn_p+fp_p)
# accurecy_p = (tp_p+tn_p)/(tp_p+tn_p+fp_p+fn_p)
# precision_p = tp_p/(tp_p+fp_p)
# f_measure_p = (2*precision_p*sensitivity_p)/(precision_p+sensitivity_p)
# tp_e = sum(acc_p_e['etp'])
# fp_e = sum(acc_p_e['efp'])
# tn_e = sum(acc_p_e['etn'])
# fn_e = sum(acc_p_e['efn'])
# sensitivity_e = tp_e/(tp_e+fn_e) #recall
# specificty_e = tn_e/(tn_e+fp_e)
# accurecy_e = (tp_e+tn_e)/(tp_e+tn_e+fp_e+fn_e)
# precision_e = tp_e/(tp_e+fp_e)
# f_measure_e = (2*precision_e*sensitivity_e)/(precision_e+sensitivity_e)
# print('the end')
# from predicted structure we got 6 returned hotspots

# selection3 = []
# for i in mix.iloc:
#     o_hotspot = i['o_hotspot']
#     o_biomuta = i['o_biomuta']
#     p_hotspot = i['m_positions']
#     pro_name = i['input_protein']
#     len_prot = protein_length(pro_name)
#     if len(o_biomuta) != 0:
#         in_hot_p = []
#         in_bio_p = []
#         non_mut_p = []
#         for j in p_hotspot:
#             if j in o_biomuta and j not in o_hotspot:
#                 in_bio_p.append(j)
#             if j in o_hotspot:
#                 in_hot_p.append(j)
#             if j not in o_biomuta and j not in o_hotspot:
#                 non_mut_p.append(j)
#         ptp = len(in_bio_p) + len(in_hot_p)
#         pfp = len(non_mut_p)
#         pfn = len(o_biomuta) - ptp
#         ptn = (len_prot - len(o_biomuta)) - len(non_mut_p)
#         selection3.append(
#             {'protein': pro_name, 'Protein_len': len_prot, 'predicted_m': p_hotspot,
#              'in_bio_m': in_bio_p, 'in_hot_m': in_hot_p, 'non_mut_m': non_mut_p, 'mtp': ptp, 'mfp': pfp, 'mfn': pfn, 'mtn': ptn})
#     else:
#         selection3.append(
#             {'protein': pro_name, 'Protein_len': len_prot, 'predicted_m': p_hotspot,
#              'in_bio_m': [], 'in_hot_m': [], 'non_mut_m': [], 'mtp': 0, 'mfp': 0, 'mfn': 0,
#              'mtn': 0})
#
# selection3 = pd.DataFrame(selection3)
# selection3.to_pickle('./pkl/acc_m.pkl')
# acc_m = pd.read_pickle('./pkl/acc_m.pkl')
# tp_p = sum(acc_m['mtp'])
# fp_p = sum(acc_m['mfp'])
# tn_p = sum(acc_m['mtn'])
# fn_p = sum(acc_m['mfn'])
# sensitivity_p = tp_p/(tp_p+fn_p) #recall
# specificty_p = tn_p/(tn_p+fp_p)
# accurecy_p = (tp_p+tn_p)/(tp_p+tn_p+fp_p+fn_p)
# precision_p = tp_p/(tp_p+fp_p)
# f_measure_p = (2*precision_p*sensitivity_p)/(precision_p+sensitivity_p)
