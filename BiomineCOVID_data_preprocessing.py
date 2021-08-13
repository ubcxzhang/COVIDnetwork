# coding: utf-8

import pandas as pd
import pickle
import math
import numpy as np
from datetime import datetime
import gc
import os


if __name__ == '__main__':
    print('Data pre-processing started:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '--------------------------------------------------------------------------------------')
    
    # load Biomine original dataset
    Biomine_data = pd.read_csv('human_3_oct_2018.bmg', sep=" ", header=None, skiprows=6)
    Biomine_data.columns = ["from","to","relation","weightProb"]
    
    # load COVID network dataset and node mapping
    COVID_proteinGene333_mapped_complete_BiomineNamingFormat = pd.read_csv('COVID_proteinGene333_mapped_complete_BiomineNamingFormat.csv', 
        header=None, index_col=0, squeeze=True, sep=',').to_dict()
    COVID_viral_gene_mapped = pd.read_csv('COVID_viral_gene_mapped.csv', header=None, index_col=0, squeeze=True, sep='\t').to_dict()
    COVID_network = pd.read_csv('COVID_network_edgelist.csv', header=None, sep='\t')
    geneList_by_cell_paper_mapped_uniprot = list(pd.read_csv("geneList_by_cell_paper_mapped.newer.uniprot.csv", header=None)[0])
    print('len of COVID_proteinGene333_mapped_complete_BiomineNamingFormat dict:', len(COVID_proteinGene333_mapped_complete_BiomineNamingFormat))
    print('len of COVID_viral_gene_mapped dict:', len(COVID_viral_gene_mapped))
    print('len of geneList_by_cell_paper_mapped_uniprot (should be 439):', len(geneList_by_cell_paper_mapped_uniprot))
    
    ## merge COVID_viral_gene_mapped and COVID_proteinGene333_mapped_complete_BiomineNamingFormat, since pandas will replace those not in dict with NaN
    COVID_viral_gene_mapped.update(COVID_proteinGene333_mapped_complete_BiomineNamingFormat)
    print('len of COVID_viral_gene_mapped dict after merge (check if == 360):', len(COVID_viral_gene_mapped)) # should be 333+27
    
    ## ====== NOTE that now COVID_viral_gene_mapped is complete map that combined COVID_proteinGene333_mapped_complete_BiomineNamingFormat
    ## map COVID network to proper naming (with 'COVID_' prefix and replace with Biomine Formatted EntrezID, UniProt, etc.)
    COVID_network_mapped = COVID_network.copy()
    COVID_network_mapped.columns = ["from","to"]
    print("Original COVID network ===============================================")
    print(COVID_network_mapped.head(10))
    COVID_network_mapped['from'] = COVID_network_mapped['from'].map(COVID_viral_gene_mapped)
    COVID_network_mapped['to'] = COVID_network_mapped['to'].map(COVID_viral_gene_mapped)
    print("Reformatted COVID network ===============================================")
    print(COVID_network_mapped.head(10))
    print("=========================================================================")
    with open('COVID_viral_and_human_node_mapped_dict_All360node.pickle', 'wb') as f4:
        pickle.dump(COVID_viral_gene_mapped, f4)
    with open('COVID_network_mapped_dataFrame_All695edge.pickle', 'wb') as f4:
        pickle.dump(COVID_network_mapped, f4)
    
    COVID_network_mapped_dataFrame_All695edge_beforeStage1 = COVID_network_mapped
    COVID_network_mapped_dataFrame_All695edge_beforeStage1["relation"] = 'COVID_before_stage1'
    COVID_network_mapped_dataFrame_All695edge_beforeStage1["weightProb"] = 1.0
    with open('COVID_network_mapped_dataFrame_All695edge_beforeStage1.pickle', 'wb') as f4:
        pickle.dump(COVID_network_mapped_dataFrame_All695edge_beforeStage1, f4)
    
    
    print('Start merging duplicated rows:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    # Merge the duplicated rows and use mean as new weightProb
    Biomine_data_unique = Biomine_data.groupby(['from','to']).agg({'relation': ', '.join,
                                                                'weightProb': 'mean' }).reset_index()
    with open('Biomine_data_unique.pickle', 'wb') as f2:
        pickle.dump(Biomine_data_unique, f2)
    
    ## find level 1 connection of COVID network
    human_node_in_COVID_network = list(COVID_proteinGene333_mapped_complete_BiomineNamingFormat.values())
    print("len of human_node_in_COVID_network:", len(human_node_in_COVID_network))
    print("how many human node in PPI network exist in geneList_by_cell_paper_mapped_uniprot:", len([i for i in human_node_in_COVID_network if i in geneList_by_cell_paper_mapped_uniprot]))
    human_node_in_COVID_network = [*geneList_by_cell_paper_mapped_uniprot, *human_node_in_COVID_network]
    print("len of human_node_in_COVID_network NOW:", len(human_node_in_COVID_network))
    human_node_in_COVID_network = list(set(human_node_in_COVID_network))
    print("len of human_node_in_COVID_network NOW after remove duplicates:", len(human_node_in_COVID_network))
    COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_edgelist = Biomine_data_unique[Biomine_data_unique[['from','to']].isin(human_node_in_COVID_network).any(1)]
    COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_nodelist = list(pd.concat([COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_edgelist['from'], COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_edgelist['to']]).unique())
    COVIDviral_nodelist = ['COVID_E','COVID_M','COVID_N','COVID_S','COVID_nsp1','COVID_nsp10','COVID_nsp11','COVID_nsp12',
    'COVID_nsp13','COVID_nsp14','COVID_nsp15','COVID_nsp2','COVID_nsp4','COVID_nsp5','COVID_nsp5_C145A','COVID_nsp6',
    'COVID_nsp7','COVID_nsp8','COVID_nsp9','COVID_orf10','COVID_orf3a','COVID_orf3b','COVID_orf6','COVID_orf7a','COVID_orf8','COVID_orf9b','COVID_orf9c']
    print('len of COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_edgelist:', len(COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_edgelist))
    print('len of COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_nodelist:', len(COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_nodelist))
    print('len of COVIDviral_nodelist:', len(COVIDviral_nodelist))
    
    node_to_skip_in_all_stages = [*geneList_by_cell_paper_mapped_uniprot, *COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_nodelist, *human_node_in_COVID_network, *COVIDviral_nodelist]
    print('len of node_to_skip_in_all_stages with duplicates (human_node_in_COVID_network + COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_nodelist + COVIDviral_nodelist)', len(node_to_skip_in_all_stages))
    ## remove duplicated nodes in node_to_skip_in_all_stages
    node_to_skip_in_all_stages = list(set(node_to_skip_in_all_stages))
    print('len of node_to_skip_in_all_stages after removing duplicates', len(node_to_skip_in_all_stages))
    with open('node_to_skip_in_all_stages.pickle', 'wb') as f2:
        pickle.dump(node_to_skip_in_all_stages, f2)
    with open('COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_edgelist.pickle', 'wb') as f2:
        pickle.dump(COVIDhuman_node_andlevel1_node_in_Biomine_data_unique_edgelist, f2)

    ## if we look for more than level 1 connection but also level 2
    COVIDhuman_node_andlevel1and2_node_in_Biomine_data_unique_edgelist = Biomine_data_unique[Biomine_data_unique[['from','to']].isin(node_to_skip_in_all_stages).any(1)]
    COVIDhuman_node_level1and2_node_in_Biomine_data_unique_nodelist = list(pd.concat([COVIDhuman_node_andlevel1and2_node_in_Biomine_data_unique_edgelist['from'], COVIDhuman_node_andlevel1and2_node_in_Biomine_data_unique_edgelist['to']]).unique())
    print('len of COVIDhuman_node_andlevel1and2_node_in_Biomine_data_unique_edgelist:', len(COVIDhuman_node_andlevel1and2_node_in_Biomine_data_unique_edgelist))
    print('len of COVIDhuman_node_level1and2_node_in_Biomine_data_unique_nodelist:', len(COVIDhuman_node_level1and2_node_in_Biomine_data_unique_nodelist))
    with open('COVIDhuman_node_andlevel1and2_node_in_Biomine_data_unique_edgelist.pickle', 'wb') as f2:
        pickle.dump(COVIDhuman_node_andlevel1and2_node_in_Biomine_data_unique_edgelist, f2)
    with open('COVIDhuman_node_level1and2_node_in_Biomine_data_unique_nodelist.pickle', 'wb') as f2:
        pickle.dump(COVIDhuman_node_level1and2_node_in_Biomine_data_unique_nodelist, f2)

    # merge Biomine unique and COVID network
    Biomine_data_unique_withCOVID = Biomine_data_unique.append(COVID_network_mapped_dataFrame_All695edge_beforeStage1)
    Biomine_data_unique_withCOVID = Biomine_data_unique_withCOVID.reset_index(drop = True) # drop means replacing old index and not add old index as new col
    
    ## still do merge to account for possible overlap in Biomine_unique and COVID_network
    ## BUT we no longer use mean, now we use max so that those Prob=1 will not be affected
    print('Start merging duplicated rows (in Biomine and COVID network): ', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    # Merge the duplicated rows and use mean as new weightProb
    Biomine_data_unique_withCOVID_almostReady4Stage1 = Biomine_data_unique_withCOVID.groupby(['from','to']).agg({'relation': ', '.join,
                                                                                'weightProb':'max' }).reset_index()
    
    ## add additional code to account for when A->B and B->A both exist
    print('Start merging looped rows (when A->B and B->A both exist):', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    Biomine_data_unique_withCOVID_almostReady4Stage1['check_string'] = Biomine_data_unique_withCOVID_almostReady4Stage1.apply(lambda row:''.join(sorted([row['from'], row['to']])), axis=1)
    # Biomine_data_unique_withCOVID_ready4Stage1[Biomine_data_unique_withCOVID_ready4Stage1.duplicated('check_string')]
    print('check_string column added:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    # with open('Biomine_data_unique_withCOVID_almostReady4Stage1_withcheck_stringColumn.pickle', 'wb') as f2:
    #     pickle.dump(Biomine_data_unique_withCOVID_almostReady4Stage1, f2)
    print('start grouping:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    Biomine_data_unique_withCOVID_Ready4Stage1 = Biomine_data_unique_withCOVID_almostReady4Stage1.groupby(['check_string']).agg({'from': lambda seriesItem: seriesItem.iloc[0], 'to': lambda seriesItem: seriesItem.iloc[0], 
                                                                                        'relation': ', '.join, 'weightProb':'max'}).reset_index()
    print('Removing check_string column:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    Biomine_data_unique_withCOVID_Ready4Stage1 = Biomine_data_unique_withCOVID_Ready4Stage1.drop(['check_string'], axis=1)
    print('len of Biomine_data_unique_withCOVID_Ready4Stage1 (remove self-loop):', len(Biomine_data_unique_withCOVID_Ready4Stage1))
    with open('Biomine_data_unique_withCOVID_Ready4Stage1.pickle', 'wb') as f2:
        pickle.dump(Biomine_data_unique_withCOVID_Ready4Stage1, f2)

    print(Biomine_data_unique_withCOVID_Ready4Stage1.head())
    
    tmp_idx = 0
    missed_node_idx = []
    for item in node_to_skip_in_all_stages:
        if item == None:
            missed_node_idx.append(tmp_idx)
            tmp_idx = tmp_idx+1
        else:
            tmp_idx = tmp_idx+1

    if len(missed_node_idx) > 0:
        print('[ERROR] possible data corruption, `None` was returned! Some of the nodes did not get mapped:')
        print([node_to_skip_in_all_stages[item] for item in missed_node_idx])
        print('Program exiting...')
        raise SystemExit
    else:
        print('Data looks good.')

    
    print('Mapping Biomine string nodes with number:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    ## map Biomine string nodes with number
    node_name_list = pd.concat([Biomine_data_unique_withCOVID_Ready4Stage1['from'], Biomine_data_unique_withCOVID_Ready4Stage1['to']]).unique()
    node_name_mapping = { node_name_list[i]: i for i in range(len(node_name_list)) }
    print("After merge and duplicates/loops removal, how many node left:", len(node_name_mapping))
    with open('node_name_mapping_Biomine_data_unique_withCOVID_Ready4Stage1.pickle', 'wb') as f4:
        pickle.dump(node_name_mapping, f4)
    Biomine_data_unique_withCOVID_Ready4Stage1['from'] = Biomine_data_unique_withCOVID_Ready4Stage1['from'].map(node_name_mapping)
    Biomine_data_unique_withCOVID_Ready4Stage1['to'] = Biomine_data_unique_withCOVID_Ready4Stage1['to'].map(node_name_mapping)
    with open('Biomine_data_unique_withCOVID_Ready4Stage1_mapped2Int.pickle', 'wb') as f2:
        pickle.dump(Biomine_data_unique_withCOVID_Ready4Stage1, f2)
    node_to_skip_in_all_stages_mapped2Int = [node_name_mapping[i] for i in node_to_skip_in_all_stages]
    assert len(node_to_skip_in_all_stages_mapped2Int) == len(node_to_skip_in_all_stages)
    with open('node_to_skip_in_all_stages_mapped2Int.pickle', 'wb') as f2:
        pickle.dump(node_to_skip_in_all_stages_mapped2Int, f2)
    
    # Output edgelist and node_to_skip list
    print('Writing files:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    ## Write processed raw Biomine to txt file
    Biomine_data_unique_withCOVID_Ready4Stage1[['from','to','weightProb']].to_csv('Biomine_data_unique_withCOVID_renamed.txt', index=False, sep='\t', header=False)
    with open('node_to_skip.txt', 'w') as f:
        for i in node_to_skip_in_all_stages_mapped2Int:
            f.write(str(i)+'\n')

    print('Data pre-processing finished:', datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '--------------------------------------------------------------------------------------')