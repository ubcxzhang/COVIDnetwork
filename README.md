# Integrative COVID-19 Biological Network Inference with Probabilistic Core Decomposition

This repo contains two parts of program: data preprocessing and peeling algorithm (PA).

## data preprocessing

Source data availability: human Biomine database can be found at [biomine.ijs.si](https://biomine.ijs.si/) and for SARS-CoV-2-human protein-protein interaction network please check [nature.com/articles/s41586-020-2286-9#data-availability](https://www.nature.com/articles/s41586-020-2286-9#data-availability).

`BiomineCOVID_data_preprocessing.py` will generate a duplicated/loop-free version of [Human Biomine + PPI network] dataset for PA. Before in Biomine database/PPI network, nodes are `string` named, they are now renamed to `integer`. 

The program outputs three important files:
- `Biomine_data_unique_withCOVID_renamed.txt`: dataset to be fed into PA
- `node_name_mapping_Biomine_data_unique_withCOVID_Ready4Stage1.pickle`: node name mapping (str -> int)
- `node_to_skip.txt`: retained nodes

## PA

The programs are already compiled with javac 1.8 (Java 8, 1.8.0_144).

`Biomine_data_unique_withCOVID_renamed.txt` should have the following format (from  to  prob):

```
0	703991	0.511
0	704014	0.6729999999999999
0	704125	0.5660000000000001
0	704129	0.6729999999999999
0	820845	0.608
1	761889	0.8
1	818497	0.8
2	1	0.8
3	704015	0.8
3	768047	0.608
...
```

Run PA using the following script:

```sh
mv Biomine_data_unique_withCOVID_renamed.txt PA_run/pcore
mv node_to_skip.txt PA_run/pcore
cd PA_run/pcore
pwd
touch Biomine_data_unique_withCOVID_renamed-proc.txt
echo "TTextProc Program start"
date
java -cp "out/production/pcore:bin:lib/*" TTextProc Biomine_data_unique_withCOVID_renamed.txt Biomine_data_unique_withCOVID_renamed-proc.txt
echo "TTextProc Program end"
date
echo "Remaining Preprocessing Program start"
date
cut -f 1-2 Biomine_data_unique_withCOVID_renamed-proc.txt > edgelistfile.txt
sort -nk 1 edgelistfile.txt | uniq > edgelistsortedfile.txt
java -cp "lib/*" it.unimi.dsi.webgraph.BVGraph -1 -g ArcListASCIIGraph dummy Biomine_data_unique_withCOVID_renamed-proc < edgelistsortedfile.txt
java -cp "out/production/pcore:bin:lib/*" GenerateWeightedGraphFromTxtLong Biomine_data_unique_withCOVID_renamed-proc Biomine_data_unique_withCOVID_renamed-proc.txt 17
echo "Remaining Preprocessing Program end"
date
echo "PA Program start"
date
java -cp "out/production/pcore:bin:lib/*" K_BZ_new_biominProject_COVID
echo "PA Program end"
date
```

The core decomposition result will be in `Biomine_data_unique_withCOVID_renamed-proc.weta-0.5-bz.txt` with the following format (node_name  coreness):

```
2	52
3	8
4	77
5	77
6	77
7	77
8	35
9	77
10	77
11	60
12	62
13	35
14	77
...
```