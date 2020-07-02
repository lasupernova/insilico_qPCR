# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 16:15:07 2019

@author: asenturk16

Modified on Thu Jul 2 19:51:02 2020

@contributor: Gabriela K. Paulus
"""
import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from difflib import SequenceMatcher
from datetime import datetime

start_time=datetime.now() #start recording time

#read subgene sequence file
df=pd.read_excel("amplicon_sequences.xlsx").set_index("name").replace('\n','', regex=True)
#read file containing primer information
primer=pd.read_excel("primer_sequences.xlsx").set_index("gene")
#create complementary and inverted sequence for primers and add to primer file
bases=[("A","T"),("T","A"),("G","C"),("C","G")]
for i,rows in primer.iterrows():
    forw=primer["forward 5'-3'"][i]
    rev=primer["reverse 5'-3'"][i]
    #generate complementary reverse sequence
    comp_rev=""
    for r in rev:
        for t in np.arange(0,4,1):
            if r==bases[t][0]:
                comp_rev+=bases[t][1]
    primer.loc[i,"reverse complementary & inverted"]=comp_rev[::-1] #insert for row i the created comp_rev string in new column
        
#check number of lines for files

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


for s in samples:
    print("Currently processing sample: "+s)
    num_lines_r1=sum(1 for line in open(s))
    num_lines_r2=sum(1 for line in open(s))

    if num_lines_r1==num_lines_r2: #check if length of both NGS files is equal
        #load first NGS file
        r1=open(s)
        #load second NGS file
        r2=open(s)
        #set name of experiment
        name=r1.name[:-28] #letters until LOO1
        #in the beginning second line of each file should be used, afterwards every 4th line until end of file
        iterations=[2]+([4]*(int(num_lines_r1/4)))
        
        hits=[] #create empty list to collect sequences where primers align
        line_count=0  #count the number of lines successfully iterated
        
        for x in iterations:       
            #read first NGS file line by line until iteration setting is true
            lines=0
            for line in r1:
                lines+=1 #increase lines number after each line
                if lines==x:
                    read1=line.replace("\n","")
                    print("read1 sequence: "+read1)
                    break
                    
            #read second NGS file line by line until iteration setting is true
            lines=0
            for line in r2:
                lines+=1
                if lines==x:
                    read2=line.replace("\n","")
                    print("read2 sequence: "+read2)
                    break     
            
            #create inverted complementary sequence of read2
            comp_read2=""
            for r in read2:
                for t in np.arange(0,4,1):
                    if r==bases[t][0]:
                        comp_read2+=bases[t][1]
                    comp_inv_read2=comp_read2[::-1]
            print("read2_inv_comp sequence: "+comp_inv_read2)
            
            ##create concatenated sequence
            #if reads are the same
            if read1==comp_inv_read2:
                sequence=read1 #then read1 is whole sequence and no further sequence info gained by read2
            #if reads share common sequence part
            else:
                match=SequenceMatcher(None, read1, comp_inv_read2).find_longest_match(0, len(read1), 0, len(comp_inv_read2))
                match_seq=read1[match.a:(match.a + match.size)] #common sequence between read1 and comp_inv_read2
                if (read1.endswith(match_seq)) & (comp_inv_read2.startswith(match_seq)): #a true intersection between read sequences means that the common sequence is at the end of read1 and at the same time at the start of comp_inv_read2
                    sequence=read1+comp_inv_read2[len(match_seq):] #the concatenated sequence is the full read1 sequence + the part of comp_inv_read2 after the common sequence               
                else: #if the common sequence only a short random sequence appearing in the middle of the read sequences then this read will be omitted
                    line_count+=x
                    continue #no concatenated sequence, go on to next reads
            
            ##check if primers are in concatenated sequence
            for i,rows in primer.iterrows():
                forw=primer["forward 5'-3'"][i]
                rev_inv_comp=primer["reverse complementary & inverted"][i]
            
                #both primers in sequence
                if (forw in sequence) & (rev_inv_comp in sequence): 
                    gene_by_primer_file=i #gene according to primer file
                    aligned_primer_type="forward & rev_inv_comp" #which primer aligned to concatenated sequence
                    aligned_seq=forw+sequence.partition(forw)[2].partition(rev_inv_comp)[0]+rev_inv_comp #sequence btw forw and rev_inv_comp in concat sequence
                    if df["seq"].str.contains(aligned_seq).any(): #check if kwr_to_card file contains aligned sequence
                        gene_by_subgene_file=[z for z in df[df["seq"].str.contains(aligned_seq)]["gene"]] #if in multiple subgenes, gene names will be collected in list
                    else:
                        gene_by_subgene_file="" #if aligned sequence not in kwr_to_card file then no entry
                    hits.append((sequence,aligned_primer_type,gene_by_primer_file,aligned_seq,gene_by_subgene_file)) #collect all important results if there is alignment of primer
                #only forward primer in sequence
                elif (forw in sequence) & (rev_inv_comp not in sequence): 
                    gene_by_primer_file=i #gene according to primer file
                    aligned_primer_type="forward" #which primer aligned to concat sequence
                    aligned_seq=forw+sequence.partition(forw)[2]
                    if df["seq"].str.contains(aligned_seq).any(): #check if kwr_to_card file contains aligned sequence
                        gene_by_subgene_file=[z for z in df[df["seq"].str.contains(aligned_seq)]["gene"]] #if in multiple subgenes, gene names will be collected in list
                    else:
                        gene_by_subgene_file="" #if aligned sequence not in kwr_to_card file then no entry
                    hits.append((sequence,aligned_primer_type,gene_by_primer_file,aligned_seq,gene_by_subgene_file)) #collect all important results if there is alignment of primer
                #only reverse primer in sequence 
                elif (rev_inv_comp in sequence) & (forw not in sequence): 
                    gene_by_primer_file=i #gene according to primer file
                    aligned_primer_type="rev_inv_comp" #which primer aligned to concat sequence
                    aligned_seq=sequence.partition(rev_inv_comp)[0]+rev_inv_comp
                    if df["seq"].str.contains(aligned_seq).any(): #check if kwr_to_card file contains aligned sequence
                        gene_by_subgene_file=[z for z in df[df["seq"].str.contains(aligned_seq)]["gene"]] #if in multiple subgenes, gene names will be collected in list
                    else:
                        gene_by_subgene_file="" #if aligned sequence not in kwr_to_card file then no entry
                    hits.append((sequence,aligned_primer_type,gene_by_primer_file,aligned_seq,gene_by_subgene_file)) #collect all important results if there is alignment of primer
                #no primer alignment
                else:
                    continue #go on to next primer
            
            line_count+=x
        
        end_time=datetime.now()
        print("Duration: {}".format(end_time - start_time))
        text_file=open("Lines_newprimer_"+name+".txt", "w")
        text_file.write("total number of lines: "+str(num_lines_r1)+"\n")
        text_file.write("number of executed lines: "+str(line_count)+"\n")
        text_file.write("duration: {}".format(end_time - start_time))
        text_file.close()
        pd.DataFrame(hits).to_excel("NGS_output_test_"+name+".xlsx", header=["concatenated sequence","primer that aligned","gene by primer file","aligned sequence","gene by subgene file"], index=False)
        
    else:
        print("Number of lines in R1 and R2 are not equal!")
    
            
                
        
  



