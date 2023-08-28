import os 
import sys 
import copy
import argparse
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib
matplotlib.rcParams["xtick.major.width"] = 0.4
matplotlib.rcParams["xtick.minor.width"] = 0.4
matplotlib.rcParams["ytick.major.width"] = 0.4
matplotlib.rcParams["ytick.minor.width"] = 0.4

_atgc_dict = {0:"A", 1:"T", 2:"G", 3:"C"}
color_dict = {"G":"#f2f059", "C":"#74b2d7", "A":"#79E5B7", "T":"#ff776c", "N":"#FFFFFF", "-":"#FFFFFF", "/":"#FFFFFF"}



def abi_to_dict(filename):
    record   = SeqIO.read(filename,'abi')
    abi_data = {"conf":[[],[]],
                "channel": {
                           "A":[[],[]],
                           "T":[[],[]],
                           "G":[[],[]],
                           "C":[[],[]],
                           },
                "_channel":{
                            "A":[[],[]],
                            "T":[[],[]],
                            "G":[[],[]],
                            "C":[[],[]], 
                           }
                }
    
    pre_pos = 0 
    pos_set = [] 
    for i, (pos, conf) in enumerate(zip(record.annotations['abif_raw']["PLOC1"], record.annotations['abif_raw']["PCON1"])):
        abi_data["conf"][0].append(i)
        abi_data["channel"]["G"][0].append(i)
        abi_data["channel"]["A"][0].append(i)
        abi_data["channel"]["T"][0].append(i)
        abi_data["channel"]["C"][0].append(i) 
        
        abi_data["conf"][1].append(conf)
        abi_data["channel"]["G"][1].append(record.annotations['abif_raw']["DATA9"][pos])
        abi_data["channel"]["A"][1].append(record.annotations['abif_raw']["DATA10"][pos])
        abi_data["channel"]["T"][1].append(record.annotations['abif_raw']["DATA11"][pos])
        abi_data["channel"]["C"][1].append(record.annotations['abif_raw']["DATA12"][pos])    
       
        step = 0.1
        for j in range(10):
            abi_data["_channel"]["G"][0].append(i+step*j)   
            abi_data["_channel"]["A"][0].append(i+step*j)
            abi_data["_channel"]["T"][0].append(i+step*j) 
            abi_data["_channel"]["C"][0].append(i+step*j)   
            try:
                abi_data["_channel"]["G"][1].append(record.annotations['abif_raw']["DATA9"][pos-5+j])   
                abi_data["_channel"]["A"][1].append(record.annotations['abif_raw']["DATA10"][pos-5+j])
                abi_data["_channel"]["T"][1].append(record.annotations['abif_raw']["DATA11"][pos-5+j]) 
                abi_data["_channel"]["C"][1].append(record.annotations['abif_raw']["DATA12"][pos-5+j])  
            except:
                abi_data["_channel"]["G"][1].append(0)   
                abi_data["_channel"]["A"][1].append(0)
                abi_data["_channel"]["T"][1].append(0) 
                abi_data["_channel"]["C"][1].append(0)  

        pos_set.append((pre_pos, pos))
        pre_pos = pos
    pos_set.append((pre_pos, record.annotations['abif_raw']["PLOC1"][-1]))
    
    return abi_data

def generate_consensusseq(abidata):
    consensus_seq = ""  
    for values in zip(abidata["channel"]["A"][1], abidata["channel"]["T"][1], abidata["channel"]["G"][1], abidata["channel"]["C"][1]):
        max_value = max(values)
        max_index = values.index(max_value)
        consensus_seq += _atgc_dict[max_index]
        position = len(consensus_seq) - 1
        base = consensus_seq[position]
        print(f"Position: {position}, Base: {base}, Max Value: {max_value}, Values: {values}, Bases: {_atgc_dict}")
     
    return (consensus_seq, consensus_seq.translate(str.maketrans("ATGC","TACG"))[::-1]) 







# TRIALS
filename = "my_chromatograms/MH3_F.ab1"

with open("my_chromatograms/MH3_F.txt", "r") as sequence_file:
    original_txt = sequence_file.read().replace("\n", "").strip()



abi_data = abi_to_dict(filename)

consensus_sequences = generate_consensusseq(abi_data)

# Unpack the results of the generate_consensusseq function
consensus_seq, reverse_complement_seq = consensus_sequences

# Print or use the generated consensus sequences
print("Consensus Sequence:", consensus_seq)

alignments = pairwise2.align.globalxx(consensus_seq, original_txt)
print("the alignment")
print(format_alignment(*alignments[0]))