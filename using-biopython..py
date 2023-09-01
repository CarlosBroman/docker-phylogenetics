from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


forward_path = "./my_chromatograms/MH3_F.ab1"
reverse_path = "./my_chromatograms/MH3_R.ab1"
reference_path = "./sequence.gb"

_atgc_dict = {0:"A", 1:"T", 2:"G", 3:"C"}

# Convert abi file into readable values
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

def find_common_area(seq1, seq2, max_consecutive_gaps=5):
    common_area = ""
    consecutive_letters_count = 0
    consecutive_gaps_count = 0

    while True:
        for char1, char2 in zip(seq1, seq2):
            if char1 != '-' and char2 != '-':
                common_area += char1
                consecutive_letters_count += 1
                consecutive_gaps_count = 0
            else:
                common_area += '-'
                consecutive_letters_count = 0
                consecutive_gaps_count += 1

            if consecutive_gaps_count >= max_consecutive_gaps:
                if len(common_area) > 50:
                    break
                else:
                    common_area=""
                    consecutive_letters_count = 0
                    consecutive_gaps_count = 0

        return common_area


# Generate sequences (forward / reverse)TCGAGCACCTATGAAGGGCGCAGCGAAGTGCGATAATCGTTGTGAATTGCAGAACTCCGTGAACCAATGGCCT
def generate_consensusseq(abidata):
    consensus_seq = ""  
    for values in zip(abidata["channel"]["A"][1], abidata["channel"]["T"][1], abidata["channel"]["G"][1], abidata["channel"]["C"][1]):
        consensus_seq += _atgc_dict[values.index(max(values))]

    return (consensus_seq) 

genbank_record = SeqIO.read(reference_path, 'genbank')

forward_seq = generate_consensusseq(abi_to_dict(forward_path))
reverse_seq = generate_consensusseq(abi_to_dict(reverse_path))
genbank_seq = str(genbank_record.seq)

forward = Seq(forward_seq)
reverse = Seq(reverse_seq).reverse_complement()

# Perform a pairwise alignment
alignments_forward_reverse = pairwise2.align.globalxx(forward, reverse, one_alignment_only=True, gap_char='-')
aligned_forward = alignments_forward_reverse[0][0]
aligned_reverse = alignments_forward_reverse[0][1]

# Find the common region allowing for up to 2 mismatches (adjust max_mismatches as needed)
common_area = find_common_area(aligned_forward, aligned_reverse, max_consecutive_gaps=5)

alignments_forward_reverse = pairwise2.align.globalxx(aligned_forward, common_area, one_alignment_only=True, gap_char='-')
common_region = alignments_forward_reverse[0][1]

# Print the common region
print("Aligned Common Region:\n", common_region)

# Print the aligned sequences
print("Aligned Forward Sequence:\n", aligned_forward)
print("Aligned Reverse Sequence:\n", aligned_reverse)