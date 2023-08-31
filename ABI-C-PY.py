from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt

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

    base_dict = {"G": 0, "A": 1, "T": 2, "C": 3}

    for i, (pos, conf) in enumerate(zip(record.annotations['abif_raw']["PLOC1"], record.annotations['abif_raw']["PCON1"])):
        abi_data["conf"][0].append(i)
        for base in base_dict:
            abi_data["channel"][base][0].append(i)
            abi_data["channel"][base][1].append(record.annotations['abif_raw'][f"DATA{9 + base_dict[base]}"][pos])
            
            for j in range(10):
                abi_data["_channel"][base][0].append(i + 0.1 * j)
                try:
                    abi_data["_channel"][base][1].append(record.annotations['abif_raw'][f"DATA{9 + base_dict[base]}"][pos - 5 + j])
                except:
                    abi_data["_channel"][base][1].append(0) 


        pos_set.append((pre_pos, pos))
        pre_pos = pos
    pos_set.append((pre_pos, record.annotations['abif_raw']["PLOC1"][-1]))
    return abi_data

def generate_values(abidata):
    intensity_values = []  # To store the maximum values

    for values in zip(abidata["channel"]["A"][1], abidata["channel"]["T"][1], abidata["channel"]["G"][1], abidata["channel"]["C"][1]):
        intensity_values.append(values)
        print(values)

    return (intensity_values)

# def plot_intensities(intensities):
    plt.figure(figsize=(8, 6))  # Create a new figure

    start_index = 50  # Replace with the starting index of your desired slice
    end_index = 70 

    sliced_positions = positions[start_index:end_index]
    sliced_values = [values[start_index:end_index] for values in [abidata["channel"]["A"][1], abidata["channel"]["T"][1], abidata["channel"]["G"][1], abidata["channel"]["C"][1]]]
    base_colors = ['blue', 'green', 'red', 'purple']  # Define colors for each base

    # Create separate line plots for each base
    for i, base_values in enumerate(sliced_values):
        base_label = _atgc_dict[i]  # Get the base label
        base_color = base_colors[i]  # Get the corresponding color
        plt.plot(sliced_positions, base_values, label=base_label, color=base_color, linewidth=2)

    # Add labels and title
    plt.xlabel("Base Position")
    plt.ylabel("Values")
    plt.title("Values vs. Base Positions for Different Bases")
    plt.legend()  # Show legend with base labels

    plt.grid(True)  # Add grid lines

    plt.show()  # Display the plot

def generate_consensusseq(abidata):
    consensus_seq = ""
    
    for values in generate_values(abidata):
        max_value = max(values)
        max_index = values.index(max_value)
        consensus_seq += _atgc_dict[max_index]
        complementary_seq = consensus_seq.translate(str.maketrans("ATGC","TACG"))[::-1]

    return (consensus_seq, complementary_seq) 


# TRIALS
filename = "my_chromatograms/MH3_F.ab1"

with open("my_chromatograms/MH3_F.txt", "r") as sequence_file:
    original_txt = sequence_file.read().replace("\n", "").strip()

abi_data = abi_to_dict(filename)

consensus_sequences = generate_consensusseq(abi_data)

# Unpack the results of the generate_consensusseq function
consensus_seq, reverse_complement_seq = consensus_sequences
alignments = pairwise2.align.globalxx(consensus_seq, original_txt)

alignment = format_alignment(*alignments[0])

print(f'\nALIGNMENT \n {alignment}')