from Bio import SeqIO
from Bio import pairwise2

import matplotlib.pyplot as plt
import patchworklib as pw
from QUEEN.queen import *

_atgc_dict = {0:"A", 1:"T", 2:"G", 3:"C"}
color_dict = {"G":"#f2f059", "C":"#74b2d7", "A":"#79E5B7", "T":"#ff776c", "N":"#FFFFFF", "-":"#FFFFFF", "/":"#FFFFFF"}

def colorbar(ax, color_dict, query, subject, char=False, fontsize=10, label=False, zero_position=0):
    if len(query) > 200:
        bars = ax.bar(list(range(len(query))), [1.0] * (len(query)), width=1.0, edgecolor="#BBBBBB", linewidth=0.1, align="edge",bottom=0.05)
    else:
        bars = ax.bar(list(range(len(query))), [1.0] * (len(query)), width=1.0, edgecolor="#BBBBBB", linewidth=0.1, align="edge",bottom=0.05)
    ax.set_xlim(0,len(query))
    ax.set_ylim(0,1.00)
    p = 0
    for bar, q, s in zip(bars, query, subject):
        color = color_dict[q]
        if q != s:
            #bar.set_edgecolor("#E83929")
            if char == True:
                if q == "/":
                    bar.set_facecolor("#BBBBBB")
                    ax.text(p+0.5,0.45,"-",va="center",ha="center",fontsize=fontsize,zorder=100,color="k")     
                else:
                    bar.set_facecolor("#993366")
                    bar.set_alpha(0.2)
                    if q == "-":
                        ax.text(p+0.5,0.45,"-",va="center",ha="center",fontsize=fontsize,zorder=100,color="k")#fontweight="bold")
                    else:
                        ax.text(p+0.5,0.45,q,va="center",ha="center",fontsize=fontsize,zorder=100,color="k")#fontweight="bold")
        else:
            bar.set_facecolor("#FFFFFF")
            if char == True:
                ax.text(p+0.5,0.45,q,va="center",ha="center",fontsize=fontsize,zorder=100)     
        p += 1 
    
    #ax.set_xticks([])
    if label == True:
        positions  = [pos - 0.5 for pos in range(0, len(query)+1, 10)]
        ticklabels = [str(int(zero_position+pos+0.5)) for pos in positions] 
        if positions[0] < 1 and zero_position == 0:
            positions[0]  = 0.5
            ticklabels[0] = 1
        ax.set_xticks(positions)
        ax.set_xticklabels(ticklabels)
    else:
        ax.set_xticks([])
    
    ax.set_zorder(2.0)
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_linewidth(0.4)  
    #ax.spines["top"].set_visible(False)
    ax.patch.set_alpha(1.0)
    ax.set_xlim(0, len(query))
    return bars

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
        consensus_seq += _atgc_dict[values.index(max(values))]
        reverse_seq = consensus_seq.translate(str.maketrans("ATGC","TACG"))[::-1]
     
    return (consensus_seq, reverse_seq)

def generate_values(abidata):
    intensity_values = []  # To store the maximum values
    positions = []
    consensus_seq = ""
    all_values= []

    for values in zip(abidata["channel"]["A"][1], abidata["channel"]["T"][1], abidata["channel"]["G"][1], abidata["channel"]["C"][1]):
        intensity_values.append(values)
        positions.append(len(intensity_values))
        max_index=values.index(max(values))
        consensus_seq += _atgc_dict[max_index]
        all_values.append(values)
    print(consensus_seq)
    print(positions)

    plt.figure(figsize=(8, 6))  # Create a new figure

    start_index = 50  # Replace with the starting index of your desired slice
    end_index = 100

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

    return (intensity_values)

def visualize_abi(pos, subject, abi_data, query, abiname=None, zero_position=0, label=False, single=False):
    values = [] 
    ax  = pw.Brick(figsize=(pos.x1 - pos.x0, 0.6)) 
    axpos = ax.get_position() 
    
    ax2 = pw.Brick(figsize=(pos.x1 - pos.x0, 0.6))
    ax2.set_position([axpos.x0, axpos.y0, abs(axpos.x0-axpos.x1), abs(axpos.y0-axpos.y1)]) 
    
    ax.set_zorder(1) 
    ax2.set_zorder(0) 
    
    ax.plot(abi_data["_channel_wgap"]["G"][0], abi_data["_channel_wgap"]["G"][1], color="#e3e14f", lw=1, zorder=0) 
    ax.plot(abi_data["_channel_wgap"]["A"][0], abi_data["_channel_wgap"]["A"][1], color="#79E5B7", lw=1, zorder=0)
    ax.plot(abi_data["_channel_wgap"]["T"][0], abi_data["_channel_wgap"]["T"][1], color="#ff776c", lw=1, zorder=0) 
    ax.plot(abi_data["_channel_wgap"]["C"][0], abi_data["_channel_wgap"]["C"][1], color="#74b2d7", lw=1, zorder=0)  
    values.extend(abi_data["_channel_wgap"]["G"][1]) 
    values.extend(abi_data["_channel_wgap"]["A"][1]) 
    values.extend(abi_data["_channel_wgap"]["T"][1])
    values.extend(abi_data["_channel_wgap"]["C"][1])
    
    try:
        top    = abs(max(values)-min(values))*1.02
        bottom = -1.0*abs(max(values)-min(values))*0.02
    except:
        top    = 1000
        bottom = -20

    i = 0 
    for p, (q,s) in enumerate(zip(query, subject.seq)):
        if q != "/" and q != s:
            ax.bar([p+0.5], [top-bottom], bottom=bottom, width=1, lw=0.0, facecolor="#993366", edgecolor="#DDDDDD", zorder=1, alpha=0.2)

        if q != "-" and q != "/":
            ax2.bar([p+0.5], [abi_data["conf"][1][i]], bottom=0, width=1, lw=0.1, facecolor="#EFEFFF", edgecolor="#DDDDDD", zorder=1)
            i += 1

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)  
    ax.spines["right"].set_visible(False) 
    ax.spines["left"].set_visible(False) 
    ax.patch.set_alpha(0.0) 
    if abiname is None:
        pass
    else:
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(abiname, rotation=0, labelpad=10, ha="left")
    
    ax2.spines["top"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)  
    ax2.spines["right"].set_visible(False) 
    ax2.spines["top"].set_linewidth(0.4)  
    ax2.spines["bottom"].set_linewidth(0.4)  
    ax2.spines["right"].set_linewidth(0.4)  
    ax2.spines["left"].set_linewidth(0.4)  
    #ax2.patch.set_alpha(0.0) 

    ax.set_xlim(0, len(subject.seq)) 
    ax.set_xticks([]) 
    ax.set_ylim(bottom, top) 
    ax.set_yticks([])
    
    #ax2.plot([0,len(subject.seq)], [20,20], lw=0.2, ls="--", color="#BBBBBB")
    ax2.set_xlim(0, len(subject.seq))
    ax2.set_xticks([])
    ax2.set_ylim(0,70)
    ax2.set_yticks([0,30,60])
    if single == True:
        ax2.set_ylabel("Quality")
    axquery = pw.Brick(figsize=(pos.x1 - pos.x0, pos.y1 - pos.y0)) 
    colorbar(axquery, color_dict, query, subject.seq, char=True, fontsize=10, label=label, zero_position=zero_position)
    return ax, ax2, axquery 

def visualize(subject, abi_data, query, abiname=None, start=0, end=None, display_quality=True):
    if start == 0 and end is None:
        axmap  = visualizemap(subject, seq=True, linebreak=len(subject.seq)+1, tick_interval=10, title="", height_scale=0.8, fontsize_nucl=10, fontsize=10) 
    else: 
        axmap   = visualizemap(subject, seq=True, start=start, end=end, linebreak=end-start, tick_interval=10, title="", height_scale=0.8, fontsize_nucl=10, fontsize=10) 
        subject = cropdna(subject, start=start, end=end, quinable=0)
        if type(query) in (list, tuple):
            new_query_list = []
            abidata_list   = [] 
            for aquery, aabi_data in zip(query, abi_data):
                conf0 = []
                conf1 = [] 
                i = 0 
                for q in aquery[:start]:
                    if q != "-" and q != "/":
                        i += 1
                    
                j = 0 
                for q in aquery[start:end]:
                    if q != "-" and q != "/":
                        conf0.append(aabi_data["conf"][0][i+j])
                        conf1.append(aabi_data["conf"][1][i+j])
                        j += 1


                new_abidata = {"conf":[conf0, conf1], "_channel_wgap":{}}
                for nucl in "ATGC":
                    positions = [p-start for p in aabi_data["_channel_wgap"][nucl][0] if start <= p <= end]
                    values    = [v for p, v in zip(aabi_data["_channel_wgap"][nucl][0], aabi_data["_channel_wgap"][nucl][1]) if start <= p <= end]
                    new_abidata["_channel_wgap"][nucl] = [positions, values]
                new_query_list.append(aquery[start:end]) 
                abidata_list.append(new_abidata)
            query    = new_query_list
            abi_data = abidata_list
        
        else:
            conf0 = []
            conf1 = [] 
            i = 0 
            for q in query[:start]:
                if q != "-" and q != "/":
                    #conf0.append(abi_data["conf"][0][i+j])
                    #conf1.append(abi_data["conf"][1][i+j])
                    i += 1
                
            j = 0 
            for q in query[start:end]:
                if q != "-" and q != "/":
                    conf0.append(abi_data["conf"][0][i+j])
                    conf1.append(abi_data["conf"][1][i+j])
                    j += 1

            new_abidata = {"conf":[conf0, conf1], "_channel_wgap":{}}
            for nucl in "ATGC":
                positions = [p-start for p in abi_data["_channel_wgap"][nucl][0] if start <= p <= end]
                values    = [v for p, v in zip(abi_data["_channel_wgap"][nucl][0], abi_data["_channel_wgap"][nucl][1]) if start <= p <= end]
                new_abidata["_channel_wgap"][nucl] = [positions, values]
            abi_data = new_abidata
            query    = query[start:end] 

    keys   = list(axmap.bricks_dict.keys())
    axmap.bricks_dict[keys[1]].set_xticks([])
 
    pos     = axmap.bricks_dict[keys[1]].get_position() 
    patches = axmap.bricks_dict[keys[1]].containers[0].patches
    for s, patch in zip(subject.seq, patches):
        patch.set_alpha(0.6) 
        patch.set_facecolor(color_dict[s]) 
        #patch.set_facecolor("k")
 
    if type(query) in (list, tuple):
        pw.param["margin"] = None
        mappos = axmap.bricks_dict[keys[0]].get_position()
        axmap.bricks_dict[keys[0]].change_plotsize([mappos.x1-mappos.x0, (mappos.y1-mappos.y0)*0.95])
        axmap = axmap.bricks_dict[keys[0]]/axmap.bricks_dict[keys[1]]

        subax_list = []  
        for i, (aabi_data, aquery) in enumerate(zip(abi_data, query)):
            if i == len(query) - 1:
                ax1, ax2, axquery  = visualize_abi(pos, subject, aabi_data, aquery, abiname=abiname[i], zero_position=start, label=True) 
            else:
                ax1, ax2, axquery  = visualize_abi(pos, subject, aabi_data, aquery, abiname=abiname[i], zero_position=start, label=False)
            
            if display_quality == True:
                ax = pw.Bricks({ax1.get_label():ax1, ax2.get_label():ax2})
            else:
                ax1.patch.set_alpha(1.0) 
                ax = ax1
            
            subax = (ax/axquery) 
            subax_list.append(subax) 
        pw.param["margin"] = 0.05
        sub_axes = pw.stack(subax_list, operator="/")
        sub_axes.set_supylabel("Quality", labelpad=4) 
        sub_axes.set_supspine()
        sub_axes._case.spines["left"].set_linewidth(0.4) 
        pw.param["margin"] = None
        ax_all = axmap/sub_axes
        pw.param["margin"] = 0.05
    else:
        pw.param["margin"] = None
        mappos = axmap.bricks_dict[keys[0]].get_position()
        axmap.bricks_dict[keys[0]].change_plotsize([mappos.x1-mappos.x0, (mappos.y1-mappos.y0)*0.95])
        axmap = axmap.bricks_dict[keys[0]]/axmap.bricks_dict[keys[1]]
        ax1, ax2, axquery = visualize_abi(pos, subject, abi_data, query, abiname=abiname, zero_position=start, label=True, single=True) 
        
        if display_quality == True:
            ax = pw.Bricks({ax1.get_label():ax1, ax2.get_label():ax2})
        else:
            ax1.patch.set_alpha(1.0) 
            ax = ax1
        ax_all = axmap/(ax/axquery) 
    return ax_all



# TRIALS
filename = "my_chromatograms/MH3_F.ab1"

with open("my_chromatograms/MH3_F.txt", "r") as sequence_file:
    reference_seq = sequence_file.read().replace("\n", "").strip()

abi_data = abi_to_dict(filename)

consensus_seq, reverse_seq = generate_consensusseq(abi_data)

visualize(consensus_seq, filename, reverse_seq)
