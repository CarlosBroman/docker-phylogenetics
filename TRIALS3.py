import os 
import sys 
import copy
import argparse
import patchworklib as pw
from QUEEN.queen import * # a framework to generate quinable and efficiently editable nucleotide sequence resources
from Bio import SeqIO
from Bio import pairwise2
import matplotlib
matplotlib.rcParams["xtick.major.width"] = 0.4
matplotlib.rcParams["xtick.minor.width"] = 0.4
matplotlib.rcParams["ytick.major.width"] = 0.4
matplotlib.rcParams["ytick.minor.width"] = 0.4

margin     = pw.param["margin"]
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

    return (consensus_seq, consensus_seq.translate(str.maketrans("ATGC","TACG"))[::-1]) 

def gl_alignment(template, query, single=False): 
    subject   = template.seq
    alignment = pairwise2.align.globalms(subject+subject, query, 2, 0, -10, -1, penalize_end_gaps=False, one_alignment_only=True)[0] 
    seqB  = alignment.seqB
    seqA  = alignment.seqA
    score = alignment.score
    for s, q in enumerate(seqB):
        if q == "-" or q == "/":
            pass 
        else:
            break 
    start = s
    
    for e, q in enumerate(seqB[::-1]):
        if q == "-" or q == "/":
            pass 
        else:
            break 
    
    end = len(2*subject) - e
    if start > len(subject):
        new_start = start - len(subject)
    else:
        new_start = start
    
    if end > len(subject): 
        new_end = end - len(subject) 
    else:
        new_end = end  

    if single == True:
        if new_start > new_end:
            nongap = joindna(cropdna(template, new_start, len(subject), quinable=0), cropdna(template, 0, new_end, quinable=0), quinable=0)
        else:
            nongap = cropdna(template, new_start, new_end) 
        
        op = 0
        flag = 0 
        gap_positions = []
        ops           = []
        for p, s in enumerate(seqA[start:-1*e]):
            if s == "-" or s == "/":
                if flag == 0:
                    gapstart    = p
                    gapstart_op = op
                    flag = 1
                else:
                    pass 

            else:
                if flag == 1:
                    gapend = p
                    gap_positions.append((gapstart, gapend))
                    ops.append(gapstart_op) 
                    flag = 0
                else:
                    pass 
                op += 1

        if len(gap_positions) == 0:
            wgap = nongap 
        else: 
            pre_op = 0
            for p, (op, gap_pos) in enumerate(zip(ops,gap_positions)):
                if p == 0:
                    wgap = joindna(cropdna(nongap, 0, op), QUEEN(seq="-"*(gap_pos[1]-gap_pos[0])))    
                else:
                    wgap = joindna(wgap, cropdna(nongap, pre_op, op), QUEEN(seq="-"*(gap_pos[1]-gap_pos[0])))    
                pre_op = op
            
            if pre_op - len(nongap.seq) == 0:
                pass
            else:
                wgap = joindna(wgap, cropdna(nongap, pre_op, len(nongap.seq)))
        
        return wgap, seqB[start:-1*e], nongap, new_start, new_end, alignment.score
    
    else:
        return seqB[start:-1*e], new_start, new_end, alignment.score    

def reform_abidata(abi_data, seqB, strand):
    # Prepare data for the reverse sequence.
    # Account for gaps in the alignment.
    if strand == -1:
        new_channel = copy.deepcopy(abi_data["_channel"]) 
        new_channel["G"][0] = abi_data["_channel"]["C"][0] 
        new_channel["G"][1] = abi_data["_channel"]["C"][1][::-1]
        new_channel["A"][0] = abi_data["_channel"]["T"][0] 
        new_channel["A"][1] = abi_data["_channel"]["T"][1][::-1]
        new_channel["T"][0] = abi_data["_channel"]["A"][0] 
        new_channel["T"][1] = abi_data["_channel"]["A"][1][::-1]
        new_channel["C"][0] = abi_data["_channel"]["G"][0] 
        new_channel["C"][1] = abi_data["_channel"]["G"][1][::-1]
        abi_data["_channel"] = new_channel
        
        new_channel = copy.deepcopy(abi_data["channel"]) 
        new_channel["G"][0] = abi_data["channel"]["C"][0] 
        new_channel["G"][1] = abi_data["channel"]["C"][1][::-1]
        new_channel["A"][0] = abi_data["channel"]["T"][0] 
        new_channel["A"][1] = abi_data["channel"]["T"][1][::-1]
        new_channel["T"][0] = abi_data["channel"]["A"][0] 
        new_channel["T"][1] = abi_data["channel"]["A"][1][::-1]
        new_channel["C"][0] = abi_data["channel"]["G"][0] 
        new_channel["C"][1] = abi_data["channel"]["G"][1][::-1]
        abi_data["channel"] = new_channel
        
        abi_data["conf"][0] = abi_data["conf"][0][::-1]
        abi_data["conf"][1] = abi_data["conf"][1][::-1]

    op   = 0 
    flag = 0 
    gap_positions = []
    ops           = []
    for p, s in enumerate(seqB):
        if s == "-" or s == "/":
            if flag == 0:
                gapstart    = p
                gapstart_op = op
                flag = 1
            else:
                pass 

        else:
            if flag == 1:
                gapend = p
                gap_positions.append((gapstart, gapend)) 
                ops.append(gapstart_op) 
                flag = 0
            else:
                pass 
            op += 1
    #print(seqB)  
    #rrint(gap_positions) 
    abi_data["_channel_wgap"] = {
                                 "A":[[],[]], 
                                 "T":[[],[]],
                                 "G":[[],[]],
                                 "C":[[],[]]
                                }
    opi    = 0 
    gapsum = 0 
    flag   = 0
    if len(ops) > 0:
        op     = ops[opi] * 10
        for i, v in enumerate(abi_data["_channel"]["G"][0]):
            #if i % 10 == 0 or i < 10:
            #    print(i, gapsum) 
            abi_data["_channel_wgap"]["G"][0].append(abi_data["_channel"]["G"][0][i] + gapsum) 
            abi_data["_channel_wgap"]["A"][0].append(abi_data["_channel"]["A"][0][i] + gapsum) 
            abi_data["_channel_wgap"]["T"][0].append(abi_data["_channel"]["T"][0][i] + gapsum)
            abi_data["_channel_wgap"]["C"][0].append(abi_data["_channel"]["C"][0][i] + gapsum) 
            
            abi_data["_channel_wgap"]["G"][1].append(abi_data["_channel"]["G"][1][i]) 
            abi_data["_channel_wgap"]["A"][1].append(abi_data["_channel"]["A"][1][i]) 
            abi_data["_channel_wgap"]["T"][1].append(abi_data["_channel"]["T"][1][i])
            abi_data["_channel_wgap"]["C"][1].append(abi_data["_channel"]["C"][1][i]) 

            if i >= op and flag == 0:
                gaplength = gap_positions[opi][1] - gap_positions[opi][0]
                for j in range(gaplength*10):
                    abi_data["_channel_wgap"]["G"][0].append(abi_data["_channel"]["G"][0][i] + gapsum + j*0.1) 
                    abi_data["_channel_wgap"]["A"][0].append(abi_data["_channel"]["A"][0][i] + gapsum + j*0.1) 
                    abi_data["_channel_wgap"]["T"][0].append(abi_data["_channel"]["T"][0][i] + gapsum + j*0.1)
                    abi_data["_channel_wgap"]["C"][0].append(abi_data["_channel"]["C"][0][i] + gapsum + j*0.1) 
                    
                    abi_data["_channel_wgap"]["G"][1].append(0) 
                    abi_data["_channel_wgap"]["A"][1].append(0) 
                    abi_data["_channel_wgap"]["T"][1].append(0)
                    abi_data["_channel_wgap"]["C"][1].append(0) 
                   
                gapsum += gaplength 
                opi    += 1
                if opi >= len(ops):
                    flag = 1
                    pass
                else:
                    op = ops[opi] * 10
    else:
        abi_data["_channel_wgap"] = copy.deepcopy(abi_data["_channel"]) 
    return abi_data         

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
    # Initialize map visualization
    axmap  = visualizemap(subject, seq=True, tick_interval=10, title="", height_scale=0.8, fontsize_nucl=10, fontsize=10) 
    keys   = list(axmap.bricks_dict.keys())
    axmap.bricks_dict[keys[1]].set_xticks([])

    # Visually represent the DNA
    pos     = axmap.bricks_dict[keys[1]].get_position() 
    patches = axmap.bricks_dict[keys[1]].containers[0].patches
    for s, patch in zip(subject.seq, patches):
        patch.set_alpha(0.6) 
        patch.set_facecolor(color_dict[s]) 
        #patch.set_facecolor("k")
 
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

def view_sanger(gbkpath, abipath, output=None, display_quality=True, output_cromatogram=None, output_fasta=None):
    # This function selects the file from the abipath, chooses the best matching sequence (forward or reverse) and calls the function visualize.
    # Also allows for choosing the outputs (figure / chromatogram data / fasta)
    template = QUEEN(record=gbkpath)
    project  = template.project    
    start = 0    
    end = len(template.seq) 

    abifile_path_list = [] 
    if os.path.isdir(abipath):
        for abifile_path in os.listdir(abipath):
            if ".ab1" == abifile_path[-4:]:
                abifile_path = abipath + "/" + abifile_path 
                abifile_path_list.append(abifile_path) 
        abifile_path_list.sort() 

    elif abipath.split(".")[-1] in ("abi", "ab1"):
        abifile_path_list.append(abipath)
    
    else:
        with open(abipath) as f:
            for line in f:
                line = line.rstrip()
                abifile_path_list.append(line) 

    abipath = abifile_path_list[0]
    abifile_name_list = [abipath.split("/")[-1]]
    abidata  = abi_to_dict(abipath)
    query    = generate_consensusseq(abidata)
    result_f = gl_alignment(template, query[0], single=True)
    result_r = gl_alignment(template, query[1], single=True)
    
    if result_f[-1] >= result_r[-1]:
        template_aligned, query_aligned, nongap, region_start, region_end, score = result_f
        strand = 1
    else:
        template_aligned, query_aligned, nongap, region_start, region_end, score = result_r
        strand = -1


    abidata = reform_abidata(abidata, query_aligned, strand)  

    ax_all  = visualize(template_aligned, abidata, query_aligned, display_quality=display_quality)
    
    if output is None:
        pass
    else:
        ax_all.savefig(output)


    if len(abifile_path_list) == 1:
        new_abidata_list   = [abidata]
        query_aligned_list = [query_aligned] 
        #template_aligned, new_abidata_list, query_aligned_list) 
    
    if output_cromatogram is None:
        pass 
    else:
        with open(output_cromatogram, "w") as f:
            for name, abidata in zip(abifile_name_list, new_abidata_list):
                name = name + ":" + str(start) + ".." + str(end) 
                row  = [name, "position"]
                row.extend(list(abidata["_channel_wgap"]["A"][0])) 
                print(*row, sep=",", file=f)  

                for nucl in "ATGC":
                    row = [name, "channel {}".format(nucl)] 
                    row.extend(list(abidata["_channel_wgap"][nucl][1])) 
                    print(*row, sep=",", file=f) 

    if output_fasta is None:
        pass 
    else: 
        with open(output_fasta, "w") as f:
            print(">{}".format(project), file=f)
            print(template_aligned.seq, file=f) 
            for name, query_seq in zip(abifile_name_list, query_aligned_list): 
                print(">{}".format(name), file=f)
                print(query_seq.replace("/","-"), file=f) 
            
    return ax_all 

# gbkpath = "sequence.gb"
# abipath = "MH3_F.ab1"

abipath = "Spec-2xU6gRNA-1.ab1"
gbkpath = "puc19_spec_2xu6grna.gb"

view_sanger(gbkpath, abipath, output="alignment.png", output_cromatogram="chromatogram.csv", output_fasta="sequences.fasta")