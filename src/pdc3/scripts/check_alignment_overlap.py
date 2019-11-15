import sys
import seq
import matplotlib.pyplot as plt
import networkx as nx
from networkx.drawing.nx_agraph import write_dot

def check_overlap(seq1,seq2):
    match = True
    once = False
    for i,j in zip(seq1,seq2):
        if i == "-" or j == "-":
            continue
        else:
            if i == j:
                once = True
                continue
            else:
                match = False
                break
    # need at least one match.
    if once == False:
        match = False
    return match

def combine_seqs(seq1,seq2):
    finalseq = ""
    for i,j in zip(seq1,seq2):
        if i == "-":
            finalseq += j
        else:
            finalseq += i
    return finalseq

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("python "+sys.argv[0]+" infile.aln outfile.aln")
        sys.exit()

    seqs = seq.read_fasta_file(sys.argv[1])
    seqs_sample = {}
    seqs_d = {}
    keep_seqs = {}
    for i in seqs:
        seqs_d[i.name] = i
        keep_seqs[i.name] = i
        spls = i.name.split("@")[0]
        try:
            seqs_sample[spls].append(i)
        except:
            seqs_sample[spls] = []
            seqs_sample[spls].append(i)

        

    for i in seqs_sample:
        G = nx.MultiGraph()
        for j in range(len(seqs_sample[i])):
            for k in range(len(seqs_sample[i])):
                if j < k:
                    ovr = check_overlap(seqs_sample[i][j].seq,seqs_sample[i][k].seq)
                    if ovr == True:
                        #print "\t",ovr,j,k,seqs_sample[i][j].name,seqs_sample[i][k].name
                        G.add_edge(seqs_sample[i][j].name,seqs_sample[i][k].name)
        if len(G) > 0:
            print(i +" "+str(len(seqs_sample[i])))
        for j in nx.connected_component_subgraphs(G):
            js = list(j)
            if len(j) == 2:
                #print "\tsimple 2 seq connect"
                #print "from"
                #print seqs_d[js[0]].seq
                #print seqs_d[js[1]].seq
                #print "to"
                print("\tcombining:",js,js[0]+"______"+js[1])
                finalse = combine_seqs(seqs_d[js[0]].seq,seqs_d[js[1]].seq)
                ns = seq.Sequence()
                ns.name = js[0]+"______"+js[1]
                ns.seq = finalse
                keep_seqs[ns.name] = ns
                print("\t\tremoving",js[0],js[1])
                del keep_seqs[js[0]]
                del keep_seqs[js[1]]
            else:
                going = True
                while going:
                    found = False
                    for k in j.edges():
                        found = True
                        ks = list(k)
                        print("\tcombining:",ks,ks[0]+"______"+ks[1])
                        finalse = combine_seqs(seqs_d[ks[0]].seq,seqs_d[ks[1]].seq)
                        ns = seq.Sequence()
                        ns.name = ks[0]+"______"+ks[1]
                        ns.seq = finalse
                        keep_seqs[ns.name] = ns
                        seqs_d[ns.name] = ns
                        x = set()
                        for m in j.neighbors(k[0]):
                            if m in k:
                                continue
                            x.add(m)
                        for m in j.neighbors(k[1]):
                            if m in k:
                                continue
                            x.add(m)
                        j.remove_node(k[0])
                        j.remove_node(k[1])
                        print("\t\tremoving",k[0],k[1])
                        del keep_seqs[k[0]]
                        del keep_seqs[k[1]]
                        for m in x:
                            if ns.name == m:
                                continue
                            if check_overlap(ns.seq,seqs_d[m].seq) == True:
                                j.add_edge(ns.name,m)
                                print("\t\t",ns.name,"->",m)
                        break
                    if found == False:
                        break
                    #for m in 
                    #j.remove_node(
                #nx.draw(G)
                #plt.show()
                #sys.exit()

    outf = open(sys.argv[2],"w")
    for i in keep_seqs:
        outf.write(keep_seqs[i].get_fasta())
    outf.close()
