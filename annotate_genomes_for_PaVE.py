#import required functions
from Bio import GenBank, SeqIO, AlignIO, Entrez
from Bio.Seq import Seq
from ORFtools import TransORF, Blast, download, annotate_E2BS, split_uppercase, annotate_spliced, annotate_E4, annotate_E8, presence, results2dict, test_start_stop
import os, glob, re, csv
from mechanize import Browser
from os import listdir
from os.path import isfile, join

to_be_annotated=["NC_004068"]
genbank_file="./files/new_viruses.gb"

if not os.path.exists("results"):
    os.makedirs("results")


with open("./results/results.csv", "w") as out, open("./results/L1.fas", "w") as fasta:
    for to_be in to_be_annotated:
        ORFs={}
        Entrez.email = 'koenraad.vandoorslaer@nih.gov'
        with open(genbank_file,"w") as GB:
                    handle=Entrez.efetch(db='nucleotide',id=to_be,rettype='gb', retmod="text")
                    GB.write(handle.read())
                    handle.close()        
        ##############################################
        ##annotate base ORFs based on BLAST###########
        ##############################################
    
        for seq_record in SeqIO.parse(genbank_file, "genbank"):
            print seq_record.id
            record= seq_record.seq
            accession=seq_record.id.strip().split(".")[0]
            table = 11                                                              
            min_pro_len = 25
            ORFs={}
            print >> out, accession+',CG,,'+record
            orf_list = TransORF(record, table, min_pro_len)
            for protein_sequence in range(0,(len(orf_list)),3):
                prot_seq = orf_list[protein_sequence]
                start = orf_list[protein_sequence+1]
                end = orf_list[protein_sequence+2]                                  
                output = Blast(accession, prot_seq, start, end, record)         
                if len(output)>1:
                    print >> out, ','.join(map(str, output))
                    ORFs[output[1]]=output[2:]
            E1=ORFs["E1"][1]
            E2=ORFs["E2"][1]
            E1_start=int(ORFs["E1"][0].split("..")[0])
            E1_stop=int(ORFs["E1"][0].split("..")[1])
            E2_stop=int(ORFs["E2"][0].split("..")[1])
            L1_stop=int(ORFs["L1"][0].split("..")[1])
            L1=ORFs["L1"][1]
            print >>fasta, ">"+accession+"\n"+L1
    
            ##############################################
            ##annotate E4 and E1^E4######################
            #############################################
                    
            E4_list = annotate_E4 (accession, E2, record)
            ORFs["E4"]=E4_list[2:]
            E4 = Seq(ORFs["E4"][1])
            E4_start=int(E4_list[2].split("..")[0])
            E4_stop= int(E4_list[2].split("..")[1])            
            E8s = annotate_E8(accession, E1, record)
            
            splices = annotate_spliced (accession, record, E1_start, E1_stop, E4, E4_start, E4_stop, E2, E2_stop, E8s)
            for splice in splices:
                SA = int(splice[2].split("+")[1].split("..")[0])
                print >>out, ",".join(map(str,splice))
    
            r=re.search("M", str((E4).translate()))
            if r:          
                results=[]
                new_E4=str(E4)[r.start()*3:]
                new_E4_start= re.search(str(new_E4),str(record)).start()
                if splices != []:
                    if new_E4_start and new_E4_start<SA:
                        results.append(accession)
                        results.append("E4")
                        results.append(str(new_E4_start+1)+".."+str(new_E4_start+len(new_E4)))
                        results.append(new_E4)
                        results.append(Seq(new_E4).translate())
                        print >>out, ",".join(map(str,results))
                else:
                    results.append(accession)
                    results.append("E4")
                    results.append(str(new_E4_start+1)+".."+str(new_E4_start+len(new_E4)))
                    results.append(new_E4)
                    results.append(Seq(new_E4).translate())
                    print >>out, ",".join(map(str,results))
    
            
            ##############################################
            ##add E2BS####################################
            ##############################################
            
            E2BS = annotate_E2BS(record)
            for x in E2BS:
                results=[]
                results.append(accession)
                results.append("E2BS")
                if int(x)+12 > len(record):
                    E2BS_stop=int(x)+12 - len(record)
                    results.append("join("+str(x+1)+".."+str(len(record))+"+1.."+str(E2BS_stop))
                    results.append(record[int(x):int(x)+12]+record[:E2BS_stop])
                else:
                    results.append(str(x+1)+".."+str(int(x)+12))
                    results.append(record[int(x):int(x)+12])
                print >>out, ",".join(map(str,results))
            
            ##############################################
            ##add URR#####################################
            ##############################################
            sizes=[]
            for O in ORFs:
                sizes.append( ORFs[O][0].split("..")[0] )
            for s in sizes:
                if L1_stop<int(s)<len(record):
                    URR_end=int(s)-1
                    break
                else:
                    URR_end=int(sorted(sizes, key=int)[0])-1
            if L1_stop != len(record) and URR_end != 0:
                URR_start= L1_stop+1
                results=[]
                results.append(accession)
                results.append("URR")
                results.append("join("+str(URR_start)+".."+str(len(record))+"+"+str(1)+".."+str(URR_end)+")")
                print >>out, ",".join(map(str,results))
            elif L1_stop == len(record):
                results=[]
                results.append(accession)
                results.append("URR")
                results.append("1.."+str(URR_end))
                print >>out, ",".join(map(str,results))
            else:
                URR_start= L1_stop+1
                results=[]
                results.append(accession)
                results.append("URR")
                results.append(str(URR_start)+".."+str(len(record)))
                print >>out, ",".join(map(str,results))

data = results2dict("./results/results.csv") #this needs to be fixed to include accession



with open("./results/log.txt","w") as log:
    for d in data:
        for x in data[d]:
            if x != "E2BS":
                if len(data[d][x]) > 2:
                    if "join" not in data[d][x][0] and data[d][x][0] != "":
                        start = int(data[d][x][0].split("..")[0])-1
                        end = int(data[d][x][0].split("..")[1])
                        if data[d]["CG"][1][start:end] != data[d][x][1]:
                            print >>log, d, x, "the coordinates are wrong"
                        elif end not in range(0, len(data[d]["CG"][1])+1):
                            print >>log, d, x, "does not fall within the linear genome range"
                        else:
                            if str(Seq(data[d][x][1]).translate()) != str(data[d][x][2]):
                                print >>log, d, x, "the translation does not match"
                    else:
                        m = re.search( "\((\d*)\.\.(\d*)\+(\d*)\.\.(\d*)\)",str(data[d][x][0]))
                        if m:
                            start = int(m.groups()[0])-1
                            SD = int(m.groups()[1])
                            SA = int(m.groups()[2])-1
                            end = int(m.groups()[3])
                            if data[d]["CG"][1][start:SD]+data[d]["CG"][1][SA:end] != data[d][x][1]:
                                print >>log, d, x, "the coordinates are wrong"
                            elif end not in range(0, len(data[d]["CG"][1])+1):
                                print >>out, d, x, "does not fall within the linear genome range"
                            else:
                                if str(Seq(data[d][x][1]).translate()) != str(data[d][x][2]):
                                    print >>out, d, x, "the translation does not match"                                
            else:
                for e2bs in data[d][x]:
                    if int(e2bs.split("..")[1]) not in range(0, len(data[d]["CG"][1])+1):
                        print >>log, d, x, e2bs, "does not fall within the linear genome range"
        
    
    accessions=[]  
    QC=test_start_stop("./results/results.csv")
    if QC!=None:
        print >>log, QC
    os.system("mafft --add ./results/L1.fas --quiet --reorder ./files/all_L1.mafft.fas > output.fas")
    AlignIO.convert("output.fas", "fasta", "output.phy", "phylip")
    os.system("phyml output.phy 0 i 1 0 GTR 4.0 e 1 1.0 BIONJ n n")
    
    
    new_L1 = []
    for seq_record in SeqIO.parse("./results/L1.fas", "fasta"):
        new_L1.append(seq_record.id)
        for nL1 in new_L1:
            presence("output.fas", nL1, log) #if doing more than a single record at a time, this needs to be edited to get a list of accession numbers and feed these into the function
    print 'all done'




onlyfiles = [ f for f in listdir("./") if isfile(join("./",f)) ]

for o in onlyfiles:
    if ".py" not in o:
        os.remove(o)

