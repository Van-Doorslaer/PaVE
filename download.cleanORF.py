#import required functions
from Bio import GenBank, SeqIO, Entrez
from Bio.Seq import Seq
from ORFtools import TransORF, Blast, grem, download, annotate_E2BS, wangcomputing, spliceview, split_uppercase, annotate_E8_E2, annotate_E4, test_start_stop, pairwise, presence, ref_table
import os, glob, re, urllib2, xlrd, csv, os.path
from mechanize import Browser

ref_table()



with open("new_viruses_duplicates_removed.txt", "rU") as f:
    for line in f:
        line=line.strip()
        handle=Entrez.efetch(db='nucleotide',id=line,rettype='gb')
        out.write(handle.read())
        handle.close()

for seq_record in SeqIO.parse("new_viruses.gb", "genbank"):
    print seq_record.id
    
    ##############################################
    ##annotate base ORFs based on BLAST###########
    #############################################
    not_new=[]
    to_be_deleted=[]
    for seq_record in SeqIO.parse("new_viruses.gb", "genbank"):
        stop=0
        record= seq_record.seq
        accession=seq_record.id.strip()
        table = 11                                                              #mammalian codon usage
        min_pro_len = 25                                                        #minimal size of the protein ()
        out = open('result.csv', 'w')                                          #open file to save results
        print >> out, accession+',CG,,,'+record
        orf_list = TransORF(record, table, min_pro_len)                     #translate virus on all three frames, calculate start and end position
        
        sizes=[]
        for protein_sequence in range(0,(len(orf_list)),3):                     #for all the identified ORFs do the following
            prot_seq = orf_list[protein_sequence]                               #the translated fragment flanked by two * is used from orf_list
            start = orf_list[protein_sequence+1]                                #the start is used from orf_list
            end = orf_list[protein_sequence+2]                                  #the end is used from orf_list
            output = Blast(accession, prot_seq, start, end, record)         #blast and annotate the ORFs starting with a M
            print >> out, ','.join(map(str, output))                            #print outfile to file
            if len(output)>1:
                sizes.append(output[2])
                if output[1] == "E1":
                    E1=output[4]
                    E1_start=output[2]
                elif output[1] == "E2":
                    E2=output[4]
                    E2_stop=int(output[3])
                    E2_start=int(output[2])
                elif output[1] == "L1":
                    L1_stop=int(output[3])
                    with open(accession+"_L1.fas", "w") as L1_fasta:
                        print >>L1_fasta, ">", accession+"_L1"
                        print >>L1_fasta, output[4]
                if int(output[2])<len(record):
                    small=output[2]
        CG=record
    
        if float(pairwise(accession, "all_L1.mafft.fas")) > 90:
            print accession +" is not a novel type"
            stop=1 ##CHANGE to stop=1
            not_new.append(accession)
    
        if stop != 1:
            ##############################################
            ##annotate E4 and E1^E4######################
            #############################################
            
            E4_list=annotate_E4(accession, E2, CG)
            E4=E4_list[4]
            E4_splice=None
            wang=0
            E4_start=E4_list[2]
            translation=str(Seq(E4).translate())
            r=re.search("M", translation)
            if r:
                E1_E4 = wangcomputing(accession, CG, E1, E4, E1_start, E4_start)
                if E1_E4 != []:
                    print "wangcomputing"
                    print >>out, ",".join(map(str,E1_E4))
                    #print E1_E4
                    E4_splice=int(E1_E4[2].split("+")[1].split("..")[0])
            
                else:
                    E1_E4 = spliceview(accession, CG, E1, E4, E1_start, E4_start)
                    print "spliceview"
                    if E1_E4 !=[]:
                        print >>out, ",".join(map(str,E1_E4))
                        E4_splice=int(E1_E4[2].split("+")[1].split("..")[0])
                        
            
                results=[]
                new_E4=str(E4)[r.start()*3:]
                #print E4
                new_E4_start= re.search(str(new_E4),str(CG)).start()
                if E4_splice != None:
                    if new_E4_start and new_E4_start<E4_splice:
                        results.append(accession)
                        results.append("E4")
                        results.append(str(new_E4_start+1))
                        results.append(str(new_E4_start+len(new_E4)))
                        results.append(new_E4)
                        print >>out, ",".join(results)
                else:
                    results.append(accession)
                    results.append("E4")
                    results.append(str(new_E4_start+1))
                    results.append(str(new_E4_start+len(new_E4)))
                    results.append(new_E4)
                    print >>out, ",".join(results)
                    
            else:
                print >>out, ",".join(map(str,E4_list))
                E1_E4 = wangcomputing(accession, CG, E1, E4, E1_start, E4_start)
                wang=1
                if E1_E4 != []:
                    #print "wangcomputing"
                    print >>out, ",".join(map(str,E1_E4))
                    E4_splice=int(E1_E4[2].split("+")[1].split("..")[0])
                    #print E1_E4
                else:
                    E1_E4 = spliceview(accession, CG, E1, E4, E1_start, E4_start)
                    #print "spliceview"
                    if E1_E4 !=[]:
                        print >>out, ",".join(map(str,E1_E4))
                        E4_splice=int(E1_E4[2].split("+")[1].split("..")[0])
    
            
            ##############################################
            ##annotate E8^E2 #############################
            ##############################################
            if E4_splice != None:
                E8_E2=annotate_E8_E2(accession, CG, E1, E2, E4_splice, E2_stop)
                print >>out, ",".join(map(str,E8_E2))
            
            ##############################################
            ##add E2BS####################################
            ##############################################
            
            E2BS = annotate_E2BS(CG)
            for x in E2BS:
                results=[]
                results.append(accession)
                results.append("E2BS")
                results.append(str(x+1))
                results.append(str(int(x)+12))
                results.append(CG[int(x):int(x)+12])
                print >>out, ",".join(map(str,results))
            
            ##############################################
            ##add URR#####################################
            ##############################################
            for s in sizes:
                if L1_stop<int(s)<len(CG):
                    URR_end=int(s)-1
                else:
                    URR_end=int(sorted(sizes, key=int)[0])-1
            URR_start= L1_stop+1
            if URR_end !=0:
                results=[]
                results.append(accession)
                results.append("URR")
                results.append("join("+str(URR_start)+".."+str(len(CG))+"+"+str(1)+".."+str(URR_end)+")")
                print >>out, ",".join(map(str,results))
            elif L1_stop == len(CG):
                results=[]
                results.append(accession)
                results.append("URR")
                results.append("1.."+str(URR_end)+")")
                print >>out, ",".join(map(str,results))
                
            else:
                results=[]
                results.append(accession)
                results.append("URR")
                results.append(str(URR_start)+".."+str(len(CG)))
                print >>out, ",".join(map(str,results))
            
        out.close()
        #remove temp files
        remove_BLAST = re.compile('.*BLAST')
        remove_TEMP = re.compile('.*temp')
        dir=os.curdir
        grem(dir,remove_BLAST)
        grem(dir,remove_TEMP)
        to_be_deleted.append(accession+"_L1.fas")
        
        
    ####################################################
    ##produce results.csv file##########################
    ####################################################
    
    result=open("result.csv","rU")   
    out=open("results.csv",'a')
    F_results=[]
    for line in result:
            line=line.rstrip().rsplit(",")
            if len(line)>1:
                if line[0] not in not_new:
                    if line[1] == "CG":
                        del line[3]
                        print >> out, ','.join(map(str, line))
                    elif "join" in line[2]:
                            print >>out, ",".join(line)
                    else:
                        line[2:4]=['..'.join(map(str,line[2:4]))]
                        print >> out, ','.join(map(str, line))
    
    out.close()
    result.close()
    to_be_deleted.append("result.csv")
    
    accessions=[]  
    with open("log.txt","w") as log:
        QC=test_start_stop("results.csv")
        if QC!=None:
            print >>log, QC
        for seq_record in SeqIO.parse("new_viruses.gb", "genbank"):
            accessions.append(seq_record.id)
            for accession in accessions:
                presence("output.fas", accession, log) #if doing more than a single record at a time, this needs to be edited to get a list of accession numbers and feed these into the function
    print 'all done'
    
for delete in to_be_deleted:
    if os.path.isfile(delete):
        os.remove(delete)
