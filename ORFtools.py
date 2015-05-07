from Bio import SeqIO, GenBank, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete2 import Tree
from collections import Counter
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from mechanize import Browser
import os, glob, re, csv, time, operator



def get_accession(query, database, rettype):
    """
    Returns a nucleotide sequence from genbank.

    :param query: the accession number
    :param database: the database to use (default=nucleotide)
    :param rettype: the return type for the sequence (native,fasta,gb,xml)

    :return: text of sequence in requested `rettype`
    :rtype: string

    """

    params = {
        'db': database,
        'tool': _toolname,
        'email': _email,
        'id': query,
        'rettype': rettype,
    }
    url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
    url = url + urllib.urlencode(params)
    data = urllib.urlopen(url).read()
    return data


def download(accession_number):
    accession=[]
    with open('file.temp','w') as temp:
	file = get_accession(accession_number, 'nucleotide','gb')
	print >> temp, file
    with open('file.temp','w') as temp:
	for record in SeqIO.parse("file.temp",'gb'):
	    print >> temp, ">"+record.description+"\n"+record.seq
    

def TransORF(seq, trans_table, min_protein_length):   #this function will translate the genome in all 6 frames and split the ORF on the stop codon *
    """translate all three frames into proteins. Next Isolate ORFs based on being flanked by *"""
    ORFs = []
    seq_len = len(seq)
    for frame in range(3):
        trans = str(seq[frame:].translate(trans_table))
        trans_len = len(trans)
        aa_start = 0
        aa_end = 0
        while aa_start < trans_len:
            aa_end = trans.find("*", aa_start)
            if aa_end == -1:   
               aa_end = trans_len
            if aa_end-aa_start >= min_protein_length:
                start = frame+aa_start*3
                end = min(seq_len,frame+aa_end*3+3)
                ORFs.append((trans[aa_start:aa_end]))
                ORFs.append(start)
                ORFs.append(end)
            aa_start = aa_end+1
    #ORFs.sort()
    return ORFs
    

def Blast(type,protein_sequence,start, end, genomic_sequence):
        result=[]
	ORF=[]
        M = re.search('M',protein_sequence)
        if M:
            query = protein_sequence[M.start():]
	    query=query+"*"
            temp = open("temp.ORF", "w")
            print >>temp, '>blasting'
            print >>temp, query
            temp.close()
            cline=blastp(query="'temp.ORF'", db = "./Blast/DB.blast.txt", evalue=0.01, outfmt=5, out=type +".BLAST")
            os.system(str(cline))
            blast_out=open(type+".BLAST")
            string=str(blast_out.read())
            DEF=re.search("<Hit_def>((.*))</Hit_def>",string)
	    if DEF:
		if DEF.group(1)=='L1':
		    real_start=start+M.start()+M.start()+M.start()
		    result.append(type)
		    result.append('L1')
		    L1_pre=genomic_sequence[(start+3*M.start()):int(end)]
		    splice='(C|T)(C|T)(A|C|G|T)(C|T)AG(A)TG'
		    spliced=re.search(splice,str(L1_pre))
		    if spliced:
			start_L1 = int(spliced.start())+6
			if start_L1 % 3 == 0:
			    if start_L1 > 600:
				L1_post=L1_pre
				result.append(str(start_L1)+".."+str(end))
				result.append(str(L1_post))
				result.append(Seq(str(L1_post)).translate())
			    else:
				L1_post=L1_pre[start_L1:]
				result.append(str(int(real_start+1)+int(start_L1))+".."+str(end))
				result.append(str(L1_post))
				result.append(Seq(str(L1_post)).translate())
			else:
			    L1_post=L1_pre
			    result.append(str(real_start+1)+".."+str(end))
			    result.append(str(L1_post))
			    result.append(Seq(str(L1_post)).translate())
		    else:
			L1_post=L1_pre
			result.append(str(real_start+1)+".."+str(end))
			result.append(str(L1_post))
			result.append(Seq(str(L1_post)).translate())
		    
		else:
		    real_start=start+M.start()+M.start()+M.start()
		    result.append(type)
		    result.append(DEF.group(1))
		    result.append(str(real_start+1)+".."+str(end))
		    result.append(genomic_sequence[int(real_start):int(end)])
		    result.append(query)

        return result

def split_uppercase(string):
    upper=[]
    lower=[]
    for i in string:
        if i.isupper():
            upper.append(i)
        else:
            str,lower.append(i)
    return "".join(upper), "".join(lower)

def annotate_E4 (type, DNA, genomic_sequence):
    result=[]
    trans=DNA[1:len(DNA)].translate()
    E4=max(trans.split("*"), key=len)
    E4_start=re.search(str(E4), str(trans)).start()
    E4_end=re.search(str(E4), str(trans)).end()
    result.append(type)
    result.append("E4")
    E4_nt=str(DNA[(E4_start*3)+1:((E4_end+1)*3)+1])
    E4_nt_start=re.search(E4_nt,str(genomic_sequence)).start()
    E4_nt_end=E4_nt_start+len(E4_nt)
    result.append(str(E4_nt_start+1)+".."+str(E4_nt_end))
    result.append(E4_nt)
    result.append(Seq(E4_nt).translate())
    return result

def annotate_spliced (accession, record, E1_start, E1_end, E4, E4_start, E4_end, E2, E2_end, E8s):
    results=[]
    with open("table","w") as output:
        browser = Browser()
        browser.open("http://wangcomputing.com/assp/index.html")
        browser.select_form(nr=0)
        browser['seqfield'] = record 
    
        ##retrieve results from website
        response = browser.submit()
        
        content = response.readlines()
        print >> output, content
    
    os.system("python html2csv.py table")
    donors = []
    acceptors = []
    for r in csv.reader(open("table.csv","rU") ):
        try:
            if "donor" in r[1]:
                SD = split_uppercase(r[2])[0]
                donors.append(SD)
            elif "acceptor in r[1]":
                SA = split_uppercase(r[2])[0]
                acceptors.append(SA.upper())
                
        except:
            continue
    
    d_end={}    
    for d in donors:
        start=re.search(d, str(record)).start()
        d_end[d]=start+len(d)
    
    E1_SD_dict={}
    E8_SD_dict={}
    for de in d_end:
        if de[-2:]=="AG" and  int(d_end[de]) in range(E1_start, E1_end+1):
            E1_SD_dict[ de ]=int(d_end[de])
        for E8 in E8s:
            E8_start = re.search(str(E8), str(record)).start()+1
            E8_end = E8_start+len(E8)
            if de[-2:]=="AG" and  int(d_end[de]) in range(E8_start, E8_end+1):
                E8_SD_dict[ de ]=int(d_end[de])
    
    sorted_E1_SD_dict = sorted(E1_SD_dict.items(), key=operator.itemgetter(1))
    sorted_E8_SD_dict = sorted(E8_SD_dict.items(), key=operator.itemgetter(1))
    
    
    a_start={}    
    for a in acceptors:
        try:
            start=re.search(a, str(record)).start()
            a_start[a]=start
        except:
            continue
    
    E4_SA_dict={}
    for sa in a_start:
        if str(record[a_start[sa]-2:a_start[sa]]) == "AG":
            if a_start[sa] in range(E4_start, E4_end+1):
                E4_SA_dict[sa] = int(a_start[sa])
    
    sorted_E4_SA_dict = sorted(E4_SA_dict.items(), key=operator.itemgetter(1))
    
    for SD in sorted_E1_SD_dict:
        for E4_SA in sorted_E4_SA_dict:
            E1_E4 = (record[E1_start-1:SD[1]]+record[E4_SA[1]:E4_end])
            coordinates=(E1_start, SD[1], E4_SA[1]+1,E4_end)
            if str(E1_E4.translate()[-20:]) == str(E4.translate()[-20:]):
                results.append([accession, "E1^E4","join(%i..%i+%i..%i)" %coordinates, E1_E4, E1_E4.translate()])
                break
        break
    try:
	for SD in sorted_E8_SD_dict:
	    E8_E2 = (record[E8_start-1:SD[1]]+record[E4_SA[1]:E2_end])
	    coordinates=(E8_start, SD[1], E4_SA[1]+1,E2_end)
	    if str(E8_E2.translate()[-20:]) == str(E2.translate()[-20:]):
		results.append([accession, "E8^E2","join(%i..%i+%i..%i)" %coordinates, E8_E2, E8_E2.translate()])
		break
    except:
	print "it is likely that %s has no E4 splice acceptor. Please check" % accession
    return results


            
def annotate_E8(accession, E1, record):
    y=0
    x=0
    count=[]
    sequences=[]
    for r in E1[1:].translate().split("*"):
        if x==0:
            y=y+len(r)+1
            count.append(y)
            if 100<count[-1]<250 and "M" in r:
                    sequence=E1[(count[-2]*3)+1+(re.search("M",str(r)).start())*3:((count[-2]*3)+1+(re.search("M",str(r)).start())*3)+len(r[re.search("M",str(r)).start():])*3].upper()
		    sequences.append(sequence)
    return sequences

def annotate_E2BS(CG):
    find="ACC......GGT"
    CG=str(CG)
    return [m.start() for m in re.finditer(find, CG+CG) if m.start() < len(CG)]

def results2dict(file):
    data={}
    accessions=[]
    
    
    for r in csv.reader(open(file,"rU")):
        if r[0] not in accessions:
            accessions.append(r[0])
    for a in accessions:
        genes={}
        E2BS_dict={}
        for r in csv.reader(open(file,"rU")):
            if a == r[0]:
                if r[1] != "E2BS":
                    genes[r[1]]=r[2:]
                else:
                    E2BS_dict [r[2]] =r[3]
                genes["E2BS"]=E2BS_dict		
        data[a]=genes
    return data


def test_start_stop(file):		
    accessions=[]
    with open(file, "rU") as f:
	for line in f:
	    line=line.strip().split(",")
	    if line[0] not in accessions:
		accessions.append(line[0])
    for accession in accessions:
	with open(file, "rU") as f:
	    for line in f:
		line=line.strip().split(",")
		if line[0]==accession:
		    if line[1]=="CG":
			CG=line[3]
		    elif "join" not in line[2] and line[1] != "URR":
			start=int(line[2].split("..")[0])-1
			end=int(line[2].split("..")[1])
			if CG[start:end] != line[3]:
			    return accession, line[1], line[2]
			        
    
def presence(alignment, accession, log):
    t=Tree("output.phy_phyml_tree.txt")
    A = t&accession
    distances={}
    for leaf in t:
        distances[(t.get_distance(A,leaf, topology_only=True))]=leaf
    neighbor= distances[sorted(distances)[1]]
    align = AlignIO.read("output.phy", "phylip")
    with open("neighbor.fas","w") as out:
        for record in align:
	    if accession in record.id:
		print >>out, ">" + accession
		print >>out, record.seq
	    if record.id in neighbor:
		print >>out, ">" + record.id
		print >>out, record.seq

    align = AlignIO.read("neighbor.fas", "fasta")
	
    A=list(align[0])
    B=list(align[1])
    count=0
    gaps=0
    for n in range(0, len(A)):
        if A[n]==B[n]:
            if A[n]!="-":
                count=count+1
            else:
                gaps=gaps+1
    identity = 100*(count/float((len(A)-gaps)))

    with open('./files/ref_table.csv', 'rU') as f:
        for line in f:
            line=line.strip().split(",")
            if "--"+line[0] in str(neighbor):
                species = line[1]

    expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E6","E7","L1","L2","CG","E2BS","URR"])
    if identity > 60.0:
	if species.split(" ")[0] == "Alphapapillomavirus":
            if identity > 70.0 and len(species.split(" "))>1:
                if species.split(" ")[1] in ["7","9","11"]:
                    expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E5_alpha","E6","E7","L1","L2","CG","E2BS","URR"])
                if species.split(" ")[1] in ["5","6"]:
                    expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E6","E7","L1","L2","CG","E2BS","URR"])
                if species.split(" ")[1] in ["2","3","4","14"]:
                    expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E5_beta","E6","E7","L1","L2","CG","E2BS","URR"])
                if species.split(" ")[1] == "10":
                    expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E5_gamma","E5_delta","E6","E7","L1","L2","CG","E2BS","URR"])
                if species.split(" ")[1] in ["8"]:
                    expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E5_delta","E6","E7","L1","L2","CG","E2BS","URR"])
                if species.split(" ")[1] == "12":
                    expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E5_epsilon","E5_zeta","E6","E7","L1","L2","CG","E2BS","URR"])
                if species.split(" ")[1] == "13":
                    expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E5","E6","E7","L1","L2","CG","E2BS","URR"])

        elif species.split(" ")[0] == "Gammapapillomavirus":
            if len(species.split(" ")) >1 and identity > 70.0:
                if species.split(" ")[1] == "6":
                    expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E7","L1","L2","CG","E2BS","URR"])
        elif species.split(" ")[0] in ["Dyopipapillomavirus","Dyodeltapapillomavirus","Omikronpapillomavirus","Upsilonpapillomavirus","Omegapapillomavirus"]:
            expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E6","L1","L2","CG","E2BS","URR"])
        elif species.split(" ")[0] in ["Kappapapillomavirus"]:
            expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E7","E6","L1","L2","E5","CG","E2BS","URR"])
        elif species.split(" ")[0] in ["Thetapapillomavirus"]:
            expected=Counter(["E1","E2","E7","E9","L1","L2","CG","E2BS","URR"])
        elif species.split(" ")[0] in ["Etapapillomavirus", "Dyoepsilonpapillomavirus"]:
            expected=Counter(["E1","E2","E6","E7","E9","L1","L2","CG","E2BS","URR"])
        elif species.split(" ")[0] == "Deltapapillomavirus":
            if identity > 70.0:
                if species.split(" ")[1] in ["1","2","3","4","5"]:
                    expected=Counter(["E1","E2","E4","E1^E4","E8^E2","E5","E6","E7","L1","L2","CG","E2BS","URR"])
        elif species.split(" ")[0] in ["Xipapillomavirus"]:
            expected=Counter(["E1","E2","E5 (E8)","E7","L1","L2", "E4","E1^E4","E8^E2","CG","E2BS","URR"])
	    
    found = []
    for line in csv.reader(open("./results/results.csv", "rU")):
            if accession == line[0]:
                found.append(line[1])
    if "E2BS" in found:
	found =filter(lambda a: a != "E2BS", found)
	found.append("E2BS")
    found = Counter(found)
    c=expected-found
    d=found-expected
    missing= list(c.elements())
    extra= list(d.elements())
     

    if extra==[] and missing !=[]:
        print >>log, "Based on %d percent sequence similarity with %s, %s was classified in %s. Based on this, the virus was expected to contain %s" %(identity, neighbor, accession, species, missing)
    elif missing==[] and extra !=[]:
        print >>log, "Based on %d percent sequence similarity with %s, %s was classified in %s. Based on this, the virus was not expected to contain %s" %(identity, neighbor, accession, species, extra)  
    elif missing!=[] and extra !=[]:
        print >>log, "based on %d percent sequence similarity with %s, %s was classified in %s. Based on this, the virus was not expected to contain %s. In addition, the virus contains an additional %s" %(identity, neighbor, accession, species, missing, extra)
    else:
	print "all is ok"
		    

    
    
