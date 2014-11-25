from Bio import SeqIO, GenBank, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ete2 import Tree
from collections import Counter
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastpCommandline as blastp
from genbankdownload import get_accession
from mechanize import Browser
import os, glob, re, urllib2, xlrd, csv, time, os.path

to_be_deleted=[]

def download(accession_number):
    accession=[]
    temp = open('file.temp','w')
    file = get_accession(accession_number, 'nucleotide','gb')
    print >> temp, file
    temp.close()
    record = SeqIO.read("file.temp",'gb')
    temp = open('file.temp','w')
    print >> temp, ">"+record.description
    print >> temp, record.seq
    temp.close()
    

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
            temp = open("temp.ORF", "w")
            print >>temp, '>blasting'
            print >>temp, query
            temp.close()
            cline=blastp(query="'temp.ORF'", db = "DB.blast.txt", evalue=0.01, outfmt=5, out=type +".BLAST")
            os.system(str(cline))
            blast_out=open(type+".BLAST")
            string=str(blast_out.read())
            DEF=re.search("<Hit_def>((.*))</Hit_def>",string)
            #print blast_out.read()
	    if DEF:
		if DEF.group(1)=='L1':
		    real_start=start+M.start()+M.start()+M.start()
		    result.append(type)
		    result.append('L1')
		    L1_pre=genomic_sequence[(start+3*M.start()):int(end)]
		    splice='(C|T)(C|T)(A|C|G|T)(C|T)AG(A)TG'
		    spliced=re.search(splice,str(L1_pre))
		    #print L1_pre
		    if spliced:
			start_L1 = int(spliced.start())+6
			if start_L1 % 3 == 0:
			    if start_L1 > 200:		#changed this from 600 to 200; this may screw up annotations of L1!
				L1_post=L1_pre
				result.append(start_L1)
				result.append(end)
				result.append(str(L1_post))
				result.append(Seq(str(L1_post)).translate())
				
			    else:
				L1_post=L1_pre[L1_start:]
				result.append(int(real_start+1)+start_L1)  #I changed this, if L1 starts acting up, this is why
				result.append(end)
				result.append(str(L1_post))
				result.append(Seq(str(L1_post)).translate())
			else:
			    L1_post=L1_pre
			    result.append(real_start+1)
			    result.append(end)
			    result.append(str(L1_post))
			    result.append(Seq(str(L1_post)).translate())
		    else:
			L1_post=L1_pre
			result.append(real_start+1)
			result.append(end)
			result.append(str(L1_post))
			result.append(Seq(str(L1_post)).translate())
		    
		else:
		    real_start=start+M.start()+M.start()+M.start()
		    result.append(type)
		    result.append(DEF.group(1))
		    result.append(real_start+1)
		    result.append(end)
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
    result.append(str(E4_nt_start+1))
    result.append(str(E4_nt_end))
    result.append(E4_nt)
    return result

def wangcomputing(accession, CG, E1, E4, E1_start, E4_start):
    results=[]
    results.append(accession)
    out=open("table","w")
    browser = Browser()
    browser.open("http://wangcomputing.com/assp/index.html")
    browser.select_form(nr=0)
    browser['seqfield'] = CG
    
    ##retrieve results from website
    response = browser.submit()
    
    content = response.readlines()
    print >> out, content
    os.system("python html2csv.py table")
    

    with open("table.csv", "rU") as f:
	acceptor={}
	donor={}
	x=0
	for line in f:
	    line=line.strip().split(",")
	    if re.search("\d", line[0]):
		if len(line)>3:
		    if "donor" in line[2]:
			donor[line[0].split("\\")[0].split("\"")[1]]=line[3].split("\"")[1]
		    if "acceptor" in line[2]:
			    if "donor" not in line[3]:
				acceptor[line[0].split("\\")[0].split("\"")[1]]=line[3].split("\"")[1]
    
    keylist = donor.keys()
    [int(k) for k in keylist]
    keylist.sort(key=int)
    for E1_key in keylist:
	if split_uppercase(donor[E1_key])[0] in E1.upper():
	    if E1_key > E1_start:
		E4_= CG[int(E1_start)-1:int(E1_key)]
		break
	    else:
		x=1
		results=[] 
		break 
    if x == 0:
	keylist = acceptor.keys()
	[int(k) for k in keylist]
	keylist.sort(key=int)
	for E4_key in keylist:
	    if split_uppercase(acceptor[E4_key])[0] in E4.upper():
		if int(E4_start)<int(E4_key)<int(E4_start)+len(E4)-1:
		    _E4= CG[int(E4_key)-1:int(E4_start)+len(E4)-1]
		    if str(Seq(E4).translate()[-20:])==str((E4_+_E4).translate()[-20:]):
			break
	if int(E1_key)>int(E1_start) and int(E4_key)<int(E4_start)+len(E4)-1:
	    results.append("E1^E4")
	    results.append("join("+str(E1_start)+".."+str(E1_key)+"+"+str(E4_key)+".."+str(int(E4_start)+len(E4)-1)+")")
	    results.append(E4_+_E4)
	else:
	    results=[]
    return results
    os.remove(table)
    os.remove(table.csv)

def spliceview(accession, CG, E1, E4, E1_start, E4_start):   
    results=[]
    results.append(accession)
    #y=0
    with open("table","w") as out:
        browser = Browser()
        browser.open("http://zeus2.itb.cnr.it/~webgene/wwwspliceview.html")
        
        browser.select_form(nr=0)
        browser['sequence'] = CG
        
        ##retrieve results from website
        response = browser.submit()
        
        content = response.readlines()
        content = [i.split('\n')[0] for i in content]
        donor_index=content.index("POSITION        EXON INTRON    SCORE")
        acceptor_index=content.index("POSITION        INTRON EXON    SCORE")
        end_index=content.index("</TEXTAREA></FORM>")
        get_donors=content[donor_index+1:acceptor_index-2]
        get_acceptors=content[acceptor_index+1:end_index-2]
        donors=[]
        acceptors=[]
        for gd in get_donors:
            for g in gd.split(" "):
                if g!="":
                    donors.append(g)
        #print donors
        for d in range(0, len(donors), 4):
            if int(E1_start)<int(donors[d])<int(E1_start)+len(E1):
                #print donors[d]
                #splice_donor=donors[d+1]
                E4_=CG[int(E1_start)-1:int(donors[d])]
                #print donors[d], E4_
                break

        for ga in get_acceptors:
            for g in ga.split(" "):
                if g!="":
                    acceptors.append(g)
        for a in range(0, len(acceptors), 4):
            #print acceptors[a]
            if int(E4_start)<int(acceptors[a])<int(E4_start)+len(E4):
                _E4=CG[int(acceptors[a])-1:int(E4_start)+len(E4)-1]
                #print acceptors[a], _E4
                #print Seq(E4_+_E4).translate()[-20:]
                #print Seq(E4).translate()[-20:]
                if str(Seq(E4_+_E4).translate()[-20:])==str(Seq(E4).translate()[-20:]):
                       y=0
                       break
            else:
                y=1
        #print y
        if y==0:
            results.append("E1^E4")
            results.append("join("+str(E1_start)+".."+str(donors[d])+"+"+str(acceptors[a])+".."+str(int(E4_start)+len(E4)-1)+")")
            results.append(E4_+_E4)
            if str(acceptors[a])>str(int(E4_start)+len(E4)-1):
                results=[]
        else:
            results=[]
            
        return results
            
def annotate_E8_E2(accession, CG, E1, E2, E4_splice, E2_stop):
    E8s={}
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
    donor=[]
    output=open("table","w")
    browser = Browser()
    browser.open("http://wangcomputing.com/assp/index.html")
    browser.select_form(nr=0)
    browser['seqfield'] = E1
    
    ##retrieve results from website
    response = browser.submit()
    
    content = response.readlines()
    print >> output, content
    #print content
    os.system("python html2csv.py table")
    with open("table.csv","rU") as f:
        for line in f:
            line=line.strip().split(",")
            if len(line)>3:
                if "donor" in line[2]:
                    donor.append(split_uppercase(line[3])[0])
    t=0
    for d in donor:
        for seq in sequences:
            if d in seq:
                results=[]
                splice = re.search(d, str(seq)).start()+len(d)
                E8=seq[:splice]
		#print E8
                E8_start=re.search(str(E8),str(CG)).start()+1
                E8_splice=E8_start+len(E8)
                coordinates="join("+str(E8_start)+".."+str(E8_splice)+"+"+str(E4_splice)+".."+str(E2_stop)+")"
                results.append(accession)
                results.append("E8^E2")
                results.append(coordinates)
                results.append(str((E8+CG[E4_splice-1:E2_stop])))
                if str((E8+CG[E4_splice-1:E2_stop]).translate()[-20:])==str(E2.translate()[-20:]):
                    if t==0:
                        return results
                        t=x+1
                        break
    os.remove("table")
    os.remove("table.csv")

def annotate_E2BS(CG):
    find="ACC......GGT"
    CG=str(CG)
    return [m.start() for m in re.finditer(find, CG+CG) if m.start() < len(CG)]

def grem(path, pattern):
    pattern=re.compile(pattern)
    for each in os.listdir(path):
        if pattern.search(each):
            name=os.path.join(path,each)
            try: os.remove(name)
            except:
                grem(name,'')
                os.rmdir(name)
		
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
			    
def pairwise(accession, existing_alignment):
    results=[]
    os.system("mafft --add "+accession+"_L1.fas --quiet --reorder " + existing_alignment + " > output.fas")
    IDs=[]
    alignment = AlignIO.read("output.fas", "fasta")
    for record in alignment:
	if record.id not in IDs:
	    IDs.append(record.id)
	    
	    
               
    alignment = AlignIO.read("output.fas", "fasta")
    for record in alignment:
	if accession+"_L1" == record.id:
	    sequence1=list(str(record.seq))
    for record in alignment:
	y=0
	z=0
	if accession+"_L1" != record.id:
	    sequence2=list(str(record.seq))
	    for x in range(0,len(sequence1)):
		if sequence1[x]==sequence2[x]:
		    if sequence1[x]!="-":
			y=y+1
		    else:
			z=z+1
	    #results.append(accession+"_L1")
	    #results.append(record.id)
	    results.append(str(100*(float(y)/float(len(sequence1)-z))))
    return sorted(results, key=float, reverse=True)[0]
    
def ref_table():   
    url = 'http://pave.niaid.nih.gov/prototype_data/Human_reference_clones.xls'
    f = urllib2.urlopen(url)
    data = f.read()
    with open("human_table.xls", "wb") as code:
	code.write(data)
    
    wb = xlrd.open_workbook('human_table.xls')
    sh = wb.sheet_by_name('Human_reference_clones.csv')
    with open('human_table.csv', 'w') as out:
	for rownum in xrange(sh.nrows):
	    entry=sh.row_values(rownum)
	    new_line=[]
	    for e in entry:
		new_line.append(e.encode('utf8', 'replace'))
	    print >>out, ",".join(new_line)
    
    
    url = 'http://pave.niaid.nih.gov/prototype_data/Animal_reference_clones.xls'
    f = urllib2.urlopen(url)
    data = f.read()
    with open("animal_table.xls", "w") as code:
	code.write(data)
    
    wb = xlrd.open_workbook('animal_table.xls')
    sh = wb.sheet_by_name('Animal_reference_clones.csv')

    with open('animal_table.csv', 'w') as out:
	for rownum in xrange(sh.nrows):
	    entry=sh.row_values(rownum)
	    new_line=[]
	    for e in entry:
		new_line.append(e.encode('utf8', 'replace'))
	    print >>out, ",".join(new_line)

    
    
    with open("ref_table.csv","a") as out:
	with open('human_table.csv', 'rU') as f:
	    for line in f:
		#print line
		if "papillomavirus"  in line:
		    line=line.rstrip().split(",")
		    print >>out, line[0]+","+line[1]
		    #print line
	with open('animal_table.csv', 'rU') as f:
	    for line in f:
		if "papillomavirus"  in line:
		    line=line.rstrip().split(",")
		    print >>out, line[0]+","+line[2]
    
    os.remove("human_table.csv")
    os.remove("human_table.xls")
    os.remove("animal_table.csv")
    os.remove("animal_table.xls")
    
def presence(alignment, accession, log):
    AlignIO.convert(alignment, "fasta", "output.phy", "phylip")
    os.system("phyml output.phy 0 i 1 0 GTR 4.0 e 1 1.0 BIONJ n n")
    t=Tree("output.phy_phyml_tree.txt")
    #accession=accession.split(".")[0]
    #print accession+"_L1"
    #print t&accession+"_L1"
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
    print align
	
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

    with open('ref_table.csv', 'rU') as f:
        for line in f:
            line=line.strip().split(",")
	    #print str(neighbor), "--"+line[0]
            if "--"+line[0] in str(neighbor):
                species = line[1]
		print species

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
    with open("results.csv", "rU") as f:
        for line in f:
            line=line.rstrip().rsplit(",")
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
        #print identity
        print >>log, "Based on %d percent sequence similarity with %s, %s was classified in %s. Based on this, the virus was expected to contain %s" %(identity, neighbor, accession, species, missing)
    elif missing==[] and extra !=[]:
        print >>log, "Based on %d percent sequence similarity with %s, %s was classified in %s. Based on this, the virus was not expected to contain %s" %(identity, neighbor, accession, species, extra)  
    elif missing!=[] and extra !=[]:
        print >>log, "based on %d percent sequence similarity with %s, %s was classified in %s. Based on this, the virus was not expected to contain %s. In addition, the virus contains an additional %s" %(identity, neighbor, accession, species, missing, extra)
    else:
	print "all is ok"
		    
to_be_deleted.append("output.phy_phyml_tree.txt")
to_be_deleted.append("output.phy_phyml_stat.txt")
to_be_deleted.append("output.phy_phyml_lk.txt")
to_be_deleted.append("output.phy")
to_be_deleted.append("output.fas")
to_be_deleted.append("neighbor.fas")
to_be_deleted.append("ref_table.csv")
to_be_deleted.append("table")
to_be_deleted.append("table.csv")
for delete in to_be_deleted:
    if os.path.isfile(delete):
	os.remove(delete)
    
    
