import warnings
import shutil
import time
import pandas as pd
import subprocess
import os
import sys
import shutil
import re
subprocess.call("pip3 install pandas", shell=True)

warnings.filterwarnings("ignore")
analy_dir = sys.path[0]

print("-------------------------------------------------------------------")
print("# This is protein sequences analysis workflow")
print("# Please fellow the instructions to specify protein family and")
print("# taxonomic group")
print("# The pipeline would automatically execute the following analysis:")
print("# 1.Multiple sequence alignment")
print("# 2.Sequence conservation of a sequence alignment visualization")
print("# 3.Motif scaning")
print("# 4.Alternative analysis")
print("-------------------------------------------------------------------")

print("\n\n\n-----------------------The analysis will start!--------------------")
# step0: download edirect
subprocess.call(
    "sh -c \"$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)\"", shell=True)
subprocess.call(
    "echo \"export PATH=\$PATH:\$HOME/edirect\" >> $HOME/.bash_profile", shell=True)
subprocess.call("export PATH=${PATH}:${HOME}/edirect", shell=True)

# define creating folder function


def mkdir(folder_name):
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)
    else:
        shutil.rmtree(folder_name)
        os.mkdir(folder_name)


mkdir('result')
# step1: download sequence
print("# Step1: Download protein sequence...")

# This function is used to download protein sequences
# parameter t is taxonomic group; parameter p is protein family
def seq_download(t, p, r):
    if r == '1':
        count = int(subprocess.check_output("esearch -db protein -query" + " " + "\"" + p + " AND" + " " + t +
                    "[Organism] NOT PARTIAL" + "\"" + " | grep '<Count>' | cut -d '>' -f 2 | cut -d '<' -f 1", shell=True))
    if r == '2':
        count = int(subprocess.check_output("esearch -db protein -query" + " " + "\"" + p + " AND" + " " + t +
                    "[Organism] NOT PREDICTED" + "\"" + " | grep '<Count>' | cut -d '>' -f 2 | cut -d '<' -f 1", shell=True))
    if r == '3':
        count = int(subprocess.check_output("esearch -db protein -query" + " " + "\"" + p + " AND" + " " +
                    t + "[Organism]" + "\"" + " | grep '<Count>' | cut -d '>' -f 2 | cut -d '<' -f 1", shell=True))
    if r != '1' and r != '2' and r != '3':
        r = '3'  # if the restriction input is not in any case we expected, then assumes that the user has chosen no restrictions
    if count == 0:  # if esearch programme get no found, then stop the analysis
        print('The following term was not found in Protein:' +
              p + "[All Fields]" + " " + t + "[Organism]")
        print('Please check your query input.')
        print('Exit process...')
        exit()
    else:
        print("\n\t**********There are " +
              str(count) + " sequences.**********\n")
        if input("Would you like to download and start analyze them?(Enter y to continue)") == 'y':
            print("\t------------------------DOWNLOADING----------------------")
            if r == '1':
                subprocess.call("esearch -db protein -query" + " " + "\"" + p + " AND" + " " + t +
                                "[Organism] NOT PARTIAL" + "\"" + " | efetch -format fasta > protein.fasta", shell=True)
            if r == '2':
                subprocess.call("esearch -db protein -query" + " " + "\"" + p + " AND" + " " + t +
                                "[Organism] NOT PREDICTED" + "\"" + " | efetch -format fasta > protein.fasta", shell=True)
            if r == '3':
                subprocess.call("esearch -db protein -query" + " " + "\"" + p + " AND" + " " +
                                t + "[Organism]" + "\"" + " | efetch -format fasta > protein.fasta", shell=True)
            print("\t---------------------SEQUENCES-DOWNLOADED----------------")
        else:
            print("Existing...")
            exit()
    return count


# ask for the searching setting
taxonomy = input(
    'Please input query taxonomic group. (For example, you can try entering \'Aves\'):')  # 'Aves'
# 'Glucose-6-phosphatase'
protein = input(
    'Please input query protein family. (For example, you can try entering \'Glucose-6-phosphatase\'):')
restriction = input('What restrictions do you want to impose?\n\
1.NOT PARTIAL\n2.NOT PREDICTED\n3.No more restriction\t')
# call the defined funcion seq_download to download protein sequence for the following analysis
seq_total_number = seq_download(taxonomy, protein, restriction)


# This function is used to check how many species in query taxonomy
# parameter sequencefile is protein sequence file name
def species_count(sequencefile,seq_count):
    protein_fasta = open(sequencefile)
    species = []  # creating a list to save species name
    for i in protein_fasta.readlines():
        if i.startswith('>'):
            if re.findall(r'\[(.*?)\]', i):
                # match [SPECIES] in annotation line of each sequence
                species.append(re.findall(r'\[(.*?)\]', i)[0])
    speciesdic = {}
    for key in species:
        speciesdic[key] = speciesdic.get(key, 0) + 1
    species = set(species)  # return member of species list
    s = []
    print("There are", len(species), "species in taxonomic group", taxonomy)
    print("Protein sequences is saved in file \'protein.fasta\'")
    if seq_count > 1000:
        new_count = 0 # sequence count of new fasta file
        n = int(len(species)/2) 
        for i in range(0,n):
            # sort species dictionary by value
            s.append(sorted(speciesdic.items(),key = lambda x:x[1],reverse = True)[i][0])
        if os.path.exists('species_select'):
            os.remove('species_select')
        with open(sequencefile,'r') as f:
            for j in f.readlines():
                if j.startswith('>'):
                    if j.split('[')[1].split(']')[0] in s:
                        with open('species_select','a') as f:
                            f.write(j.split('>')[1].split()[0]) # write accession of selected species
                            f.write('\n')
                            new_count += 1 
            subprocess.call('/localdisk/data/BPSM/Assignment2/pullseq -i protein.fasta -n species_select > protein_sub.fasta',shell=True)
            os.remove('protein.fasta')
            os.rename('protein_sub.fasta','protein.fasta')
        print("The clustalo will align sequences of "+ str(n) +" most most frequent species.")
        print("There are total "+ str(new_count) +" species.")
        os.remove('species_select')
        return new_count
    else:
        print("The clustalo will align all downloaded sequences.")
        return seq_count

# call the defined funcion species_count to calculate species amount
final_seq_count = species_count('protein.fasta',seq_total_number)
print("------------------------------Step1 Finish--------------------------\n\n")

# step2: clustalo multiple sequence alignment
# multiple sequence alignment
print("------------------------------Step2 Start---------------------------")
print("# Step2: Multiple Sequence Alignment...\n")
print("\tStep2.1: Clustalo is running...\n")
subprocess.call(
    'clustalo -i protein.fasta -o protein_msa.fasta --threads 10', shell=True)
subprocess.call('clustalo -i  protein.fasta -o protein.msf --outfmt msf --threads 10',
                shell=True)  # use parameter --outfmt to change output format

mkdir('result/msa')

# This function is used to present multiple sequences alignment result
# taking msf file and mode(input from 'msastyle') as parameters


def msa_result(msf, program):
    if program == '1':  # if user input 1, then run showalign program
        print("\n\t\tEMBOSS \'showalign\' program is running to present multiples sequence alignment....")
        subprocess.call('showalign -sequence ' + msf +
                        ' -outfile ./result/msa/showalign.txt &> log', shell=True)
        print('\t\t====EMBOSS_showalign program finish==============')
        print('\t\t====The output is saved in file \'showalign\'====')
    elif program == '2':  # if user input 2, then run prettyplot program
        print("\n\t\tEMBOSS \'prettyplot\' program is runing to present multiple sequences alignment....")
        subprocess.call('prettyplot -sequences ' + msf +
                        ' -ratio=0.59 -docolour -graph svg &> log', shell=True)
        print('\t\t====EMBOSS_prettyplot program finish...')
        print('\t\t====File \'prettyplot.svg\' file is a visualization of multiple sequences alignment.')
    else:  # if user input anything other than 1 and 2, then run both programs
        print("\n\t\tEMBOSS \'showalign\' program is running to present multiple sequences alignment....")
        subprocess.call('showalign -sequence ' + msf +
                        ' -outfile ./result/msa/showalign.txt &> log', shell=True)
        print('\t\t====EMBOSS_showalign program finish')
        print('\t\t====The output is saved in file \'showalign\'')
        print("\n\t\tEMBOSS \'prettyplot\' program is runing to present multiple sequences alignment....")
        subprocess.call('prettyplot -sequences ' + msf +
                        ' -ratio=0.59 -docolour -graph svg &> log', shell=True)
        print('\t\t====EMBOSS_prettyplot program finish...')
        print('\t\t====File \'prettyplot.svg\' file is a visualization of multiple sequences alignment.')
    # run infoalign program and display head of the result
    print("\n\t\tEMBOSS \'infoalign\' is running to display basic information about multiple sequences alignment")
    subprocess.call('infoalign ' + msf +
                    ' -out ./result/msa/infoalign -nousa &> log', shell=True)
    print('\t\t====EMBOSS_infoalign program finish.\n')
    print('\t\tHere list the frist 5 rows of infoalign result:')
    info = pd.read_csv('./result/msa/infoalign', header=0, sep='\t')
    print(info.iloc[:, 0:5].head())
    print('\n\t\t====The output is saved in file \'infoalign\'')


# present multiple sequences alignment result
print("\tStep2.2: Display a multiple sequence alignment in pretty format...\n")
# define msa result present style via user input
msastyle = str(input("\tPlease select multiple sequence alignment presenting programme (please enter 1 or 2 or 3):\n\
    \t\t1.showalign program\n\t\t2.prettyplot program\n\
    \t\t3.I have no idea, please perform both.\n\
    \t\tPlease select one of the mode above:\n\
    \t\t*If you enter anything else, the pipeline will execute both for you."))
# call the defined funcion msa_result to present multiple sequences alignment result
msa_result('protein.msf', msastyle)
print("------------------------------Step2 Finish--------------------------\n\n")

# step3: Conservation Visualization via plotcon
# step3.1
# This funciton is used to perform blastp
# taking multiple sequences(in fasta format) and aligned sequences as input
# perform blastp program: generated consensus sequence is query sequence and database is produced from multiple sequences


def blast_process(pro_fasta, pro_msa_fasta):
    print("\t\tStep3.1.1: Generating consensus sequence...\n")
    subprocess.call('cons -sequence ' + pro_msa_fasta +
                    ' -outseq protein_msa_consensus.fasta &> log', shell=True)
    print("\t\tDone!\n")
    print("\t\tStep3.1.2: Making blast database...\n")
    subprocess.call('makeblastdb -in ' + pro_fasta +
                    ' -dbtype prot -out protein_db', shell=True)
    print("\t\tDone!\n")
    print("\t\tStep3.1.3: Performing blastp...\n")
    subprocess.call(
        'blastp -db protein_db -query protein_msa_consensus.fasta -outfmt 7 -max_hsps 1 > blast.out', shell=True)
    print("\t\tDone!\n")


print("------------------------------Step3 Start---------------------------")
print("# Step3: Conservation Visualization ...\n")
# call the defined funcion blast_process to present multiple sequences alignment result
print("\tStep3.1: Calculation sequence similarity...\n")
blast_process('protein.fasta', 'protein_msa.fasta')

# step3.2
# This function is used pull protein sequence
# take blastp output file and the number of sequence to pull as parameters
def pull_pro_seq(blastout, number):
    # perfrom analysis when the input number within (3,seq_total_number)
    if number > 3 and number < seq_total_number:
        hsp = pd.read_csv(blastout, comment='#', sep='\t', header=None)
        # sort the blastp output dataframe by bit score
        hsp.sort_values(by=[11], ascending=False)[1][0:number].to_csv(
            'header_name', index=False, header=None)
        subprocess.call(
            '/localdisk/data/BPSM/Assignment2/pullseq -i protein_msa.fasta -n header_name > extract_msa.fasta', shell=True)
        subprocess.call(
            '/localdisk/data/BPSM/Assignment2/pullseq -i protein.fasta -n header_name > ./result/msa/extract.fasta', shell=True)
        print('\tSequences to be analyzed are saved in file \'extract_msa.fasta\'')
        os.remove('header_name')
    else:
        exit()


print("\tStep3.2: Pull Sequences...\n")
print('\t', 'There are', str(final_seq_count), 'sequence in your fasta file.')
seq_number = input('\tPlease enter the number of sequences you would like to perform their conservation\
    (Please input an intege from 3 to the total number of your download sequence):')
if seq_number.isdigit():
    # call the defined funcion pull_pro_seq to pull sequences
    pull_pro_seq('blast.out', int(seq_number))
else:
    print('Please enter an integer!')
    exit()
print('\tDone!\n')


# step3.3: plotcon
print("\tStep3.3: Plot the level of conservation between the protein sequences.\n")
win_size = input(
    "\tPlease specify the window size to perform plotcon program(must be an integer):")
if win_size.isdigit():
    subprocess.call(
        'clustalo -i ./result/msa/extract.fasta -o ./result/msa/extract_msa.fasta --threads 10', shell=True)
    subprocess.call('plotcon -sequence ./result/msa/extract_msa.fasta -winsize ' +
                    win_size + ' -graph svg', shell=True)
    subprocess.call('mv plotcon.svg ./result', shell=True)
else:
    print('Please enter an integer!')
    exit()
print("Step3 Finish. Please check the file \'plotcon.svg\'")
print("------------------------------Step3 Finish--------------------------\n\n")

# step4: motif scaning
print("------------------------------Step4 Start---------------------------")
print("# Step4: Motif Scaning...\n")

# This function is used to split fasta file
# take multiple sequences file as input
def fasta_split(pro_fasta):
    mkdir('fasta')
    # extract each protein's sequence in a dictionary
    sequence = {}
    with open(pro_fasta, 'r') as f:
        for lines in f:
            if '>' in lines:
                protein_accesion = lines.replace('>', '').split(' ')[0]
                temp_seqname = open("seqname", "w")
                temp_seqname.write(protein_accesion)
                subprocess.call('/localdisk/data/BPSM/Assignment2/pullseq -i ' + pro_fasta +
                                ' -n seqname > ' + sys.path[0] + '/fasta/' + protein_accesion, shell=True)
    print('\tDone!')
    print('\tSequences are saved in folder \'fasta\'. Each sequence saved in one file.\n')
    os.remove('seqname')


print('\tStep4.1: Extract each protein sequence to seperate file')
# call the defined funcion fasta_split to split fasta files
fasta_split('protein.fasta')

# Step4.2: patmatmotif scan motif for each protein fasta file
mkdir('motif')  # create file folder for motif scanning result
print('\tStep4.2: Motif scanning...')
# use a loop to scan each fasta
for i in os.listdir(sys.path[0] + "/" + 'fasta' "/"):
    subprocess.call("patmatmotifs -sequence " + sys.path[0] + "/fasta/" + i +
                    " -outfile " + sys.path[0] + "/motif/" + i + " &> log", shell=True)
print('\tDone!')
print("\tThe motif scanning output for each sequence are saved in \"motif\" folder,\n")

# Step4.3: Present motif scanning result, including the number of protein which get motif hit and motif name
# This function is designed to present motif scanning result


def motif_name():
    motif_name = []
    for i in os.listdir(sys.path[0] + "/motif/"):
        motif_result = open(sys.path[0] + "/motif/" + i)
        flag = True  # define flag variable to flag the status of motif hit
        motif = ''  # variable for motif scanning result file
        for lines in motif_result:
            if lines.startswith('# HitCount: 0'):
                flag = False  # if there is no motif hit, change flag to false
            elif lines.startswith('Motif'):
                motif = motif + lines.strip('\s')
                # match motif name via regular expression
                motif_name.append(re.compile(
                    r'Motif = (\S+)').match(lines).group(1))
        if(flag):
            with open(sys.path[0] + '/result/motif_name.txt', 'a') as f:
                f.write('The protein ' + i +
                        ' contains the following motif(s):' + '\n')
                f.write(motif)
    print('There are', len(set(motif_name)), 'motifs:')
    for i in set(motif_name):
        print(i)
        time.sleep(0.5)


print('\tStep4.3: Accessing motif name...')
motif_name()
print('\tDone!')
print("\tPlease access \'motif_name.txt\' file to get motif scanning summary.")
print("------------------------------Step4 Finish--------------------------\n\n")

# step5: wildcard option
print("------------------------------Step5 Start---------------------------")
print("# Step5: wildcard option\n")
mkdir('result/wildcard')
if input('Would you like to perform \'hmoent\' analysis(Enter \'y\' to start analyze)?') == 'y':
    # hmoment: Calculate and plot hydrophobic moment for protein sequence(s)
    subprocess.call(
        'hmoment -seqall protein.fasta -outfile ./result/wildcard/hmoment.result', shell=True)
if input('Would you like to perform \'iep\' analysis(Enter \'y\' to start analyze)?') == 'y':
    # ipe: Calculate the isoelectric point of proteins
    subprocess.call(
        'iep -sequence protein.fasta -graph svg -outfile ./result/wildcard/iep', shell=True)
if input('Would you like to perform \'charge\' analysis(Enter \'y\' to start analyze)?') == 'y':
    # charge: Draw a protein charge plot
    subprocess.call(
        'charge -seqall protein.fasta -outfile ./result/wildcard/charge -graph svg', shell=True)
if input('Would you like to perform \'antigenic\' analysis(Enter \'y\' to start analyze)?') == 'y':
    # antigenic:Find antigenic sites in proteins
    subprocess.call(
        'antigenic -sequence protein.fasta -outfile ./result/wildcard/antigenic -minlen 6', shell=True)
print("------------------------------Step5 Finish---------------------------")

# restruct file tree
os.remove('log')
srcfile = analy_dir+'/motif'
dstfile = analy_dir+'/result/motif'
shutil.move(srcfile, dstfile)
srcfile = analy_dir+'/fasta'
dstfile = analy_dir+'/result/fasta'
shutil.move(srcfile, dstfile)
mkdir('result/blast')
subprocess.call('rm *db.*', shell=True)
subprocess.call('mv *ms* ./result/msa', shell=True)
subprocess.call('mv blast.out ./result/blast', shell=True)
subprocess.call('mv protein.fasta ./result', shell=True)
if os.path.exists("prettyplot.svg"):
    subprocess.call('mv prettyplot.svg ./result/msa', shell=True)
print("------------------------------Analysis Finish---------------------------")
print("--------------Result files are saved in \'result\' folder---------------")
