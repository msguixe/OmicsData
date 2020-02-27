###############################################################################
# IMPORTS
###############################################################################
import sys
import time
from Bio import SeqIO
import mysql.connector
from mysql.connector import (connection)
from mysql.connector import errorcode
from Bio.SubsMat import MatrixInfo
from decimal import *

###############################################################################
# OPTIONS
###############################################################################

##READ_FILE_FORMAT##########################################
formatalign = "fasta"

##GAP_PENALTIES### Values can be adjusted by user from 0 to 100 (Default: Po = 100 and Pe = 0.2)
Po = 10
Pe = 0.2

###############################################################################
# FUNCTIONS
###############################################################################
def readfilealign(formatalign, intent, option):
    '''
    Read a file that contain a fasta aligment.
    file --> list
    '''
    filename = input('Write the name of the file that contain the MSA: \n > ')
    try:
        records = list(SeqIO.parse((open(filename + '.'+ formatalign,'rU')), formatalign)) #Turns a sequence file into an iterator returning SeqRecords.
    except:
        print ("An error occur! Format or type are indicated wrong!")
        sys.exit()
    return records

def mysql_identification(seqs, x):
    '''
    Creates the identifications for each sequence alignment in pairs.
    '''
    split = seqs[x].id.split('|') #['sp', 'P04440', 'DPB1_HUMAN']
    Code_Uniprot = split[1] #'P04440'
    Identification = split[2] #'DPB1_HUMAN'
    return Code_Uniprot, Identification

def mysql_database(cursor, n_database):
    '''
    Creates the database where table of results will be stored.
    '''
    mysql_database ='CREATE DATABASE {}'.format(n_database)
    cursor.execute(mysql_database)

def mysql_table(n_table, n_database, paswrd):
    '''
    Connect to database (if not exists creates new one) and creates
    the table of results
    '''
    cnx = mysql.connector.connect(user='root', password= paswrd, host='localhost')
    cursor = cnx.cursor()
    try:
        cnx.database = n_database
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_BAD_DB_ERROR:
            mysql_database(cursor, n_database)
            cnx.database = n_database
        else:
            print(err)
            #exit(1)
    mysql_table = '''CREATE TABLE {} (Entry varchar(10) NOT NULL, Entry_code_1 varchar(100) NOT NULL, Entry_name_1 varchar(100) NOT NULL, Entry_code_2 varchar(100) NOT NULL, Entry_name_2 varchar(100) NOT NULL, Identity float(10) NOT NULL, Normalised_Global_Similarity_Score float(10) NOT NULL, PRIMARY KEY (Entry)) ENGINE = InnoDB'''.format(n_table) #Create table function in mysql
    try:
        print("Creating table {}: ".format(n_table), end='') #Print a message indicating that table is created.
        cursor.execute(mysql_table) #Execute the function in mysql
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
            print("Already exists.") #In case that the talbe already exists print this message
        else:
            print(err.msg)
    else:
        print("OK") #If table is created without problem print this.

    cursor.close()
    cnx.close()

def identity(seq1, seq2):
    '''
    Function that returns the normalised total score similarity between two sequences.
    str --> float
    '''
    ident_score = 0
    for i in range(len(seq1)): #Compare each aminoacid in the two seq with position i
        amino1 = seq1[i]
        amino2 = seq2[i]
        if amino1 == amino2:
            ident_score += 1
    getcontext().prec = 5 #Module decimal: Provides support for fast correctly-rounded decimal floating point arithmetic.
    identity = (Decimal(ident_score)/Decimal(len(seq1)))*100 #Is the identity formula. Elements need the word Decimal(float)
    return identity

def score(seq1, seq2, n_matrix):
    '''
    Function that returns the score between two sequences using a substitution matrix.
    str --> int
    '''
    for key in list(n_matrix): #This part add to the dictionary the scores if the compare is inverted. example --> matrix = {('K','L'):-2, ('L','K'):-2,...}
        n_matrix[(key[1],key[0])] = n_matrix[key]
    score = 0
    for i in range(len(seq1)):
        amino1 = seq1[i]
        amino2 = seq2[i]
        l_key = [] #For clean list
        l_key.append(amino1) #Add the amino in pos 0
        l_key.append(amino2) #Add the amino in pos 1 l_key = ['K','L']
        if amino1 != '-' and amino2 !='-': #Discard alignments like: M -, - M, - -
            score += n_matrix[tuple(l_key)]#Transform list to tuple for search the value how a key in substituion matrix
        else:
            continue
    return score

def gap_penalty (seqs):
    '''
    Returns the total number of opening gaps and extension gaps.
    list --> int
    '''
    j = 0
    ext_gap = 0
    gap = 0
    num_gaps_len = 0
    num_gaps = []
    while j < len(seqs[0].seq): #j iterates through the positions of each sequence (same length)
        i = 0
        while i < len(seqs):
            if seqs[i].seq[j] == '-': #If program found a gap
                gap = 1
                break
            i += 1
        if gap == 1: #In case the program found a gap, start to count extension and opening gap
            num_gaps_len  += 1 #
            ext_gap += 1 #Gap extension
        else:
            if num_gaps_len  != 0: #In case gap extension stops, add the len of the gap
                num_gaps.append(num_gaps_len)
                num_gaps_len  = 0
        j += 1
        gap = 0
        if j == len(seqs[0].seq) and num_gaps_len != 0: #All seq have the same len. If the gap extension is on the end.
            num_gaps.append(num_gaps_len)
    return len(num_gaps), ext_gap

def similarity (seq1, seq2, n_matrix, seqs, Po, Pe,same):
    '''
    Function that returns the normalised total score similarity between two sequences using a substitution matrix.
    str,float,boolean --> float
    '''
    sim_score = score(seq1, seq2, n_matrix)
    equ_score = score(seq1, seq1, n_matrix)
    num_gaps, ext_gap = gap_penalty(seqs)
    getcontext().prec = 3
    if same == True: #Change the operation, depends if sequences are equal or not
        similarities = Decimal(sim_score)/Decimal(equ_score) #Eliminates the gap penalty
    else:
        similarities = Decimal(sim_score - num_gaps*Po - ext_gap*Pe)/Decimal(equ_score) 	
    return similarities

def mysql_results(k, n_table, n_database, paswrd,Entry_code_1, Entry_name_1, Entry_code_2, Entry_name_2, identities, similarities):
    '''
    Function that save the results on mysql format.
    int,str,float --> mysql database
    '''
    cnx = mysql.connector.connect(user='root', password= paswrd, host='localhost') #Connect to mysql
    cursor = cnx.cursor() #Is an executor of commands on mysql from python.
    try:
        Entry = '{} -'.format(k) #Each aligment comparation in pairs get a number (PRIMARY KEY)
        add_result = ('INSERT INTO {}.{}(Entry, Entry_code_1, Entry_name_1, Entry_code_2, Entry_name_2, Identity, Normalised_Global_Similarity_Score) VALUES (%s, %s, %s, %s, %s, %s, %s)'.format(n_database,n_table)) #Insert data organization function for mysql
        data_result = (Entry, Entry_code_1, Entry_name_1, Entry_code_2, Entry_name_2, identities, similarities) #Results that go on the table
        cursor.execute(add_result, data_result) #Insert the information in the table created
    except mysql.connector.IntegrityError as err: #If occur some error
        print('Error: {}'.format(err))
    cnx.commit() #Make sure data is committed to the database
    cursor.close()
    cnx.close()

###############################################################################
# MAIN PROGRAM
###############################################################################
print ('----{ WELCOME TO ISC }---- \n \n || INFORMATION TO USER ||')
print('\n - ISC perhaps a calculation of the identity and similarity between two sequences.')
print('\n - ISC can use a multiple sequence alignment file as input.')
print('\n - Gap penalties parameters can be adjusted by user from 0 to 100, editing the main program file. (Default: Popen = 100 and Pextension  = 0.2)')
print('\n - ISC use python 3.5.2 under Ubuntu 16.04.1 and mysql server version 5.7.16-0ubuntu0.16.04.1 (Ubuntu).\n') #Info about the program.
time.sleep(20)
print ('ISC: Sir/Lady, I will ask some things for do my work. Wait \n')
time.sleep(1)

option = 0
intent = 1
while option == 0:
    try:
        seqs = readfilealign(formatalign, intent, option)
        option = 1
    except:
        print ('Wrong name of the file or is not in the directory! You use {}/3 tries.'.format(intent))
        intent += 1
        if intent == 4:
            sys.exit('File not found!.') #In case that user don't introduce well the name.

paswrd = input('Mysql password: \n > ') #Password of the user in mysql
n_database = input('Mysql database name: \n > ') #Name of the database where you like save the table with results
n_table = input('Mysql table name (contain results organized): \n > ') #Table that contain the results
mysql_table(n_table, n_database, paswrd)

option = 0
intent = 1
while option == 0:
    select = input('Select the substitution matrix (blosum62/pam250): \n > ') #User have two choices
    if select == 'blosum62':
        n_matrix = MatrixInfo.blosum62
        option = 1
    elif select == 'pam250':
        n_matrix = MatrixInfo.pam250
        option = 1
    elif intent == 4:
        sys.exit('Matrix name not corresponding with the options!') #If user don't introduce a valid name of matrix
    else:
        print ('Matrix {} doesn\'t exist!  You use {}/3 tries.'.format(select, intent))
        intent += 1

k = 0
for i in range(len(seqs)): #Compare each sequence with all sequences in pairs (i and j)
    for j in range(i,len(seqs)):
        if seqs[i].id == seqs[j].id: #In case that the two compared sequences are the same do a different operation
            same = True
        else:
            same = False
        identities = identity(seqs[i].seq,seqs[j].seq)
        similarities  = similarity (seqs[i].seq, seqs[j].seq, n_matrix, seqs, Po, Pe,same)
        Entry_code_1, Entry_name_1 = mysql_identification (seqs,i)
        Entry_code_2, Entry_name_2 = mysql_identification (seqs,j)
        mysql_results(k, n_table, n_database, paswrd,Entry_code_1, Entry_name_1, Entry_code_2, Entry_name_2, identities, similarities)
        k += 1 #Enumerates each comparation between two sequences. k --> Entry

print ('You found the results on mysql/' + n_database + '/' + n_table)
print ('Thank you for your patient')