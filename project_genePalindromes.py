'''
Pseudocode
1.	User defined function – request input sequence (either manual sequence, FASTA or GENBANK file) and minimum palindrome length from user 
2.	Function is created to read and extract sequence in FASTA and GenBank files
3.	Compute the reverse complement for the sequence 
4.  Compute short sequences from the main sequence
5.  Extract all the palindromes from the short sequences
6.  Compute classification of the palindromes into normal palindromes and spacer palindromes

    Normal Palindrome:
    if the palindrome is common with its reverse complement sequence and in even number,  

    Spacer Palindrome:
    if the palindrome is not exactly common with its reverse complement, contains '_' and in odd/even number

7.	Extract and output normal gene palindromes
8.	Extract and output non-repeating palindromes with an intervening spacer region
'''
#import library 
import re 
from prettytable import PrettyTable


#function to request sequence format choice from user
def fInput(): 
    fChoice = int(input("Please choose your input type:\n1 for Raw sequence\n2 for FASTA file\n3 for Genbank file\n"))
    if fChoice == 1: #choice 1 for raw sequence
        seq = input("\nPlease enter your raw sequence: ")
    
    elif fChoice == 2: #choice 2 for fasta
        fname = input("\nPlease enter your FASTA filename: ")
        seq = readFASTA(fname) #call readFASTA function

    elif fChoice == 3: #choice 3 for genbank
        fname = input("\nPlease enter your Genbank filename: ")
        seq = readGB(fname) #call readGB function
    
    else: #fileChoice not in [1,2,3]: #print error message for wrong input
        print("\nInvalid choice. Please try again.","\n")

    seq = seq.upper() #convert sequence to uppercase
    print("\nYour input sequence is:\n",seq,"\n") #print extracted sequence for confirmation
    
    
    #request user to input the desired minimum length for palindrome(s)
    minLength = int(input("Please enter the desired minimum length for your palindrome(s): ",))
    
    
    return seq,minLength #return sequence which is extracted and converted and the min palindrome length

#read in FASTA and GenBank files
def readFASTA(fname):
    with open(fname, 'r+') as FASTAfile:
        seq = '' 

        for line in FASTAfile.readlines():
            if line.startswith('>'):
                pass
            else:
                line = line.strip() #remove spaces at the beginning and at the end of the string
                line = re.sub('\n','',line) #remove \n by replacing it with nothing 
                seq += line #update seq variable with extarcted lines
        return seq #return extracted seq 

def readGB(fname):
    with open(fname, 'r+') as GBfile:
        seq = ''

        for line in GBfile.readlines():
            seqlines = re.search(r"^\s+\d+\s+(.+)",line) #regular expression to match seq lines

            if seqlines:
                g1 = seqlines.group(1) #extract the first sub group
                g1 = re.sub('[\s]','',g1) #remove any whitespace character with nothing
                seq += g1 #update seq variable with extarcted lines
        return seq #return extarcted seq 
 

#Compute the reverse complement seq
def reverseComplement(seq):
    Bases = {'A':'T','C':'G','G':'C','T':'A','_':'_'} #dictionary for base pairs
    comp = ''.join([Bases[i] for i in seq]) #concatenate all complement bases
    rcomp = comp[::-1] #reverses the order
    
    return rcomp  #return reverse complement seq


#Compute short seqs from the main seq
def shortSeq(seq,minLength): 
    seqLength = len(seq) #compute the sequence length
    shortseq = [] #empty list to store all short seqs
    
    #extract short seq and store it in a list
    for i in range(seqLength,minLength-1,-1): #loop from the reverse of sequence until minimum palindrome length
        for k in range(seqLength-i+1): #loop in the length of short sequence
            shortseq.append(seq[k:i+k]) #append the short sequence into list

    return shortseq #return all the short seqs extracted from the sequence

#compute for all palindromes
def allPalindromes(shortseqs): 
    revcompsseq = "" #empty string to store reverse complement of short seq
    allPalindromeMatches = [] #empty list to store all palindromes matches

    for sseq in shortseqs: #for every shortseq in short seqs 
        revcompsseq = reverseComplement(sseq) #compute the reverse complement for shortseq
        
        #check if the first index of both shortseq and its reverse comp matches 
        #skip to the next iteration of shortseq if not match
        if (sseq[0] != revcompsseq[0]): 
            continue 

        #if the first and second index of both shortseq and its reverse comp matches, append it to the all palindrome matches list   
        elif (sseq[0] == revcompsseq[0] and sseq[1] == revcompsseq[1]): 

            #to check if the short seq is already present in the list
            if sseq not in allPalindromeMatches: 
                allPalindromeMatches.append(sseq) #if not, append it to list 
            else: 
                break #break if short seq already present in list                                    
              
    if len(allPalindromeMatches): #true if list is not empty
        return(allPalindromeMatches) #return list that contains all palindrome matches
    else: 
        pass #pass if list is empty 

#compute classification of the palindromes into normal palindromes and spacer palindromes
def classifyPalindrome(palindromes):
    normalP = [] #empty list to store normal palindromes
    spacerP = [] #empty list to store spacer palindromes
    cntnP = 0    #initiate count for normal palindromes
    cntsP = 0    #initiate count for spacer palindromes


    for palindrome in palindromes: 
        #if the palindrome is equal to its reverse complement seq, length is in even number, '_' absent in palindrome and palindrome not in normalPalindrome list
        if palindrome == reverseComplement(palindrome) and (len(palindrome) % 2) == 0 and '_' not in palindrome and palindrome not in normalP:
            normalP.append(palindrome) #append palindrome to nP list
        
        else:
            spacerP.append(palindrome) #else append palindrome to sP list
    
    #output for normal Palindrome
    if len(normalP): #returns true if normalP present
        print("\nNo. of non-spacer palindromes found:", len(normalP))

        # create the table header
        NormalP_table = PrettyTable(["No","Sequence","Indexes","Length"])

        for n in normalP: 
            cntnP += 1         

            match = (re.search(str(n),seq))            

            indexes = str(match.start()+1) +"-" +str(match.end())

            NormalP_table.add_row([str(cntnP),str(n),indexes,str(len(n))]) #add the output into the table
        print(NormalP_table)
        
    else: 
        print("There is no non-spacer palindromes found.\n")
        
        
    #output for spacer Palindrome 
    if len(spacerP): #returns true if spacerP present
        print("\nNo. of palindromes with an intervening spacer region found: ", len(spacerP))
        
        
        # create the table header
        SpacerP_table = PrettyTable(["No","Sequence","Indexes","Length"])

        for n in spacerP: 
            cntsP += 1
            
            match = (re.search(str(n),seq))
            
            indexes = str(match.start()+1) +"-" +str(match.end())
            
            SpacerP_table.add_row([str(cntsP),str(n),indexes,str(len(n))]) #add the output into the table
        print(SpacerP_table)
         
        
    else:
        print("There is no palindromes with an intervening spacer region found.\n")

seq,minLength = fInput() #call function fInput()
revCom = reverseComplement(seq) #call function reverseComplement()
shortseqs = shortSeq(seq,minLength) #call function shortSeq()
palindromes = allPalindromes(shortseqs) #call function allPalindromes()
classifyPalindrome(palindromes) #call function classifyPalindrome()


