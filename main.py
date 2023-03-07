#!/usr/bin/python
__author__ = "David Zhang"
__email__ = "david.zhang.ddz5@yale.edu"
__copyright__ = "Copyright 2023"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm

import argparse 
import numpy as np
import pandas as pd

### Read in arguments from python, require an input sequence text file, a similarity matrix text file, open gap and extension gap scores 
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()

# Implement a trace class to keep track of back trace after finding optimal alignment 
class Trace():
    DIAGONAL = 3
    LEFT = 2
    UP = 1
    STOP = 0

### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
    ### calculation
    
    # Obtain similarity matrix for alignment and mismatches 

    simMatrix = calcSimMatrix(scoreFile) 

    # Sequences to align 
    seqTxt = open(inputFile, "r")
    seqs = seqTxt.readlines() 

    if (len(seqs) != 2):
        print("Please provide two sequences")
        return 1
    
    # First sequence provided 
    refSeq = seqFilter(seqs[0])

    # Secondary sequence 
    sseq = seqFilter(seqs[1]) 
    scoreMatrix, backTrace, max_score, max_index = SW(simMatrix, refSeq, sseq, openGap, extGap) 

    # Write to file 

    output = open("output-ddz5.txt", "w") 
    # Delete any existing contents in file 
    output.seek(0)
    output.truncate(0)

    output.write("-----------\n")
    output.write("|Sequences|\n")
    output.write("-----------\n")
    
    output.write("sequence1\n")
    output.write(refSeq)
    output.write("\n")
    output.write("sequence2\n")
    output.write(sseq)
    output.write("\n")

    output.write("--------------\n")
    output.write("|Score Matrix|\n")
    output.write("--------------\n")

    # Dimensions of score matrix 
    col = len(refSeq) + 1
    row = len(sseq) + 1

    output.write("\t\t")

    # Row sequence 
    for i in range(len(refSeq)):
        output.write(refSeq[i])
        output.write("\t")

    output.write("\n\t")

    # Similarity score after alignment  

    # Row 1 (all zeros)

    for j in range(col): 
        output.write(str(int(scoreMatrix[0, j])))
        output.write("\t")

    # Similarity scores for all remaining rows 
    for i in range(1, row):
        # Column sequence 
        output.write("\n" + sseq[i - 1] + "\t")
        # Ouput scores 
        for j in range(col): 
            output.write(str(int(scoreMatrix[i, j])))
            output.write("\t")

    # Alignment Results
    output.write("\n----------------------\n")
    output.write("|Best Local Alignment|\n")
    output.write("----------------------\n")

    output.write("Alignment Score:")
    output.write(str(int(max_score)) + "\n")

    output.write("Alignment Results:\n")

    # Print out traceback

    (max_i, max_j) = max_index

    final_aligned_seq1 = ""
    final_aligned_seq2 = "" 

    # sseq is number of rows 
    # Non-matching zone of sequences (i.e. base pairs that come after "()")
    final_aligned_seq1 += ")"
    final_aligned_seq1 += sseq[max_i:]
    
    # refSeq is number of cols 
    final_aligned_seq2 += ")" 
    final_aligned_seq2 += refSeq[max_j:]

    # Matching zone of sequences 
    current_aligned_seq1 = ""
    current_aligned_seq2 = ""

    while backTrace[max_i, max_j] != Trace.STOP and max_i > 0 and max_j > 0: 
        if backTrace[max_i, max_j] == Trace.DIAGONAL:
            # Due to the insertion of extra row and column when calculating scoreMatrix, our indicies are one-shifted 
            current_aligned_seq1 += sseq[max_i - 1]
            current_aligned_seq2 += refSeq[max_j - 1]
            max_i -= 1
            max_j -= 1

        # Horizontal reference (refSeq) aligns to a gap in vertical (sseq) 
        elif backTrace[max_i, max_j] == Trace.LEFT: 
            current_aligned_seq1 += "-"
            current_aligned_seq2 += refSeq[max_j - 1]
            max_j -= 1

        # Vertical reference (sseq) aligns to a gap in horizontal (refSeq)
        elif backTrace[max_i, max_j] == Trace.UP: 
            current_aligned_seq1 += sseq[max_i - 1]
            current_aligned_seq2 += "-"
            max_i -= 1
    
    # Reverse the sequences 
    current_aligned_seq1 = current_aligned_seq1[::-1]
    current_aligned_seq2 = current_aligned_seq2[::-1]

    # Append the non-matching zone with matching zone  
    final_aligned_seq1 = current_aligned_seq1 + final_aligned_seq1
    final_aligned_seq2 = current_aligned_seq2 + final_aligned_seq2
 
    # Append front non-matching zone  
    final_aligned_seq2 = refSeq[0:max_j] + "(" + final_aligned_seq2
    final_aligned_seq1 = sseq[0:max_i] + "(" + final_aligned_seq1

    # Post-Processing Spaces
    # Spaces in front of shorter head 
    if (max_i < max_j):
        while (max_i != max_j): 
            final_aligned_seq1 = " " + final_aligned_seq1
            max_j -= 1
    elif (max_i > max_j):
        while (max_i != max_j): 
            final_aligned_seq2 = " " + final_aligned_seq2
            max_i -= 1

    # Spaces behind shorter tail 
    
    diff = max(len(final_aligned_seq1), len(final_aligned_seq2)) - min(len(final_aligned_seq1), len(final_aligned_seq2))

    if (len(final_aligned_seq1) < len(final_aligned_seq2)):
        while diff != 0:
            final_aligned_seq1 += " "
            diff -= 1
    else: 
        while diff != 0:
            final_aligned_seq2 += " "
            diff -= 1

    output.write(final_aligned_seq2 + "\n")
    
    # Output matches 

    seq_index = 0
    match = False 
    # Should be same length but iterate across both sequences and compare 
    while (seq_index != max(len(final_aligned_seq1), len(final_aligned_seq2))): 
        # Check for when entering match zone (denoted by "()")
        if (final_aligned_seq1[seq_index] == "("):
            match = True
            output.write(" ")
        else: 
            # Check for end of matching zone 
            if (final_aligned_seq1[seq_index] == ")"): 
                match = False
                output.write(" ")
            # Matching zone 
            elif match:
                if (final_aligned_seq1[seq_index] == final_aligned_seq2[seq_index]):
                    output.write("|")
                else:
                    output.write(" ")
            else:
                output.write(" ")
        seq_index += 1
    
    output.write("\n")
    output.write(final_aligned_seq1 + "\n")    
    output.close()    

### Calculates the score matrix given two sequences and a similarity matrix 
def SW(simMatrix, seq1, seq2, openGap, extGap):
    # seq1 is horizontal, seq2 is vertical  
    col = len(seq1) + 1
    row = len(seq2) + 1
    scMatrix = np.zeros(shape = (row, col))
    backTrace = np.zeros(shape = (row, col)) 

    # Initialize first row and column of backtrace matrix 
    backTrace[0] = [Trace.LEFT] * col
    backTrace[:, 0] = [Trace.UP] * row
    backTrace[0, 0] = Trace.STOP

    # Initialize index for cell with max score 
    max_score = -1 
    max_index = (-1, -1) 

    # Calculate score for each cell according to SW algorithm 
    
    for i in range(1, row):
        for j in range(1, col): 
            # Calculate all possible scores and take the max 

            # indicies are one ahead of sequence base pairs due to extra row and column 
            match_value = simMatrix[seq1[j - 1].upper()][seq2[i - 1].upper()]

            # Value coming from mis(match) and the previous diagonal cell 
            diagonal = scMatrix[i - 1, j - 1] + match_value 

            # Calculate the maximum score across all possible horizontal gaps leading to cell 
            
            # Keep track of max gap score 
            maxLeft = 0
            # Keep track of how far left
            numLeft = 0
            for x in range(0, j): 
                gap = openGap + extGap * (j - 1 - x) 
                leftScore = scMatrix[i, x] + gap
                if (leftScore >= maxLeft): 
                    maxLeft = leftScore
                    # Keep track of how far left 
                    numLeft = j - x

            # Calculate the maximum score across all possible vertical gaps leading to cell  

            # Keep track of max gap score 
            maxUp = 0
            # Keep track of how far up
            numUp = 0
            for y in range(0, i):                                                                                                                         
                gap = openGap + extGap * (i - 1 - y)
                upScore = scMatrix[y, j] + gap 
                if (upScore >= maxUp): 
                    maxUp = upScore
                    numUp = i - y
             
            scMatrix[i, j] = max(0, diagonal, maxLeft, maxUp)

            # Implement trace back

            # Creating backTrace 
            if scMatrix[i, j] == 0:
                backTrace[i, j] = Trace.STOP

            elif scMatrix[i, j] == diagonal: 
                backTrace[i, j] = Trace.DIAGONAL

            elif scMatrix[i, j] == maxLeft:
                # How large of a horizontal gap 
                for step in range(numLeft): 
                    backTrace[i, j - step] = Trace.LEFT

            elif scMatrix[i, j] == maxUp: 
                # How large of a vertical gap 
                for step in range(numUp): 
                    backTrace[i - step, j] = Trace.UP
            
            # Updating max score cell 
            if scMatrix[i, j] >= max_score: 
                max_score = scMatrix[i, j]
                max_index = (i, j)

    return scMatrix, backTrace, max_score, max_index

### Returns the similarity matrix as specified by user input
def calcSimMatrix(scoreFile): 
    # Read in input scoreFile 
    scoreTxt = open(scoreFile, "r")

    lines = scoreTxt.readlines()

    # Obtain all keys from first column of matrix
    colkeys = [] 
    for char in lines[0]:
        if char.isalpha():
            colkeys.append(char.upper())

    # Create score dictionary to index into 
    simDict = {}

    # Iterate through each row in matrix
    for row in lines[1:]:
        # Only process valid rows (i.e. row length must be greater than or equal to the number of colkeys)

        if (len(row) >= len(colkeys)):
            vals = {}

            # First occurrence of a letter denotes the beginning of the scores for that letter's row 
            rowkey = ""
            fflag = True 

            # Number of scores in each row 
            count = 0 

            # Flag to check for negative values 
            neg = False 

            # Flag to check for beginning of an integer 
            iflag = False 

            val = "" 

            # Process each char in respective row 
            for char in row:
                # First (and only, if valid input) letter in the row 
                if char.isalpha() and fflag:
                    rowkey = char
                    fflag = False 
                
                # Negative value 
                elif char == "-": 
                    neg = True 

                # Start of an integer 
                elif char.isdigit(): 
                    val += char
                    iflag = True 

                # A space marks the end of an integer 
                elif char.isspace() and iflag: 
                    # Negative values 
                    if (neg):
                        vals[colkeys[count]] = -1 * int(val) 
                        neg = False 
                    else: 
                        vals[colkeys[count]] = int(val) 
                    count += 1
                    iflag = False 
                    val = ""
    
            simDict[rowkey.upper()] = vals  
    return simDict

### Return the properly formatted sequence 
def seqFilter(seq):
    finalSeq = ""
    # Only append alphabetical chars 
    for char in seq:
        if char.isalpha(): 
            finalSeq += char.upper()
    return finalSeq

### Run your Smith-Waterman Algorithm
if __name__ == "__main__":
    runSW(args.input, args.score, args.opengap, args.extgap)