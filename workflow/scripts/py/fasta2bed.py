import pandas as pd
import sys

inFasta = snakemake.input[0] # take fasta as command argument

def fastaParser(fasta):
    headers = []
    with open(fasta) as f:
        header = None
        for line in f:
            if line.startswith('>'): # identifies fasta header line
                headers.append(line[1:-1]) # append all of the line that isnt >
                header = line[1:] # in reset header
    newHeader = (header.replace(':',',') for header in headers) # format to be accepted later
    newnewHeader = (header.replace('-',',') for header in newHeader) # format to accept later
    bedHead = (header.split(',') for header in newnewHeader) # separate by comma from format above
    return bedHead

fastaParsed = fastaParser(inFasta) # run the function

headers = list(fastaParsed) # go from generator to list
bedfile = pd.DataFrame(headers) # create dataframe which will output bedfile

bedfile.to_csv(snakemake.output[0], sep='\t', index=False, header=None) # output bedfile format
