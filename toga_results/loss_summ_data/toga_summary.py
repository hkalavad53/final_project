import sys
from collections import defaultdict
dataDict = defaultdict(dict)
specList=[]

# Check if command line arguments are provided
if len(sys.argv) < 2:
    print("Usage: python script.py file1.tsv,file2.tsv,file3.tsv")
    sys.exit(1)

# Get the comma-separated list of file paths from command line arguments
fileList = sys.argv[1].split(',')

for infile in fileList:
    spec=infile.split('.')[0]
    specList.append(spec)
    for line in open(infile):
        f=line.strip('\n').split('\t')
        if f[0]=="GENE":
            gene,result=f[1],f[2]
            dataDict[gene][spec]= result

with open('output.tsv', 'w') as outfile:
    header = "GENE\t" + "\t".join(specList)
    print(header, file=outfile)
    for gene in dataDict:
        row = gene
        for species in specList:
            result = dataDict[gene].get(species)
            if result is None:
                result = 'NA'
            row += "\t" + result
        print(row, file=outfile)


