import pandas as pd
from sys import argv

# list of transcript IDs to include
ids_in = argv[1]

# gtf file
gtf_in = argv[2]

# save transcript ids to list
transcript_ids = []
for line in open(ids_in):
    transcript_ids.append(line.strip())

# read gtf data to a datafram
gtf_data = pd.read_csv(gtf_in, names = ['Scaff', 'Prot','Feature', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute'], sep = '\t')

# split gene info to get transcript id
gtf_data['Transcript_id'] = (gtf_data['Attribute'].str.split(' ',  expand = True)[1] + "_" + gtf_data['Attribute'].str.split(' ',  expand = True)[5]).str.replace(';','').str.replace('"','')

# subset gtf for only cds
cds_data = gtf_data[gtf_data['Feature']=='CDS']


for trans in transcript_ids:

    # get current transcript from cds data
    transcript_data = cds_data[cds_data['Transcript_id']==trans]

    # only work with transcripts on a single scaffold
    if len(transcript_data.Scaff.unique()) == 1:
        scaff = transcript_data.iloc[0]['Scaff']
        strand = transcript_data.iloc[0]['Strand']

        # get start and stop positions
        start = transcript_data.iloc[0]['Start'] - 1
        end = transcript_data.iloc[-1]['End']

        
        name = transcript_data.iloc[0]['Transcript_id']
        score = transcript_data.iloc[0]['Score']
        strand = transcript_data.iloc[0]['Strand']
        thickstart = start
        thickend = start
        itemRgb = 0
        blockcount = len(transcript_data)

        blocksizes = ""
        blockstarts = ""

        # get data on each exon
        for exon in range(0, len(transcript_data)):

            # get size and start of each exon    
            ex_size = str(transcript_data.iloc[exon]['End'] - (transcript_data.iloc[exon]['Start']-1))
            ex_start = str((transcript_data.iloc[exon]['Start']-1) - start)

            # add exon sizes to blocksize
            blocksizes += ex_size
            blocksizes += ","

            # add exon starts to blockstarts
            blockstarts += ex_start
            blockstarts += ","

        print(scaff, start, end, name, score, strand, thickstart,thickend,itemRgb, blockcount, blocksizes, blockstarts, sep = '\t')
