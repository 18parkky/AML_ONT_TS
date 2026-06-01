import os, time, gzip, pickle
import pysam, pandas as pd, intervaltree

def load_vcf( PATH_vcf, isGzip=True ):
    skiplines = 0 
    if isGzip == True: 
        open_func = gzip.open 
    else:
        open_func = open
    with open_func(PATH_vcf, 'rt') as vcf:
        for line in vcf:
            if line[:2] == '##':
                skiplines += 1
            else:
                break 
    vcf = pd.read_csv(PATH_vcf, sep='\t', skiprows=skiplines)
    return vcf

def loadFromPickle(dir_pickle):
    with open(dir_pickle, 'rb') as handle:
        unserialized_pickle = pickle.load(handle)
    return unserialized_pickle

def parse_info(info_str, keys=['ALLELEID', 'CLNSIG', 'CLNHGVS', 'CLNVI', 'RS']):
    info_dict = dict(
        entry.split('=', 1)
        for entry in info_str.split(';')
        if '=' in entry
    )
    return pd.Series({ key : info_dict.get(key, 'N/A') for key in keys})

def elapsedTime( starting_time ):
    return time.time()-starting_time

def strip_chr_prefix(chrom):
    return chrom[3:] if chrom.startswith("chr") else chrom

def segment_chromosome_batch(chrom_list, path_bam, name, dir_out, segment_length_bp):
    bamfile = pysam.AlignmentFile(path_bam, "rb")

    for chrom in chrom_list:
        output_name = f"{name}_{chrom}.segmented.fq.gz"
        output_path = os.path.join(dir_out, output_name)
        out_fh = gzip.open(output_path, 'wt')

        num_reads = 0

        for read in bamfile.fetch(chrom):
            if read.is_secondary or read.is_supplementary:
                continue

            seq = read.query_sequence
            qual = read.query_qualities

            if not seq or not qual or len(seq) < 2 * segment_length_bp:
                continue

            for idx, i in enumerate(range(0, len(seq), segment_length_bp)):
                seg_seq = seq[i:i + segment_length_bp]
                seg_qual = qual[i:i + segment_length_bp]

                if len(seg_seq) < segment_length_bp:
                    continue

                seg_qual_str = ''.join(chr(q + 33) for q in seg_qual)
                header = f"@{read.query_name}_segment{idx+1}"
                out_fh.write(f"{header}\n{seg_seq}\n+\n{seg_qual_str}\n")
                num_reads += 1

        out_fh.close()
        if num_reads == 0: os.remove(output_path)

    bamfile.close()
    
def processCigarString( cigarstring ):
    list_str_idx = list()
    for idx, elem in enumerate(cigarstring):
        try:
            int(elem)
        except ValueError:
            list_str_idx.append(idx)
    
    anchor_idx = -1
    largest_ins, largest_del = 0, 0
    largest_sclip = 0
    
    for str_idx in list_str_idx:
        cigaroperation = cigarstring[str_idx]
        operationMagnitude = int(cigarstring[anchor_idx+1:str_idx])
        anchor_idx = str_idx
        
        if cigaroperation == 'I':
            if largest_ins < operationMagnitude: largest_ins = operationMagnitude
        elif cigaroperation == 'D':
            if largest_del < operationMagnitude: largest_del = operationMagnitude
        elif cigaroperation == 'S':
            if largest_sclip < operationMagnitude: largest_sclip = operationMagnitude
    
    return largest_ins, largest_del, largest_sclip

def SummarizeAlignment_targetGenes( PATH_bamfile, BED_slice, process_name, DIR_temp ):
    bamfile = pysam.AlignmentFile( PATH_bamfile, 'rb', )

    # ~26 minutes for bc11, no multiprocessing
    AlignmentSummary = list()

    for tup in BED_slice.itertuples():    
        
        chrom       = str(tup.chrom)
        start_pos   = int(tup.start)
        end_pos     = int(tup.end)
        
        for read in bamfile.fetch(chrom, start_pos, end_pos, multiple_iterators=True):
            if read.is_secondary == True or read.is_supplementary == True:
                continue 
            
            OriginalReadName    = read.query_name.split('_')[0]
            SegmentNumber       = int(read.query_name.split('_')[1].split('segment')[1])
            ReferenceName       = read.reference_name
            ReferencePositions  = read.get_reference_positions()
            CigarString         = read.cigarstring
            MAPQ = read.mapping_quality
            flag = read.flag
            LargestIns, LargestDel, LargestSclip = processCigarString(CigarString)
            
            AlignmentSummary.append( [OriginalReadName, SegmentNumber, ReferenceName, 
                                        ReferencePositions[0], ReferencePositions[-1], MAPQ, flag,
                                        LargestIns, LargestDel, LargestSclip, CigarString, ] ) 
        
    AlignmentSummary = pd.DataFrame(AlignmentSummary, 
                                    columns=['ReadName', 'SegmentNumber', 'ReferenceName', 
                                            'AlignedStartPos', "AlignedEndPos", 'MAPQ', 'flag',
                                            'LargestIns', 'LargestDel', 'LargestSclip', 'cigar',])
    # AlignmentSummary.sort_values(['ReadName', 'SegmentNumber'], inplace=True)
    if len(AlignmentSummary) == 0: return 
    
    AlignmentSummary.to_csv(f'{DIR_temp}/temp.{process_name}.tsv', sep='\t', index=False)
    return 

def AnnotateAlignedRegion( AlignmentSummary, bed, padding=10000 ):
    
    # Build tree
    trees = {} 
    for tup in bed.itertuples():
        chrom = tup.chrom
        try:
            trees[chrom]
        except KeyError:
            trees[chrom] = intervaltree.IntervalTree()
            
        trees[chrom][int(tup.start)-padding:int(tup.end)+1+padding] = tup.gene
    
    # Annotate
    annot = list()
    for tup in AlignmentSummary.itertuples():
        try:
            SearchRes = sorted(trees[tup.ReferenceName][tup.AlignedStartPos])    
        except KeyError:
            annot.append('N/A')
            continue
        
        if len(SearchRes) == 0:
            annot.append( 'N/A' )
        elif len(SearchRes) == 1:
            annot.append( SearchRes[0][2])
        else:
            genes = "::".join([ IT[2] for IT in SearchRes ])
            annot.append(genes)
            
    AlignmentSummary['AlignedRegion'] = annot
    return 

def has_large_indel(cigar_tuples, allowed_gap):
    max_op = None
    max_len = 0

    for op, length in cigar_tuples:
        if op in (1, 2, 3):  # 1=I, 2=D, 3=N
            if length > allowed_gap and length > max_len:
                max_len = length
                max_op = {1: 'insertion', 2: 'deletion', 3: 'refskip'}[op]

    if max_op is not None:
        return max_op, max_len, True
    else:
        return None, None, False
    
def LoadGencodeGFF( PATH_to_GencodeGFF, featureOfInterest='exon', skiprows=7, ):
    GencodeGFF = pd.read_csv(PATH_to_GencodeGFF, sep='\t', skiprows=skiprows, header=None)
    GencodeGFF.columns = ['sequence', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']

    GencodeGFF.dropna(inplace=True)
    if featureOfInterest == None:
        pass 
    else:
        GencodeGFF = GencodeGFF[(GencodeGFF['feature']==featureOfInterest)]

    col1, col2, col3, col4 = list(), list(), list(), list()

    for tup in GencodeGFF.itertuples():
        try:
            attributes = tup.attributes.split(';')
        except:
            print(tup)
            raise ValueError
        
        ENSG, ENST = None, None
        GeneName, TranscriptName = None, None
        
        for attribute in attributes:
            if 'gene_id' in attribute:
                ENSG = attribute.strip().split(" ")[1].split('.')[0][1:]
                col1.append(ENSG.split('.')[0][1:])
                
            elif 'transcript_id' in attribute:
                ENST = attribute.strip().split(" ")[1].split('.')[0][1:]
                col2.append(ENST.split('.')[0][1:])

            elif 'gene_name' in attribute:
                GeneName = attribute.strip().split(" ")[1].upper()[1:-1]
                col3.append(GeneName)

            elif 'transcript_name' in attribute:
                TranscriptName = attribute.strip().split(" ")[1][1:-1]
                col4.append(TranscriptName)

        if ENSG == None:
            col1.append(None)
        if ENST == None:
            col2.append(None)
        if GeneName == None:
            col3.append(None)
        if TranscriptName == None:
            col4.append(None)
        
    GencodeGFF['ENSG'] = col1
    GencodeGFF['ENST'] = col2
    GencodeGFF['GeneName'] = col3
    GencodeGFF['TranscriptName'] = col4
    
    return GencodeGFF

def processCigarString( cigarstring ):
    list_str_idx = list()
    for idx, elem in enumerate(cigarstring):
        try:
            int(elem)
        except ValueError:
            list_str_idx.append(idx)
    
    anchor_idx = -1
    largest_ins, largest_del = 0, 0
    largest_sclip = 0
    
    for str_idx in list_str_idx:
        cigaroperation = cigarstring[str_idx]
        operationMagnitude = int(cigarstring[anchor_idx+1:str_idx])
        anchor_idx = str_idx
        
        if cigaroperation == 'I':
            if largest_ins < operationMagnitude: largest_ins = operationMagnitude
        elif cigaroperation == 'D':
            if largest_del < operationMagnitude: largest_del = operationMagnitude
        elif cigaroperation == 'S':
            if largest_sclip < operationMagnitude: largest_sclip = operationMagnitude
    
    return largest_ins, largest_del, largest_sclip

def SummarizeAlignment_fullscan( PATH_bamfile, chrom, region, process_name, DIR_temp ):
    # ~26 minutes for bc11, no multiprocessing
    AlignmentSummary = list()
    
    start_pos   = region[0]
    end_pos     = region[1]
    
    bamfile = pysam.AlignmentFile( PATH_bamfile, 'rb' )
    
    for read in bamfile.fetch(chrom, start_pos, end_pos, multiple_iterators=True):
        if read.is_secondary == True or read.is_supplementary == True:
            continue 
        
        OriginalReadName    = read.query_name.split('_')[0]
        SegmentNumber       = int(read.query_name.split('_')[1].split('segment')[1])
        ReferenceName       = read.reference_name
        ReferencePositions  = read.get_reference_positions()
        CigarString         = read.cigarstring
        MAPQ = read.mapping_quality
        flag = read.flag
        LargestIns, LargestDel, LargestSclip = processCigarString(CigarString)
        
        AlignmentSummary.append( [OriginalReadName, SegmentNumber, ReferenceName, 
                                    ReferencePositions[0], ReferencePositions[-1], MAPQ, flag,
                                    LargestIns, LargestDel, LargestSclip, CigarString, ] ) 
    
    AlignmentSummary = pd.DataFrame(AlignmentSummary, 
                                    columns=['ReadName', 'SegmentNumber', 'ReferenceName', 
                                            'AlignedStartPos', "AlignedEndPos", 'MAPQ', 'flag', 
                                            'LargestIns', 'LargestDel', 'LargestSclip', 'cigar',])
    if len(AlignmentSummary) == 0: return 
    
    AlignmentSummary.to_csv(f'{DIR_temp}/temp.{chrom}.{process_name}.tsv', sep='\t', index=False)
    return 

def SummarizeAlignment_targetGenes( PATH_bamfile, BED_slice, process_name, DIR_temp ):
    bamfile = pysam.AlignmentFile( PATH_bamfile, 'rb', )

    # ~26 minutes for bc11, no multiprocessing
    AlignmentSummary = list()

    for tup in BED_slice.itertuples():    
        
        chrom       = str(tup.chrom)
        start_pos   = int(tup.start)
        end_pos     = int(tup.end)
        
        for read in bamfile.fetch(chrom, start_pos, end_pos, multiple_iterators=True):
            if read.is_secondary == True or read.is_supplementary == True:
                continue 
            
            OriginalReadName    = read.query_name.split('_')[0]
            SegmentNumber       = int(read.query_name.split('_')[1].split('segment')[1])
            ReferenceName       = read.reference_name
            ReferencePositions  = read.get_reference_positions()
            CigarString         = read.cigarstring
            MAPQ = read.mapping_quality
            flag = read.flag
            LargestIns, LargestDel, LargestSclip = processCigarString(CigarString)
            
            AlignmentSummary.append( [OriginalReadName, SegmentNumber, ReferenceName, 
                                        ReferencePositions[0], ReferencePositions[-1], MAPQ, flag,
                                        LargestIns, LargestDel, LargestSclip, CigarString, ] ) 
        
    AlignmentSummary = pd.DataFrame(AlignmentSummary, 
                                    columns=['ReadName', 'SegmentNumber', 'ReferenceName', 
                                            'AlignedStartPos', "AlignedEndPos", 'MAPQ', 'flag',
                                            'LargestIns', 'LargestDel', 'LargestSclip', 'cigar',])
    # AlignmentSummary.sort_values(['ReadName', 'SegmentNumber'], inplace=True)
    if len(AlignmentSummary) == 0: return 
    
    AlignmentSummary.to_csv(f'{DIR_temp}/temp.{process_name}.tsv', sep='\t', index=False)
    return 