import os, time, glob, math, argparse, logging, subprocess, multiprocessing
import pysam, pandas as pd, numpy as np
import utility, config

def count_lines(filepath):
    return int(subprocess.check_output(['wc', '-l', filepath]).split()[0])

def segment_FASTQ_mainframe(PATH_FASTQ, segment_length_bp, barcode, threads, DIR_segmented_FASTQ):
    nlines = count_lines(PATH_FASTQ)
    raw_chunk_size = nlines // threads
    nlines_per_chunk = (raw_chunk_size // 4) * 4  # ensure it's a multiple of 4
    is_input_gzipped = PATH_FASTQ.endswith('.gz')

    cmd = f'gunzip -c {PATH_FASTQ} | split -l {nlines_per_chunk} --numeric-suffixes=1 --suffix-length=3 --additional-suffix=.fq - {DIR_segmented_FASTQ}/{barcode}_chunk_'
    subprocess.call(cmd, shell=True)

    # Segment reads using multiprocessing
    processes = list()
    for idx, PATH_FASTQ_chunk in enumerate(sorted(glob.glob(f'{DIR_segmented_FASTQ}/{barcode}*.fq'))):
        p = multiprocessing.Process( target=segmentFQ, 
                                    args=[PATH_FASTQ_chunk, f'{barcode}_chunk_{idx+1}', DIR_segmented_FASTQ, segment_length_bp] )
        p.start()
        processes.append(p)
        
    for process in processes:
        process.join()

def segmentFQ( PATH_FASTQ, filename, DIR_OUT, segment_length_bp ):
    with pysam.FastqFile(PATH_FASTQ) as OrigFq:

        PATH_segmented_fq = f'{DIR_OUT}/{filename}.segmented.fq'
        SegmentedFq = open( PATH_segmented_fq, mode='w' )
        
        for entry in OrigFq:
            readname = entry.name
            sequence = entry.sequence
            qualities = entry.quality
            if len(sequence) < 2*segment_length_bp: continue # Don't write reads that can't produce more than 1 segments
            
            for idx, i in enumerate(range(0, len(sequence), segment_length_bp)):
                segmented_sequence  = sequence[i:i+segment_length_bp]
                segmented_qualities = qualities[i:i+segment_length_bp]
                if len(segmented_sequence) < segment_length_bp: continue 
                
                SegmentedFq.write(f'@{readname}_segment{idx+1}\n')
                SegmentedFq.write(f'{segmented_sequence}\n')
                SegmentedFq.write('+\n')
                SegmentedFq.write(f'{segmented_qualities}\n')
                
        SegmentedFq.close()
        # os.remove(PATH_FASTQ)
        return
        
def main():
    parser = argparse.ArgumentParser(description="Gene fusion calling using in-house scripts")

    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    # Required arguments
    required.add_argument('-d', '--DIR_MINKNOW_OUT',      
                        help='Directory of demulitplexed, basecalled outputs', 
                        required=True,
                        )
    required.add_argument('-s', '--SAMPLE',      
                        help='Name of the sample', 
                        required=True,
                        )
    required.add_argument('-b', '--LIST_BARCODES',      
                        help='Comma-separated list of barcodes, e.g., barcode01,barcode02,barcode03', 
                        required=True,
                        )

    # Optional arguments
    optional.add_argument('-l', '--segment_length_bp',      
                        help="Length of the resulting segments. Lower values increase sensitivity but increases runtime (default: 400)", 
                        required=False, type=int, default=400,
                        )
    optional.add_argument('-mmi', '--minimap2_idx',      
                        help=f'PATH to minimap2 index (default: {config.mm2_idx})', 
                        required=False, type=str,
                        default=config.mm2_idx,
                        ) 
    optional.add_argument('--gzip',  
                        help="Apply gzip compression the output segmented FASTQ file (default: False)", 
                        action='store_true'
                        ) 
    optional.add_argument('-t', '--threads',
                        help=f'Number of threads to use (default: {config.threads})',
                        required=False, type=int,
                        default=config.threads,
                        )
    optional.add_argument('-o', '--DIR_OUT',
                        help='Directory to write output files to (default: working directory)',
                        required=False, default=os.getcwd(),
                        )
    
    start_time = time.time()
    args = vars(parser.parse_args())

    DIR_MINKNOW_OUT     = args['DIR_MINKNOW_OUT']
    SAMPLE              = args['SAMPLE']
    LIST_BARCODES       = args['LIST_BARCODES']
    minimap2_idx        = args['minimap2_idx']
    segment_length_bp   = args['segment_length_bp']
    threads             = args['threads']
    DIR_OUT             = args['DIR_OUT']

    DIR_OUT = f'{DIR_OUT}/{SAMPLE}'

    # Create output directories
    if os.path.exists(DIR_OUT)==False: os.mkdir(DIR_OUT)
    DIR_GeneFusion_out = f'{DIR_OUT}/3_GeneFusion'
    if os.path.exists(DIR_GeneFusion_out)==False: os.mkdir(DIR_GeneFusion_out)

    # Create log file
    PATH_log = f"{DIR_OUT}/3_GeneFusion.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'\t{k}\t:\t{v}')

    # Find FASTQ file paths for each barcode
    list_DIR_FASTQ_OUT  = sorted([ d for d in glob.glob(f'{DIR_MINKNOW_OUT}/barcode*') if os.path.basename(d) in LIST_BARCODES ])
    LIST_BARCODES       = sorted([ e.strip() for e in LIST_BARCODES.split(',') ])
    
    # Check presence of FASTQ files before runninng
    for DIR_FASTQ_out in list_DIR_FASTQ_OUT:
        barcode = os.path.basename(DIR_FASTQ_out)
        PATH_FASTQ = glob.glob(f'{DIR_FASTQ_out}/*.fastq.gz')
        if len(PATH_FASTQ)==0:
            logging.error(f'FASTQ file not found for {barcode}!')
            raise FileNotFoundError
        elif len(PATH_FASTQ)==1:
            logging.info(f'FASTQ file found for {barcode}')
        else:
            logging.warning(f'Multiple FASTQ files found for {barcode}! Merging into one.')
            PATH_FASTQ = f'{DIR_FASTQ_out}/{barcode}.fastq.gz'
            cmd = f'cat *fastq.gz > {PATH_FASTQ}'
            subprocess.call(cmd, shell=True)

    # Segment reads
    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f'Segmenting reads (elapsed time: {round(elapsed_time, 2)}s)')
    for DIR_FASTQ_out in list_DIR_FASTQ_OUT:
        barcode = os.path.basename(DIR_FASTQ_out)
        DIR_segmented_FASTQ = f'{DIR_GeneFusion_out}/{barcode}'
        if os.path.exists(DIR_segmented_FASTQ)==False: os.mkdir(DIR_segmented_FASTQ)
        
        PATH_FASTQ = glob.glob(f'{DIR_FASTQ_out}/*.fastq.gz')[0]
        segment_FASTQ_mainframe(PATH_FASTQ, segment_length_bp, barcode, threads, DIR_segmented_FASTQ)
        
        # Merge segmented reads into a single file
        cmd = f'ls {DIR_segmented_FASTQ}/{barcode}_chunk_*segmented.fq* | sort | xargs cat >> {DIR_segmented_FASTQ}/{barcode}.segmented.fq'
        subprocess.call(cmd, shell=True)
        # Delete chunks # f'{barcode}_chunk_{idx+1}',
        for PATH_FASTQ_chunk in glob.glob(f'{DIR_segmented_FASTQ}/{barcode}_chunk*.fq*'):
            os.remove(PATH_FASTQ_chunk)
    
    # Align segments to genome
    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f'Aligning segments to genome (elapsed time: {round(elapsed_time, 2)}s)')
    for barcode in LIST_BARCODES:
        PATH_segmented_FASTQ = f'{DIR_GeneFusion_out}/{barcode}/{barcode}.segmented.fq'
        cmd = f'minimap2 -a {minimap2_idx} {PATH_segmented_FASTQ} -t {threads} | samtools view -Sb | samtools sort -@ {threads} > {DIR_GeneFusion_out}/{barcode}/{barcode}.segmented.bam'
        subprocess.call(cmd, shell=True)
        cmd = f'samtools index {DIR_GeneFusion_out}/{barcode}/{barcode}.segmented.bam -@ {threads}'
        subprocess.call(cmd, shell=True)
        os.remove(f'{DIR_GeneFusion_out}/{barcode}/{barcode}.segmented.fq')

    # Summarize alignments & Calling gene fusions
    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f'Summarizing alignment results and calling gene fusions (elapsed time: {round(elapsed_time, 2)}s)')
    TargetGenes = pd.read_csv(config.PATH_BED, sep='\t', )
    TargetGenes.columns = ['chrom', 'start', 'end', 'gene']            
    chunkSize  = math.ceil( len(TargetGenes) / threads )
    chunks      = [ TargetGenes[i:i+chunkSize] for i in range(0, len(TargetGenes), chunkSize) ]

    if threads > len(TargetGenes):
        logging.warning(f"Number of target genes ({len(TargetGenes)}) is less than the number of threads to use. {len(TargetGenes)} threads will be used.")
        threads = len(TargetGenes)
        chunkSize  = math.ceil( len(TargetGenes) / threads )
        chunks      = [ TargetGenes[i:i+chunkSize] for i in range(0, len(TargetGenes), chunkSize) ]

    for barcode in LIST_BARCODES:
        elapsed_time = utility.elapsedTime( start_time )
        logging.info(f'Processing {barcode} (elapsed time: {round(elapsed_time, 2)}s)')
        PATH_segmented_BAM = f'{DIR_GeneFusion_out}/{barcode}/{barcode}.segmented.bam'
        processes = list()
        for idx, chunk in enumerate(chunks):
            p = multiprocessing.Process( target=utility.SummarizeAlignment_targetGenes, 
                                        args=[PATH_segmented_BAM, chunk, idx, f'{DIR_GeneFusion_out}/{barcode}'] )
            p.start()
            processes.append(p)

            for process in processes: process.join()

        AlignmentSummary_merged = list()
        for PATH_temp in glob.glob(f'{DIR_GeneFusion_out}/{barcode}/temp.*.tsv'):
            AlignmentSummary_merged.append( pd.read_csv(PATH_temp, sep='\t') )
            os.remove( PATH_temp )

        AlignmentSummary_merged = pd.concat(AlignmentSummary_merged)

        # Remove duplicate entries (Needs to be fixed?)
        AlignmentSummary_merged['ReadSegment_ID'] = [ f'{tup.ReadName}_{tup.SegmentNumber}' for tup in AlignmentSummary_merged.itertuples() ]
        duplicates = set()
        for ReadSegment_ID, edf in AlignmentSummary_merged.groupby('ReadSegment_ID'):
            if len(edf)>1: duplicates.add(ReadSegment_ID)
        AlignmentSummary_merged = AlignmentSummary_merged[~(AlignmentSummary_merged['ReadSegment_ID'].isin(duplicates))]
        # AlignmentSummary_merged.to_csv(f'{DIR_GeneFusion_out}/{barcode}/{barcode}.AlignmentSummary.tsv', sep='\t', index=False)
        
        # Calling fusions
        AlignmentSummary_hq = AlignmentSummary_merged[(AlignmentSummary_merged['LargestDel']<=config.maximum_insertion) & 
                                                      (AlignmentSummary_merged['LargestIns']<=config.maximum_deletion) & 
                                                      (AlignmentSummary_merged['LargestSclip']<=config.maximum_softclip) & 
                                                      (AlignmentSummary_merged['MAPQ']>=config.minimum_MAPQ)].copy()
        utility.AnnotateAlignedRegion(AlignmentSummary_hq, TargetGenes, ) 
        AlignmentSummary_hq.dropna().to_csv(f'{DIR_GeneFusion_out}/{barcode}/{barcode}.AlignmentSummary_annotated.tsv', sep='\t', index=False)

        GeneFusionDetections_ReadLevel = list() 

        for ReadName, edf in AlignmentSummary_hq.groupby('ReadName', observed=True):
            edf_valid = edf[(edf['AlignedRegion']!='N/A')]

            if len(set(edf_valid['AlignedRegion']))==2:
                gene1 = edf_valid.iloc[0].AlignedRegion

                for tup in edf_valid.itertuples():
                    gene = tup.AlignedRegion
                    
                    if gene != gene1:
                        gene2 = gene
                        gene2_segment = tup.SegmentNumber
                        breakpoint2 = tup.AlignedStartPos
                        break 

                    gene1_segment   = tup.SegmentNumber
                    breakpoint1     = tup.AlignedEndPos
                    
                GeneFusionDetections_ReadLevel.append( [ReadName, '::'.join( sorted([gene1, gene2]) ), gene1, gene2, 
                                                        breakpoint1, breakpoint2] )
        GeneFusionDetections_ReadLevel = pd.DataFrame(GeneFusionDetections_ReadLevel, columns=['ReadName', 'GeneFusion', 'Gene1', 'Gene2', 'Breakpoint1', 'Breakpoint2',]) # 'n_genes'

        GeneFusionDetections = list()
        known_fusions = pd.read_csv(config.PATH_GF_WHITELIST, sep='\t')
        list_known_fusions = list( known_fusions['fusion'] )
        for GeneFusion, edf in GeneFusionDetections_ReadLevel.groupby('GeneFusion', observed=True):
            if len(edf)==1: continue
            if GeneFusion in list_known_fusions: known = True 
            else: known = False
            GeneFusionDetections.append( [GeneFusion, known, len(edf), int(np.std( edf['Breakpoint1'] )), int(np.std( edf['Breakpoint2'] ))] )
            
        GeneFusionDetections = pd.DataFrame(GeneFusionDetections, columns=['GeneFusion', 'known', 'n', 'std1' ,'std2'])

        GeneFusionDetections.sort_values(['known', 'n'], ascending=False, inplace=True)
        GeneFusionDetections.reset_index(inplace=True, drop=True)

        GeneFusionDetections_ReadLevel.to_csv(f'{DIR_GeneFusion_out}/{barcode}/{barcode}.GeneFusion_ReadLevel.tsv', sep='\t', index=False)
        GeneFusionDetections.to_csv(f'{DIR_GeneFusion_out}/{barcode}/{barcode}.GeneFusion.tsv', sep='\t', index=False)

        os.remove(f'{DIR_GeneFusion_out}/{barcode}/{barcode}.segmented.bam')
        os.remove(f'{DIR_GeneFusion_out}/{barcode}/{barcode}.segmented.bam.bai')

if __name__ == "__main__":
    main()