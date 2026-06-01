import os, glob, argparse, subprocess, logging, time
import config, utility

def main():
    parser = argparse.ArgumentParser(description="Align basecalled reads to reference genome using minimap2")

    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    # Required arguments
    required.add_argument('-d', '--DIR_MINKNOW_OUT',      
                        help='Directory of MinKNOW outputs, e.g., ~/minknow/pass', 
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
    optional.add_argument('-t', '--threads',      
                        help=f'Number of threads to use (default: {config.threads})', 
                        required=False, type=int,
                        default=config.threads,
                        ) 
    optional.add_argument('-mmi', '--minimap2_idx',      
                        help=f'PATH to minimap2 index (default: {config.mm2_idx})', 
                        required=False, type=str,
                        default=config.mm2_idx,
                        ) 
    optional.add_argument('-o', '--DIR_OUT',
                        help='Directory to write output files to (default: working directory)', 
                        required=False, default=os.getcwd(),
                        )
    
    start_time = time.time()
    args = vars(parser.parse_args())

    DIR_MINKNOW_OUT = args['DIR_MINKNOW_OUT']
    SAMPLE          = args['SAMPLE']
    LIST_BARCODES   = args['LIST_BARCODES']
    threads         = args['threads']
    minimap2_idx    = args['minimap2_idx']
    DIR_OUT         = args['DIR_OUT']

    DIR_OUT = f'{DIR_OUT}/{SAMPLE}'
    if os.path.exists(DIR_OUT)==False: os.mkdir(DIR_OUT)
    DIR_mm2_out = f'{DIR_OUT}/1_minimap2_out'
    if os.path.exists(DIR_mm2_out)==False: os.mkdir(DIR_mm2_out)

    PATH_log = f"{DIR_OUT}/1_alignment.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'\t{k}\t:\t{v}')

    list_DIR_FASTQ_OUT  = sorted([ d for d in glob.glob(f'{DIR_MINKNOW_OUT}/barcode*') if os.path.basename(d) in LIST_BARCODES ])
    LIST_BARCODES        = sorted([ e.strip() for e in LIST_BARCODES.split(',') ])

    for DIR_FASTQ_OUT_per_barcode in list_DIR_FASTQ_OUT:
        PATH_FASTQ=glob.glob(f'{DIR_FASTQ_OUT_per_barcode}/*.f*q.gz')[0]
        barcode = os.path.basename(DIR_FASTQ_OUT_per_barcode)
        PATH_mm2_out=f'{DIR_mm2_out}/{barcode}.bam'
        if os.path.exists(PATH_mm2_out)==False:
            cmd = f'minimap2 -a {minimap2_idx} -t {threads} {PATH_FASTQ} | samtools view -Sb | samtools sort > {PATH_mm2_out}'
            logging.info(f'Running minimap2 for {barcode}')
            subprocess.call(cmd,shell=True)
        else:
            logging.info(f'BAM already exists for {barcode}')
        PATH_BAM_idx = f'{PATH_mm2_out}.bai'
        if os.path.exists(PATH_BAM_idx)==False:
            logging.info(f'BAI doesn\'t exist for {barcode}, creating BAI')
            cmd = f'samtools index {PATH_mm2_out} -@ {threads}'
            subprocess.call(cmd, shell=True)

    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f"Finished alignment (Total elapsed time: {round(elapsed_time, 2)}s)")

if __name__ == "__main__":
    main()
    