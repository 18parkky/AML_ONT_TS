import os, time, argparse, subprocess, logging, pathlib
import config, utility

script_dir = pathlib.Path(__file__).resolve().parent
def main():
    parser = argparse.ArgumentParser(description="Process basecalled ONT data")

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
    minimap2_idx    = args['minimap2_idx']
    threads         = args['threads']
    DIR_OUT         = args['DIR_OUT']

    # Create output directories
    if os.path.exists(DIR_OUT)==False: os.mkdir(DIR_OUT)
    if os.path.exists(f'{DIR_OUT}/{SAMPLE}')==False: os.mkdir(f'{DIR_OUT}/{SAMPLE}')
    else:
        logging.error(f'{DIR_OUT}/{SAMPLE} already exists!')
        raise FileExistsError

    # Create log file
    PATH_log = f"{DIR_OUT}/{SAMPLE}/RunPipeline.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'\t{k}\t:\t{v}')
        
    # Run 1_Alignment.py
    cmd = f'python {script_dir}/1_Alignment.py -d {DIR_MINKNOW_OUT} -s {SAMPLE} -b {LIST_BARCODES} -t {threads} -mmi {minimap2_idx} -o {DIR_OUT}'
    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f'Running 1_Alignment.py (elapsed time: {round(elapsed_time, 2)}s)')
    logging.info(cmd)
    subprocess.call(cmd, shell=True)
    
    # Run 2_CallSNV.py
    DIR_ALIGNMENT_OUT = f'{DIR_OUT}/{SAMPLE}/1_minimap2_out'
    cmd = f'python {script_dir}/2_CallSNV.py -d {DIR_ALIGNMENT_OUT} -s {SAMPLE} -b {LIST_BARCODES} -t {threads} -o {DIR_OUT}'
    logging.info(f'Running 2_CallSNV.py (elapsed time: {round(elapsed_time, 2)}s)')
    logging.info(cmd)
    subprocess.call(cmd, shell=True)

    # Run 3_GeneFusion.py
    cmd = f'python {script_dir}/3_GeneFusion.py -d {DIR_MINKNOW_OUT} -s {SAMPLE} -b {LIST_BARCODES} -t {threads} -o {DIR_OUT}'
    logging.info(f'Running 3_GeneFusion.py (elapsed time: {round(elapsed_time, 2)}s)')
    logging.info(cmd)
    subprocess.call(cmd, shell=True)

    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f"Finished pipeline (Total elapsed time: {round(elapsed_time, 2)}s)")

if __name__ == "__main__":
    main()
    