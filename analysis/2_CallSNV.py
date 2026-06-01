import os, time, argparse, subprocess, logging
import utility, config
import pandas as pd

def run_ClairS_TO( PATH_BAM, DIR_ClairSTO_out, threads, start_time):
    ''' 
    description: Run ClairS-TO using subprocess
    input: Various parameters for ClairS-TO
    output: None
    '''
    logging.info(f"Starting ClairS-TO",)
    uid = os.getuid()
    gid = os.getgid()
    
    cmd = f'docker run --user {uid}:{gid} ' \
        f'-v {PATH_BAM}:{PATH_BAM}:ro ' \
        f'-v {PATH_BAM}.bai:{PATH_BAM}.bai:ro ' \
        f'-v {DIR_ClairSTO_out}:{DIR_ClairSTO_out}:rw ' \
        f'-v {config.PATH_BED}:{config.PATH_BED}:ro ' \
        f'-v {config.PATH_REF}:{config.PATH_REF}:ro ' \
        f'-v {config.PATH_REF}.fai:{config.PATH_REF}.fai:ro ' \
        f'hkubal/clairs-to:latest /opt/bin/run_clairs_to ' \
        f'--tumor_bam_fn {PATH_BAM} --ref_fn {config.PATH_REF} ' \
        f'--threads {threads} --output_dir {DIR_ClairSTO_out} ' \
        f'--bed_fn {config.PATH_BED} ' \
        f'{" ".join(config.ClairS_TO_param)}'
    logging.info(cmd,)
    subprocess.call(cmd, shell=True)
    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f"Finished ClairS-TO (elapsed time: {round(elapsed_time, 2)}s)",)
    return None

def run_Sniffles2( PATH_BAM, PATH_Sniffles2_VCF_OUT, threads, start_time, ):
    ''' 
    description: Run Sniffles2 using subprocess
    input: Various parameters for Sniffles2
    output: None
    '''
    logging.info(f"Starting Sniffles2",)        

    cmd = f'sniffles -i {PATH_BAM} ' \
        f'-v {PATH_Sniffles2_VCF_OUT} ' \
        f'--reference {config.PATH_REF} ' \
        f'-t {threads} ' \
        f'{" ".join(config.Sniffles2_param)}'
    logging.info(cmd,)
    subprocess.call(cmd, shell=True)
    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f"Finished Sniffles2 (elapsed time: {round(elapsed_time, 2)}s)",)
    return None
    
def run_bcftools( PATH_BAM, PATH_bcftools_call_VCF_OUT, threads, start_time ):
    ''' 
    description: Run bcftools mpileup & call using subprocess
    input: Various parameters for bcftools
    output: None
    '''
    logging.info(f"Starting bcftools mpileup/call",)
    cmd = f'bcftools mpileup ' \
        f'-f {config.PATH_REF} ' \
        f'--regions-file {config.PATH_BED} ' \
        f'--threads {threads} ' \
        f'{" ".join(config.bcftools_mpileup_param)} ' \
        f'{PATH_BAM} | ' \
        f'bcftools call -o {PATH_bcftools_call_VCF_OUT} {" ".join(config.bcftools_call_param)}'
    logging.info(cmd,)
    subprocess.call(cmd, shell=True),
    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f"Finished bcftools mpileup/call (elapsed time: {round(elapsed_time, 2)}s)",)
    return None

def merge_VCFs(PATH_vcf_merged, DIR_ClairSTO_out, PATH_Sniffles2_VCF_OUT, PATH_bcftools_call_VCF_OUT, sample_name, start_time):
    ''' 
    description: Merge VCF files created from each tool - ClairS-TO, Sniffles2, bcftools - and write merged VCF to disk.
    input: Outout Path/Directory of each tool and output directory to write the merged VCF file to
    output: None
    '''
    logging.info(f"Merging VCF files")
    vcf_columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'VALUE']
    vcf_merged = list()
    
    # 1. ClairS-TO
    PATH_clairsto_snv_vcf = f'{DIR_ClairSTO_out}/snv.vcf.gz'
    PATH_clairsto_snv_indel = f'{DIR_ClairSTO_out}/indel.vcf.gz'
    
    try:
        vcf_clairsto_snv = utility.load_vcf(PATH_clairsto_snv_vcf)
        vcf_clairsto_snv.columns = vcf_columns
        vcf_clairsto_snv['SampleID'] = sample_name
        vcf_clairsto_snv['tool'] = 'ClairS-TO_SNV'
        vcf_merged.append( vcf_clairsto_snv )
    except FileNotFoundError:
        logging.info(f"ClairS-TO did not call any SNVs")

    try:
        vcf_clairsto_indel = utility.load_vcf(PATH_clairsto_snv_indel)
        vcf_clairsto_indel.columns = vcf_columns
        vcf_clairsto_indel['SampleID'] = sample_name
        vcf_clairsto_indel['tool'] = 'ClairS-TO_indel'
        vcf_merged.append( vcf_clairsto_indel )
    except FileNotFoundError:
        logging.info(f"ClairS-TO did not call any SVs")

    # 2. Sniffles2
    try:
        vcf_sniffles2 = utility.load_vcf(PATH_Sniffles2_VCF_OUT)
        vcf_sniffles2.columns = vcf_columns
        vcf_sniffles2['SampleID'] = sample_name
        vcf_sniffles2['tool'] = 'Sniffles2'
        vcf_merged.append(vcf_sniffles2)
    except FileNotFoundError:
        logging.info(f"Sniffles did not call any SVs")

    # 3. bcftools
    try:
        vcf_bcftools = utility.load_vcf(PATH_bcftools_call_VCF_OUT)
        vcf_bcftools.columns = vcf_columns
        vcf_bcftools['SampleID'] = sample_name
        vcf_bcftools['tool'] = 'bcftools'
        vcf_merged.append(vcf_bcftools)
    except FileNotFoundError:
        logging.info(f"bcftools did not call any SNVs")

    if len(vcf_merged)==0:
        logging.warning(f'No variants called by any tool for {sample_name}')
        return
    
    vcf_merged = pd.concat(vcf_merged)
    vcf_merged.reset_index(inplace=True, drop=True)
    
    vcf_merged.to_csv(PATH_vcf_merged, 
                        sep='\t', 
                        compression='gzip', 
                        index=False)
    
    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f"Finished merging VCFs (elapsed time: {round(elapsed_time, 2)}s)")    
    return

def annotate_VCF(PATH_vcf_merged, DIR_SnpEff_out, DIR_SnpSift_out, sample_name, start_time):
    ''' 
    description: Given Path to the merged VCF file, predict & annotate variants using SnpEff and SnpSift.
    input: Path to merged VCF file, binary executables to Java, SnpSift, SnpEff jar files, and other SnpEff/SnpSift parameters.
    output: None
    '''
    logging.info(f"Annotating variants with SnpEff/SnpSift")
    PATH_SnpEff_out = f'{DIR_SnpEff_out}/{sample_name}_vcf_merged.snpEff.vcf.gz'
    PATH_SnpSift_out = f'{DIR_SnpSift_out}/{sample_name}_vcf_merged.snpEff.snpSift.vcf.gz' 

    cmd = f'{config.PATH_JAVA} -jar {config.PATH_SNPEFF_JAR} {config.SnpEff_db} {PATH_vcf_merged} | gzip > {PATH_SnpEff_out}'
    logging.info(cmd)
    subprocess.call(cmd, shell=True),
    
    cmd = f'{config.PATH_JAVA} -jar {config.PATH_SNPSIFT_JAR} annotate {config.PATH_DBSNP_VCF} {PATH_SnpEff_out} | gzip > {PATH_SnpSift_out}'
    logging.info(cmd,)
    subprocess.call(cmd, shell=True),

    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f"Finished variant annotations (elapsed time: {round(elapsed_time, 2)}s)")
    return None

def filter_variants_by_pop_freq( vcf, dict_Variant_to_Pop_AF, population_AF_threshold ):
    ''' 
    description: Given VCF and Dictionary that maps variants to population 
                AF (according to gnomAD), filter out variants with high population AF.
    input: VCF (Pandas DataFrame), variants-to-population dictionary (dict), frequency threshold (float)
    output: filtered VCF (Pandas DataFrame)
    '''
    logging.info(f"Filtering out variants with population AF above {population_AF_threshold}")
    
    vcf['gnomAD_AF'] = [
        dict_Variant_to_Pop_AF.get(identifier, 0)
        for identifier in vcf['Identifier']
    ]
    return vcf[(vcf['gnomAD_AF']<=population_AF_threshold)].copy() 

def convert_CLNSIG_values( clnsig, dbSNP_Number_to_CLNSIG ):
    ''' 
    description: Given CLNSIG annotation value (from dbSNP VCF), convert the dbSNP CLNSIG numerical values to actual interpretable CLNSIG values.
    input: CLNSIG value obtained from SnpSift annotation using dbSNP VCF (string) and numeric-to-string CLNSIG mapper (dictionary)
    output: interpretable CLNSIG value (string)
    '''
    dict_clnsig = dict()
    for e1 in clnsig.split(','):
        for e2 in e1.split('|'):
            try: dict_clnsig[int(e2)] = dbSNP_Number_to_CLNSIG[e2]
            except: pass
    dict_clnsig = dict(sorted(dict_clnsig.items(), key=lambda x : x[0], reverse=True))
    if len(dict_clnsig)==0: return 'NotAvailable'
    for k, v in dict_clnsig.items():
        return v

def parse_info_field(vcf_info_string):
    ''' 
    Parse VCF INFO field into dictionary
    '''
    return dict(
        entry.split('=', 1)
        for entry in vcf_info_string.split(';')
        if '=' in entry
    )

def extract_allele_metrics(vcf):
    ''' 
    Extract AF and DP based on tool type
    '''
    af_list, dp_list = [], []
    
    for tup in vcf.itertuples():
        info_dict = parse_info_field(tup.INFO)
        
        if tup.tool in ['ClairS-TO_SNV', 'ClairS-TO_indel']:
            values = tup.VALUE.split(':')
            dp_list.append(int(values[2]))
            af_list.append(float(values[3]))
            
        elif tup.tool == 'Sniffles2':
            support = info_dict.get('SUPPORT', '0')
            af = info_dict.get('AF', '0')
            dp_list.append(int(support) if support else 0)
            af_list.append(float(af) if af else 0.0)
            
        elif tup.tool == 'bcftools':
            dp4 = list(map(int, info_dict.get('DP4', '0,0,0,0').split(',')))
            total = sum(dp4)
            dp_list.append(total)
            af_list.append((dp4[2] + dp4[3]) / total if total > 0 else 0)
    
    return af_list, dp_list

def extract_dbsnp_annotations(vcf, info_dicts):
    ''' 
    Extract dbSNP-related fields
    '''
    return {
        'ALLELEID': [info_dicts[id].get('ALLELEID', 'NotAvailable') for id in vcf.Identifier],
        'CLNSIG': [convert_CLNSIG_values(info_dicts[id].get('CLNSIG', 'NotAvailable'), 
                                         config.dbSNP_Number_to_CLNSIG) for id in vcf.Identifier],
        'CLNHGVS': [info_dicts[id].get('CLNHGVS', 'NotAvailable') for id in vcf.Identifier],
        'CLNVI': [info_dicts[id].get('CLNVI', 'NotAvailable') for id in vcf.Identifier],
        'RS': [info_dicts[id].get('RS', 'NotAvailable') for id in vcf.Identifier],
    }

def extract_snpeff_annotations(vcf, info_dicts, target_genes):
    ''' 
    Extract SnpEff ANN field annotations
    '''
    target_gene_set = set(target_genes['gene'])
    
    annotations = {
        'VarType': [], 'GeneName': [], 'VarImpact': [],
        'MANE': [], 'AaChange': [], 'LOF': []
    }
    
    for identifier in vcf.Identifier:
        ann = info_dicts[identifier].get('ANN', 'NotAvailable')
        
        if ann == 'NotAvailable':
            for key in annotations:
                annotations[key].append('NotAvailable')
            continue
        
        var_types, gene_names, impacts, mane, aa_changes = set(), set(), set(), set(), set()
        
        for entry in ann.split(','):
            fields = entry.split('|')
            if len(fields) < 11:  # ANN should have at least 11 fields
                continue
            
            var_types.add(fields[1])
            impacts.add(fields[2])
            mane.add(fields[6])
            if fields[3] in target_gene_set:
                gene_names.add(fields[3])
            if fields[10]:
                aa_changes.add(fields[10])
        
        annotations['VarType'].append(','.join(var_types))
        annotations['GeneName'].append(','.join(gene_names))
        annotations['VarImpact'].append(','.join(impacts))
        annotations['MANE'].append(','.join(mane))
        annotations['AaChange'].append(','.join(aa_changes))
        annotations['LOF'].append(info_dicts[identifier].get('LOF', 'NotAvailable'))
    
    return annotations

def append_INFO_entries(vcf, target_genes):
    ''' 
    Enrich VCF with parsed annotations from INFO field
    '''
    # Parse INFO fields once
    info_dicts = {
        tup.Identifier: parse_info_field(tup.INFO)
        for tup in vcf.itertuples()
    }
    
    # Extract metrics
    vcf['AF'], vcf['DP'] = extract_allele_metrics(vcf)
    
    # Extract dbSNP annotations
    for key, values in extract_dbsnp_annotations(vcf, info_dicts).items():
        vcf[key] = values
    
    # Extract SnpEff annotations
    for key, values in extract_snpeff_annotations(vcf, info_dicts, target_genes).items():
        vcf[key] = values
    
    return vcf

def append_INFO_entries_depr( vcf, TargetGenes,):
    ''' 
    description: Given VCF append several INFO entires as columns to the VCF for each INFO entry of interest
    input: VCF (Pandas DataFrame)
    output: VCF (Pandas DataFrame)
    '''
    dict_VariantIdentifier_to_INFOdict = dict()
    for tup in vcf.itertuples():
        info_dict = dict(
            entry.split('=', 1)
            for entry in tup.INFO.split(';')
            if '=' in entry
        )
        dict_VariantIdentifier_to_INFOdict[tup.Identifier] = info_dict
    
    # 1. Allele frequency (AF) and Depth (DP)
    temp_AF, temp_DP = list(), list()
    for tup in vcf.itertuples():
        if tup.tool == 'ClairS-TO_SNV' or tup.tool == 'ClairS-TO_indel':
            temp_DP.append( int(tup.VALUE.split(':')[2]) )
            temp_AF.append( float(tup.VALUE.split(':')[3]) )
                
        elif tup.tool == 'Sniffles2':
            info_dict = dict(
                entry.split('=', 1)
                for entry in tup.INFO.split(';')
                if '=' in entry
            )
            try: temp_DP.append( int(info_dict.get('SUPPORT')) )
            except TypeError: temp_DP.append( 0 ) 
            try: temp_AF.append( float(info_dict.get('AF')) )
            except TypeError: temp_AF.append( 0 ) 
        
        elif tup.tool == 'bcftools':
            info_dict = dict(item.split('=') for item in tup.INFO.split(';') if '=' in item)
            dp4 = list(map(int, info_dict.get('DP4', '0,0,0,0').split(',')))
            num_ref = dp4[0] + dp4[1]
            num_alt = dp4[2] + dp4[3]
            total = num_ref + num_alt      
            temp_DP.append( total )
            try: temp_AF.append( num_alt/total ) 
            except ZeroDivisionError: temp_AF.append( 0 )
    vcf['AF'] = temp_AF
    vcf['DP'] = temp_DP
    
    # 2. ALLELEID, CLNSIG, CLNHGVS, CLNVI, RS 
    vcf['ALLELEID'] = [ dict_VariantIdentifier_to_INFOdict[Identifier].get('ALLELEID', 'NotAvailable') for Identifier in vcf.Identifier ]
    vcf['CLNSIG']   = [ convert_CLNSIG_values(dict_VariantIdentifier_to_INFOdict[Identifier].get('CLNSIG', 'NotAvailable'), config.dbSNP_Number_to_CLNSIG) for Identifier in vcf.Identifier ]
    vcf['CLNHGVS']  = [ dict_VariantIdentifier_to_INFOdict[Identifier].get('CLNHGVS', 'NotAvailable') for Identifier in vcf.Identifier ]
    vcf['CLNVI']    = [ dict_VariantIdentifier_to_INFOdict[Identifier].get('CLNVI', 'NotAvailable') for Identifier in vcf.Identifier ]
    vcf['RS']   = [ dict_VariantIdentifier_to_INFOdict[Identifier].get('RS', 'NotAvailable') for Identifier in vcf.Identifier ]
    
    Column_Mane_transcript, Column_GeneName, Column_VarImpact, Column_VarType, Column_AaChange, Column_LOF = list(), list(), list(), list(), list(), list()

    for tup in vcf.itertuples():
        
        annot_entries = dict_VariantIdentifier_to_INFOdict[tup.Identifier].get('ANN','NotAvailable')
        
        # If 'ANN' field is missing or 'NotAvailable', append 'NotAvailable' for both columns and continue
        if annot_entries == 'NotAvailable':
            Column_Mane_transcript.append('NotAvailable')
            Column_GeneName.append('NotAvailable')
            Column_VarImpact.append('NotAvailable')
            Column_AaChange.append('NotAvailable')
            Column_VarType.append('NotAvailable')
            Column_LOF.append('NotAvailable')
            continue 
        
        # Process the annotation entries
        ManeTranscript = set()
        VarType = set()
        GeneName = set()
        VarImpact = set()
        AaChange = set()
        for entry in annot_entries.split(','):
            fields = entry.split('|')
            VarType.add(fields[1])   # Variant type is at index 1
            VarImpact.add(fields[2])  # Impact is at index 2
            ManeTranscript.add(fields[6])   # MANE transcript is at index 6
            if fields[3] in list(TargetGenes['gene']):
                GeneName.add(fields[3])  # Gene name is at index 3
            if fields[10]!='':
                AaChange.add(fields[10])  # AaChange is at index 10
        
        # Append comma-separated values for each gene and variant impact
        Column_VarType.append(','.join(VarType))
        Column_GeneName.append(','.join(GeneName))
        Column_Mane_transcript.append(','.join(ManeTranscript))
        Column_VarImpact.append(','.join(VarImpact))
        Column_AaChange.append(','.join(AaChange))

        # Collect LOF prediction
        Column_LOF.append( info_dict.get('LOF', 'NotAvailable') )
        
    vcf['VarType'] = Column_VarType
    vcf['GeneName'] = Column_GeneName
    vcf['VarImpact'] = Column_VarImpact
    vcf['MANE'] = Column_Mane_transcript
    vcf['AaChange'] = Column_AaChange
    vcf['LOF'] = Column_LOF
    
    return vcf

def CallSNV(PATH_BAM, threads, start_time, sample_name, dict_VariantIdentifier_to_popAF, TargetGenes, DIR_OUT ):

    # ClairS-TO configuration
    DIR_ClairSTO_out    = f"{DIR_OUT}/ClairS_TO_{sample_name}"
    if os.path.exists( DIR_ClairSTO_out ) == False: os.mkdir( DIR_ClairSTO_out )

    # Sniffles2 configuration
    DIR_Sniffles2_out = f"{DIR_OUT}/Sniffles2_{sample_name}"
    if os.path.exists( DIR_Sniffles2_out ) == False: os.mkdir( DIR_Sniffles2_out )

    # bcftools configuration
    DIR_bcftools_out = f"{DIR_OUT}/bcftools_{sample_name}"
    if os.path.exists( DIR_bcftools_out ) == False: os.mkdir( DIR_bcftools_out )

    # SnpEff & SnpSift configurations
    DIR_SnpEff_out = f"{DIR_OUT}/SnpEff_{sample_name}"
    DIR_SnpSift_out = f"{DIR_OUT}/SnpSift_{sample_name}"

    if os.path.exists(DIR_SnpEff_out)==False:   os.mkdir(DIR_SnpEff_out)
    if os.path.exists(DIR_SnpSift_out)==False:  os.mkdir(DIR_SnpSift_out)
    
    ### ClairS-TO ###
    run_ClairS_TO( PATH_BAM, DIR_ClairSTO_out, threads, start_time)

    ### Sniffles2 ###
    PATH_Sniffles2_VCF_OUT = f'{DIR_Sniffles2_out}/{sample_name}.sniffles.vcf.gz'
    run_Sniffles2( PATH_BAM, PATH_Sniffles2_VCF_OUT, threads, start_time, )

    ### bcftools ###
    PATH_bcftools_call_VCF_OUT = f'{DIR_bcftools_out}/{sample_name}.bcftools.vcf.gz'
    run_bcftools( PATH_BAM, PATH_bcftools_call_VCF_OUT, threads, start_time )

    ### Merge VCF files ###
    PATH_vcf_merged = f'{DIR_OUT}/{sample_name}_vcf_merged.vcf.gz'
    merge_VCFs(PATH_vcf_merged, DIR_ClairSTO_out, PATH_Sniffles2_VCF_OUT, PATH_bcftools_call_VCF_OUT, sample_name, start_time)
    if not os.path.exists(PATH_vcf_merged):
        logging.warning(f'No merged VCF for {sample_name}. Skipping annotation.')
        return

    ### Annotate with SnpEff and SnpSift ###
    annotate_VCF(PATH_vcf_merged, DIR_SnpEff_out, DIR_SnpSift_out, sample_name, start_time)

    PATH_SnpSift_out = f'{DIR_SnpSift_out}/{sample_name}_vcf_merged.snpEff.snpSift.vcf.gz' 
    vcf_annotated = utility.load_vcf(PATH_SnpSift_out)
    
    ### Filtering out variants with population frequency > threshold ###
    vcf_annotated['Identifier'] = [ f'{list(tup)[1]}:{tup.POS}:{tup.REF}:{tup.ALT}' for tup in vcf_annotated.itertuples() ]
    vcf_annotated = filter_variants_by_pop_freq( vcf_annotated, dict_VariantIdentifier_to_popAF, config.population_af, )
    vcf_annotated = append_INFO_entries( vcf_annotated, TargetGenes )
    
    vcf_annotated = vcf_annotated[(vcf_annotated['DP']>=config.minimum_coverage)].copy()
    vcf_annotated.reset_index(inplace=True, drop=True)
    
    PATH_FinalTotal_TSV = f'{DIR_OUT}/{sample_name}_annotated_variants_total.tsv.gz'
    vcf_annotated.to_csv(PATH_FinalTotal_TSV, sep='\t', index=False, compression='gzip')

    if len(vcf_annotated)==0:
        logging.warning('Patient does not have any variants')
        elapsed_time = utility.elapsedTime( start_time )
        logging.info(f"Finished SNV/SV calling pipeline (Total elapsed time: {round(elapsed_time, 2)}s)")
        return 

    # (vcf_annotated['LOF']!='NotAvailable')
    # Apply different variant condition for different tool
    # ClairS-TO
    target_gene_list = list(TargetGenes['gene'])

    condition0 = (vcf_annotated['tool'].isin( ['ClairS-TO_SNV', 'ClairS-TO_indel'] ))
    condition1 = ( (vcf_annotated['gnomAD_AF']<=config.population_af) & (vcf_annotated['GeneName'].isin( target_gene_list )))
    condition2 = ( ~(vcf_annotated['AaChange'].isin(['', 'NotAvailable'])) )
    condition3 = ( ~(vcf_annotated['CLNSIG'].isin(['Benign', 'NotAvailable'])) | vcf_annotated['VarImpact'].str.contains(r'\bHIGH\b', case=True, na=False))
    vcf_annotated_final_1 = vcf_annotated[ condition0 & condition1 & condition2 & condition3 ].sort_values('SampleID')[config.relevant_columns].reset_index(drop=True)

    # Sniffles2
    condition0 = (vcf_annotated['tool'].isin( ['Sniffles2',] ))
    condition1 = ( (vcf_annotated['gnomAD_AF']<=config.population_af)  & (vcf_annotated['GeneName'].isin( target_gene_list )))
    condition2 = ( ~(vcf_annotated['AaChange'].isin(['', 'NotAvailable'])) )
    vcf_annotated_final_2 = vcf_annotated[ condition0 & condition1 & condition2 ].sort_values('SampleID')[config.relevant_columns].reset_index(drop=True)

    # bcftools
    condition0 = (vcf_annotated['tool'].isin( ['bcftools'] ))
    condition1 = (  (vcf_annotated['gnomAD_AF']<=config.population_af)  &
                    (vcf_annotated['AF']>=0.5) &
                    (vcf_annotated['GeneName'].isin( target_gene_list )))
    condition2 = ( ~(vcf_annotated['AaChange'].isin(['', 'NotAvailable'])) & (vcf_annotated['CLNSIG'].isin(['Likely pathogenic', 'Pathogenic', 'Pathogenic/Likely pathogenic'])) )
    vcf_annotated_final_3 = vcf_annotated[ condition0 & condition1 & condition2 ].sort_values('SampleID')[config.relevant_columns].reset_index(drop=True)
        
    vcf_annotated_final = pd.concat([vcf_annotated_final_1, vcf_annotated_final_2, vcf_annotated_final_3])
    # vcf_annotated_final = vcf_annotated[ condition1 & condition2 ].sort_values('SampleID')[relevant_columns].reset_index(drop=True)
    vcf_annotated_final.sort_values('SampleID', inplace=True)
    vcf_annotated_final.reset_index(inplace=True, drop=True)
    vcf_annotated_final['AF'] = [float(AF)*100 for AF in vcf_annotated_final['AF']]
    PATH_Final_TSV = f'{DIR_OUT}/{sample_name}_annotated_variants_final.tsv'
    vcf_annotated_final.to_csv(PATH_Final_TSV, sep='\t', index=False,)
    return 

def main():
    parser = argparse.ArgumentParser(description="SNV/SV calling pipeline, using bcftools, ClairS-TO, and Sniffles2")

    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    # Required arguments
    required.add_argument('-d', '--DIR_ALIGNMENT_OUT',      
                        help='Directory of alignment outputs', 
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
    optional.add_argument('-o', '--DIR_OUT',
                        help='Directory to write output files to (default: working directory)',
                        required=False, default=os.getcwd(),
                        )
    
    start_time = time.time()
    args = vars(parser.parse_args())

    DIR_ALIGNMENT_OUT   = args['DIR_ALIGNMENT_OUT']
    SAMPLE              = args['SAMPLE']
    LIST_BARCODES       = args['LIST_BARCODES']
    threads             = args['threads']
    DIR_OUT             = args['DIR_OUT']

    DIR_OUT = f'{DIR_OUT}/{SAMPLE}'

    # Create output directories
    if os.path.exists(DIR_OUT)==False: os.mkdir(DIR_OUT)
    DIR_CallSNV_out = f'{DIR_OUT}/2_CallSNV_out'
    if os.path.exists(DIR_CallSNV_out)==False: os.mkdir(DIR_CallSNV_out)

    # Create log file
    PATH_log = f"{DIR_OUT}/2_CallSNV.log"
    logging.basicConfig(filename=PATH_log, level=logging.INFO)
    logging.info(f"Listing inputs:")
    for k, v in args.items():
        logging.info(f'\t{k}\t:\t{v}')

    # Load necessary reference files for variant preprocessing
    TargetGenes = pd.read_csv(config.PATH_BED, sep='\t', header=None)
    TargetGenes.columns = ['chrom', 'start', 'end', 'gene']
    dict_VariantIdentifier_to_popAF = utility.loadFromPickle(config.gnomAD)

    # Find BAM file paths for each barcode
    LIST_BARCODES = sorted([ e.strip() for e in LIST_BARCODES.split(',') ])

    dict_barcode_to_PATH_BAM = dict()
    for barcode in LIST_BARCODES:
        PATH_BAM_e = f'{DIR_ALIGNMENT_OUT}/{barcode}.bam'
        PATH_BAI_e = f'{DIR_ALIGNMENT_OUT}/{barcode}.bam.bai'
        if os.path.exists(PATH_BAM_e)==False:
            logging.error(f'BAM not found for {barcode}! Aborting')
            return
        elif os.path.exists(PATH_BAI_e)==False:
            logging.warning(f'BAI not found for {barcode}. Creating BAI.')
            cmd = f'samtools index {PATH_BAM_e} -@ {threads}'
            subprocess.call(cmd, shell=True)
        
        dict_barcode_to_PATH_BAM[barcode] = PATH_BAM_e

    for barcode, PATH_BAM in dict_barcode_to_PATH_BAM.items():
        logging.info(f'Processing {barcode}',)
        DIR_OUT_s = f'{DIR_CallSNV_out}/{barcode}'
        if os.path.exists(DIR_OUT_s)==False: os.mkdir(DIR_OUT_s)

        CallSNV(PATH_BAM, threads, start_time, barcode, dict_VariantIdentifier_to_popAF, TargetGenes, DIR_OUT_s )

    elapsed_time = utility.elapsedTime( start_time )
    logging.info(f"Finished SNV/SV calling pipeline (Total elapsed time: {round(elapsed_time, 2)}s)")
    

if __name__ == "__main__":
    main()
    