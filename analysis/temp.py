def CallSNV(PATH_BAM, PATH_BED, PATH_REF, threads, start_time, sample_name, dict_VariantIdentifier_to_popAF, 
            TargetGenes, PATH_JAVA, PATH_SNPEFF_JAR, PATH_SNPSIFT_JAR, PATH_DBSNP_VCF, database, DIR_OUT, 
            minimum_cov, elapsed_time ):
  
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
    run_ClairS_TO( PATH_BAM, DIR_ClairSTO_out, PATH_BED, PATH_REF, threads, start_time)

    ### Sniffles2 ###
    PATH_Sniffles2_VCF_OUT = f'{DIR_Sniffles2_out}/{sample_name}.sniffles.vcf.gz'
    run_Sniffles2( PATH_BAM, PATH_REF, PATH_Sniffles2_VCF_OUT, threads, start_time, )

    ### bcftools ###
    PATH_bcftools_call_VCF_OUT = f'{DIR_bcftools_out}/{sample_name}.bcftools.vcf.gz'
    run_bcftools( PATH_BAM, PATH_BED, PATH_REF, PATH_bcftools_call_VCF_OUT, threads, start_time )

    ### Merge VCF files ###
    PATH_vcf_merged = f'{DIR_OUT}/{sample_name}_vcf_merged.vcf.gz'
    merge_VCFs(PATH_vcf_merged, DIR_ClairSTO_out, PATH_Sniffles2_VCF_OUT, PATH_bcftools_call_VCF_OUT, DIR_OUT, sample_name, start_time)

    ### Annotate with SnpEff and SnpSift ###
    annotate_VCF(PATH_vcf_merged, PATH_JAVA, PATH_SNPEFF_JAR, PATH_SNPSIFT_JAR, PATH_DBSNP_VCF, database, DIR_SnpEff_out, DIR_SnpSift_out, sample_name, start_time)

    PATH_SnpSift_out = f'{DIR_SnpSift_out}/{sample_name}_vcf_merged.snpEff.snpSift.vcf.gz' 
    vcf_annotated = utility.load_vcf(PATH_SnpSift_out)
    
    ### Filtering out variants with population frequency > threshold ###
    vcf_annotated['Identifier'] = [ f'{list(tup)[1]}:{tup.POS}:{tup.REF}:{tup.ALT}' for tup in vcf_annotated.itertuples() ]
    vcf_annotated = filter_variants_by_pop_freq( vcf_annotated, dict_VariantIdentifier_to_popAF, 0.01 )
    vcf_annotated = append_INFO_entries( vcf_annotated, TargetGenes )
    
    vcf_annotated = vcf_annotated[(vcf_annotated['DP']>=minimum_cov)].copy()
    vcf_annotated.reset_index(inplace=True, drop=True)
    
    PATH_FinalTotal_TSV = f'{DIR_OUT}/{sample_name}_annotated_variants_total.tsv.gz'
    vcf_annotated.to_csv(PATH_FinalTotal_TSV, sep='\t', index=False, compression='gzip')

    if len(vcf_annotated)==0:
        logging.warning('Patient does not have any variants')
        logging.info(f"Finished SNV/SV calling pipeline (Total elapsed time: {round(elapsed_time, 2)}s)")
        return 

    # (vcf_annotated['LOF']!='NotAvailable')
    # Apply different variant condition for different tool
    # ClairS-TO
    condition0 = (vcf_annotated['tool'].isin( ['ClairS-TO_SNV', 'ClairS-TO_indel'] ))
    condition1 = ( (vcf_annotated['gnomAD_AF']<=0.01) & (vcf_annotated['GeneName'].isin( list(TargetGenes['gene']) )))
    condition2 = ( ~(vcf_annotated['AaChange'].isin(['', 'NotAvailable'])) )
    condition3 = ( ~(vcf_annotated['CLNSIG'].isin(['Benign', 'NotAvailable'])) | vcf_annotated['VarImpact'].str.contains(r'\bHIGH\b', case=True, na=False))
    vcf_annotated_final_1 = vcf_annotated[ condition0 & condition1 & condition2 & condition3 ].sort_values('SampleID')[config.relevant_columns].reset_index(drop=True)

    # Sniffles2
    condition0 = (vcf_annotated['tool'].isin( ['Sniffles2',] ))
    condition1 = ( (vcf_annotated['gnomAD_AF']<=0.01)  & (vcf_annotated['GeneName'].isin( list(TargetGenes['gene']) )))
    condition2 = ( ~(vcf_annotated['AaChange'].isin(['', 'NotAvailable'])) )
    vcf_annotated_final_2 = vcf_annotated[ condition0 & condition1 & condition2 ].sort_values('SampleID')[config.relevant_columns].reset_index(drop=True)

    # bcftools
    condition0 = (vcf_annotated['tool'].isin( ['bcftools'] ))
    condition1 = (  (vcf_annotated['gnomAD_AF']<=0.01)  &
                    (vcf_annotated['AF']>=0.025) &
                    (vcf_annotated['GeneName'].isin( list(TargetGenes['gene']) )))
    condition2 = ( ~(vcf_annotated['AaChange'].isin(['', 'NotAvailable'])) & (vcf_annotated['CLNSIG'].isin(['Likely pathogenic', 'Pathogenic', 'Pathogenic/Likely pathogenic'])) )
    vcf_annotated_final_3 = vcf_annotated[ condition0 & condition1 & condition2 ].sort_values('SampleID')[config.relevant_columns].reset_index(drop=True)
        
    vcf_annotated_final = pd.concat([vcf_annotated_final_1, vcf_annotated_final_2, vcf_annotated_final_3])
    # vcf_annotated_final = vcf_annotated[ condition1 & condition2 ].sort_values('SampleID')[relevant_columns].reset_index(drop=True)
    vcf_annotated_final.sort_values('SampleID', inplace=True)
    vcf_annotated_final.reset_index(inplace=True, drop=True)

    vcf_annotated_final['AF'] = [float(AF)*100 for AF in vcf_annotated_final['AF']]
    vcf_annotated_final.reset_index(inplace=True, drop=True)
    PATH_Final_TSV = f'{DIR_OUT}/{sample_name}_annotated_variants_final.tsv'
    vcf_annotated_final.to_csv(PATH_Final_TSV, sep='\t', index=False,)
    return 