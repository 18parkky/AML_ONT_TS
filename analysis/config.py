### Parameters ### 
# PATH_SNPEFF_JAR='/home/nanopore/programs/snpEff/snpEff.jar'
# PATH_SNPSIFT_JAR='/home/nanopore/programs/snpEff/SnpSift.jar'
# PATH_REF='/home/nanopore/Documents/AML/reference/hg38.fa'
# PATH_DBSNP_VCF='/home/nanopore/Documents/AML/reference/GCF_000001405.40.ROI_SnpSiftCompatible.vcf'
# PATH_BED='/home/nanopore/Documents/AML/reference/TargetGenesAndFusions.bed'
# mm2_idx='/home/nanopore/Documents/AML/reference/hg38_map_ont.mmi'
# gnomAD='/home/nanopore/Documents/AML/reference/Variant_to_AF.pickle'

DIR_reference       ='/home/18parkky/workspace/AML_project/260416/AML/reference'
PATH_SNPEFF_JAR     =f'{DIR_reference}/snpEff/snpEff.jar'
PATH_SNPSIFT_JAR    =f'{DIR_reference}/snpEff/SnpSift.jar'
PATH_REF            =f'{DIR_reference}/hg38.fa'
PATH_DBSNP_VCF      =f'{DIR_reference}/GCF_000001405.40.ROI_SnpSiftCompatible.vcf'
PATH_BED            =f'{DIR_reference}/TargetGenesAndFusions.bed'
PATH_GF_WHITELIST   =f'{DIR_reference}/known_genefusions.tsv'
mm2_idx             =f'{DIR_reference}/hg38_map_ont.mmi'
gnomAD              =f'{DIR_reference}/Variant_to_AF.pickle'
PATH_GeneFusion_prev_results=f'{DIR_reference}/GeneFusion_260602.tsv'

minimum_coverage=10
PATH_JAVA='java'
SnpEff_db='GRCh38.mane.1.5.refseq'
population_af=0.01
threads=18

maximum_insertion=10
maximum_deletion=10
maximum_softclip=10
minimum_MAPQ=60

DIR_JAFFAL_reference = f'{DIR_reference}/JAFFAL_ref'

relevant_columns = ['SampleID', 'GeneName', 'MANE', 'AaChange', 'CLNHGVS', 'CLNSIG', 'VarType', 'VarImpact',  'RS', 'LOF', 'AF', 'DP', 'ALLELEID', 'Identifier', 'tool', ]

ClairS_TO_param = [
    '--platform ont_r10_dorado_hac_4khz',
    '--snv_min_af 0.01',
    '--indel_min_af 0.01',
]

Sniffles2_param = [
    '--minsvlen 10',
    '--mosaic',
    '--mosaic-af-max 0.5',
    '--mosaic-af-min 0.01',
    '--mosaic-qc-invdup-min-length 10',
    '--mosaic-include-germline',
    '--cluster-binsize 1',
]

bcftools_mpileup_param = [
    '-d 99999', # max_depth
    '--min-BQ 30',
    '--max-BQ 60',
    '--skip-indels',
]

bcftools_call_param = [
    '-m',
    '-Oz',
    '--keep-alts',
]

dbSNP_Number_to_CLNSIG = {
    '0' : 'Uncertain significance',
    '1' : 'Not provided',
    '2' : 'Benign',
    '3' : 'Likely benign',
    '4' : 'Likely pathogenic',
    '5' : 'Pathogenic',
    '6' : 'Drug response',

    '8' : 'Confers sensitivity',
    '9' : 'Risk factor', 
    '10' : 'Association',
    '11' : 'Protective',

    '12' : 'Conflicting interpretations of pathogenicity',
    '13' : 'Affects',
    '14' : 'Association not found',
    '15' : 'Benign/Likely benign',
    '16' : 'Pathogenic/Likely pathogenic',
    '17' : 'Conflicting data from submitters',
    '18' : 'Pathogenic, low penetrance',
    '19' : 'Likely pathogenic, low penetrance',
    '20' : 'Established risk allele',

    '21' : 'Likely risk allele',
    '22' : 'Uncertain risk allele',

    '225' : 'other',
    
    'NotAvailable' : 'NotAvailable',
}