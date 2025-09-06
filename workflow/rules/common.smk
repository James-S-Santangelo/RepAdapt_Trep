def get_population_saf_files(wildcards):
    all_saf_files = ANGSD_SAF_TARGET_FILES
    pop1 = wildcards.pop_comb.split('_')[0]
    pop2 = wildcards.pop_comb.split('_')[1]
    regex1 = r'(?<=_)' + re.escape(wildcards.city) + '_' + re.escape(pop1) + '_' + re.escape(wildcards.chrom) + r'(?=\.saf\.idx)'
    regex2 = r'(?<=_)' + re.escape(wildcards.city) + '_' + re.escape(pop2) + '_' + re.escape(wildcards.chrom) + r'(?=\.saf\.idx)'
    saf1 = [x for x in all_saf_files if re.search(regex1, os.path.basename(x), re.IGNORECASE)]
    saf2 = [x for x in all_saf_files if re.search(regex2, os.path.basename(x), re.IGNORECASE)]
    return saf1 + saf2

def get_population_saf_files_random100Mb(wildcards):
    all_saf_files = ANGSD_SAF_RANDOM100MB_TARGET_FILES
    pop1 = wildcards.pop_comb.split('_')[0]
    pop2 = wildcards.pop_comb.split('_')[1]
    regex1 = r'^' + re.escape(wildcards.city) + '_' + re.escape(pop1) + '_random100Mb' + r'(?=\.saf\.idx)'
    regex2 = r'^' + re.escape(wildcards.city) + '_' + re.escape(pop2) + '_random100Mb' + r'(?=\.saf\.idx)'
    saf1 = [x for x in all_saf_files if re.search(regex1, os.path.basename(x), re.IGNORECASE)]
    saf2 = [x for x in all_saf_files if re.search(regex2, os.path.basename(x), re.IGNORECASE)]
    return saf1 + saf2
