rule create_bam_list_all_samples:
    """
    Create text file with paths to BAMS for all samples
    """
    input:
        BAM_PATH
    output:
        f'{PROGRAM_RESOURCE_DIR}/bam_lists/all_samples_bams.list'
    run:
        glue_bams = [bam for bam in glob.glob(f"{input}/final/*.bam") if not os.path.basename(bam).startswith('s_')]
        tor_bams = [bam for bam in glob.glob(f"{input}/toronto_bams/*.bam")]

        all_bams = glue_bams + tor_bams
        with open(output[0], 'w') as f:
            for bam in all_bams:
                sample = os.path.basename(bam).split('_merged')[0]
                if sample in SAMPLES:
                    f.write('{0}\n'.format(bam))


rule create_bam_list_by_city_and_population:
    """
    Create text file with paths to BAMS for each population in each city 
    """
    input:
        bams = rules.create_bam_list_all_samples.output
    output:
        f'{PROGRAM_RESOURCE_DIR}/bam_lists/by_city/{{city}}/{{city}}_{{population}}_bams.list'
    params:
        samples = config['samples']
    run:
        df = pd.read_table(params.samples, sep = '\t')
        df_sub = df[(df['city'] == wildcards.city) & (df['pop'] == int(wildcards.population))]
        samples_city_population = df_sub['sample'].tolist()
        bams = open(input[0], 'r').readlines()
        with open(output[0], 'w') as f:
            for bam in bams:
                sample = os.path.basename(bam).split('_merged')[0]
                if sample in samples_city_population:
                    f.write('{0}'.format(bam))

rule angsd_snps_allSamples:
    """
    Identify SNPs across all samples using ANGSD
    """
    input:
        bams = rules.create_bam_list_all_samples.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        gls = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.beagle.gz',
        mafs = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.mafs.gz',
        snp_stats = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.snpStat.gz',
        hwe = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.hwe.gz',
        pos = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.pos.gz',
        counts = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps.counts.gz'
    log: f"{LOG_DIR}/angsd_snps_allSamples/{{chrom}}_angsd_snps.log"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/snps/allSamples/{{chrom}}/{PREFIX}_{{chrom}}_allSamples_snps'
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 24000,
        runtime = 1440
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doGlf 2 \
            -baq 2 \
            -ref {input.ref} \
            -doCounts 1 \
            -dumpCounts 3 \
            -minQ 30 \
            -minMapQ 30 \
            -remove_bads 1 \
            -skipTriallelic 1 \
            -uniqueOnly 1 \
            -only_proper_pairs 1 \
            -dosnpstat 1 \
            -doHWE 1 \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule create_sites_file:
    """
    Create position file of ngsParalog
    """
    input:
        rules.angsd_snps_allSamples.output.pos
    output:
        f"{PROGRAM_RESOURCE_DIR}/angsd_sites/{PREFIX}_{{chrom}}.sites"
    shell:
        """
        zcat {input} | tail -n +2 | cut -f1,2 > {output}
        """

rule index_snps:
    """
    Index SNPs sites files for ANGSD
    """
    input:
        sites = rules.create_sites_file.output 
    output:
        idx = f"{PROGRAM_RESOURCE_DIR}/angsd_sites/{PREFIX}_{{chrom}}.sites.idx",
        bin = f"{PROGRAM_RESOURCE_DIR}/angsd_sites/{PREFIX}_{{chrom}}.sites.bin"
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    shell:
        """
        angsd sites index {input}
        """

rule angsd_alleleCounts_freq_byPopulation:
    input:
        bams = rules.create_bam_list_by_city_and_population.output,
        sites = rules.create_sites_file.output,
        sites_idx = rules.index_snps.output,
        ref = rules.copy_ref.output,
        ref_idx = rules.samtools_index_ref.output
    output:
        mafs = f'{ANGSD_DIR}/snps/{{city}}/{{population}}/{PREFIX}_{{city}}_{{population}}_{{chrom}}_snps.mafs.gz',
        pos = f'{ANGSD_DIR}/snps/{{city}}/{{population}}/{PREFIX}_{{city}}_{{population}}_{{chrom}}_snps.pos.gz',
        counts = f'{ANGSD_DIR}/snps/{{city}}/{{population}}/{PREFIX}_{{city}}_{{population}}_{{chrom}}_snps.counts.gz',
    log: f'{LOG_DIR}/angsd_alleleCounts_freqs_byPopulation/{{city}}_{{population}}_{{chrom}}_snps.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.938'
    params:
        out = f'{ANGSD_DIR}/snps/{{city}}/{{population}}/{PREFIX}_{{city}}_{{population}}_{{chrom}}_snps'
    threads: 2
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        runtime = lambda wildcards, attempt: attempt * 60
    shell:
        """
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -baq 2 \
            -ref {input.ref} \
            -doCounts 1 \
            -dumpCounts 3 \
            -doMaf 1 \
            -minQ 20 \
            -minMapQ 30 \
            -sites {input.sites} \
            -anc {input.ref} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

ANGSD_GEA_TARGET_FILES = []
for city, pop in zip(CITIES, POPULATIONS):
    for chrom in CHROMOSOMES:
        mafs = f'{ANGSD_DIR}/snps/{city}/{pop}/{PREFIX}_{city}_{pop}_{chrom}_snps.mafs.gz'
        ANGSD_GEA_TARGET_FILES.append(mafs)
        pos = f'{ANGSD_DIR}/snps/{city}/{pop}/{PREFIX}_{city}_{pop}_{chrom}_snps.pos.gz'
        ANGSD_GEA_TARGET_FILES.append(pos)
        counts = f'{ANGSD_DIR}/snps/{city}/{pop}/{PREFIX}_{city}_{pop}_{chrom}_snps.counts.gz'
        ANGSD_GEA_TARGET_FILES.append(counts)

rule angsd_gea_allele_frequencies_done:
    input:
        ANGSD_GEA_TARGET_FILES
    output:
        f"{ANGSD_DIR}/angsd_gea_allele_frequencies.done"
    shell:
        """
        touch {output}
        """
