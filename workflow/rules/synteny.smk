rule run_jcvi_ortholog:
    input:
        bed=(
            'results/synteny/{species1}.bed',
            'results/synteny/{species2}.bed'
        ),
        cds=(
            'results/synteny/{species1}.cds',
            'results/synteny/{species2}.cds'
        )
    output: 'results/synteny/{species1}.{species2}.anchors'
    log: 'log/{species1}_{species2}_jcvi.log'
    conda: '../envs/jcvi.yaml'
    threads: 4
    params:
        seq_type=config['cds']['type'],
        species1=lambda w: w.species1,
        species2=lambda w: w.species2
    shell:
        """
        cd results/synteny
        python -m jcvi.compara.catalog ortholog \\
            --cpus={threads} \\
            --dbtype={params.seq_type} \\
            {params.species1} \\
            {params.species2}
        """

rule chrom_dotplot:
    input:
        anchors='results/synteny/{species1}.{species2}.anchors' \
                .format(species1=config['species']['source']['name'], \
                    species2=config['species']['target']['name']),
        source_bed='results/synteny/{species1}.bed',
        target_bed=lambda w: [f'results/synteny/{x}.bed' \
                              for x in get_species() \
                              if x != w.species1],
        source_fai=lambda w: get_genome_fasta_index(w.species1),
        target_fai=lambda w: get_genome_fasta_index(*[x for x in get_species() if x != w.species1])
    output:
        png='results/synteny/plots/{species1}_{seq}_vs_{species2}.png'
    log: 'log/{species1}_{seq}_vs_{species2}_chrom_dotplot.log'
    conda: '../envs/python.yaml'
    script: '../scripts/chrom_dotplot.R'
