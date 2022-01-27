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
