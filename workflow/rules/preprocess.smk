rule gff_to_bed:
    input: lambda w: config['gff'][w.species]
    output: 'results/bed/{species}.bed'
    log: 'log/{species}_gff_to_bed.log'
    params:
        feature=config['gff_to_bed']['feature'],
        name_attr=config['gff_to_bed']['name_attr']
    script: '../scripts/gff_to_bed.py'

checkpoint clean_fasta_headers:
    input: lambda w: config['cds'][w.species]
    output: 'results/cds/{species}_cds_clean_headers.fasta'
    log: 'log/{species}_clean_fasta_headers.log'
    params:
        name_regex=lambda w: config['clean_fasta_headers']['name_regex'][w.species]
    script: '../scripts/clean_fasta_headers.py'

rule link_bed:
    input: 'results/bed/{species}.bed'
    output: 'results/synteny/{species}.bed'
    log: 'log/{species}_link_bed.log'
    shell:
        """
        cd $(dirname {output})
        ln -s ../../{input} $(basename {output})
        """

rule link_cds:
    input: get_clean_cds
    output: 'results/synteny/{species}.cds'
    log: 'log/{species}_link_cds.log'
    shell:
        """
        cd $(dirname {output})
        ln -s ../../{input} $(basename {output})
        """
