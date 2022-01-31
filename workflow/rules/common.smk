def get_species():
    return [
        config['species']['source']['name'],
        config['species']['target']['name']
    ]

def get_clean_cds(wildcards):
    regex = config['clean_fasta_headers']['name_regex'][wildcards.species]
    if regex is None or len(regex) == 0:
        return config['cds'][wildcards.species]
    else:
        return 'results/cds/{species}_cds_clean_headers.fasta'

def get_genome_fasta(species):
    return config['genome'][species]

def get_genome_fasta_index(species):
    return '{}.fai'.format(get_genome_fasta(species))
