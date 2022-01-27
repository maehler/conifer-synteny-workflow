import re

feature = snakemake.params.feature
name_attr = snakemake.params.name_attr

with open(snakemake.input[0]) as f,
        open(snakemake.output[0], 'w') as of:
    for line in f:
        if line.startswith('#'):
            continue
        line = line.strip().split()
        if line[2] != feature:
            continue
        name = re.findall(r'{}=([^;]+);?'.format(name_attr), line[8])
        assert(len(name) == 1)
        print('\t'.join((line[0], line[3], line[4], name[0])), file=of)
