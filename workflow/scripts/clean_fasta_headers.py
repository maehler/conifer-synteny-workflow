import re

class TranscriptNameError(Exception):
    pass

name_regex = re.compile(snakemake.params['name_regex'])

with open(snakemake.input[0]) as f, open(snakemake.output[0], 'w') as of:
    for line in f:
        if not line.startswith('>'):
            of.write(line)
            continue
        name = name_regex.search(line)
        if name is None:
            raise TranscriptNameError(
                f'name matching "{name_regex.pattern}" not found in header: "{line.strip()}"'
            )
        of.write(f'>{name[0]}\n')
