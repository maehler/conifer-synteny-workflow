import click

class BedParsingError(Exception):
    pass

def parse_bed(filename):
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            if len(line) < 4:
                raise BedParsingError('too few columns')

            seq = line[0]
            name = line[3]

            try:
                start = int(line[1])
                end = int(line[2])
            except ValueError:
                raise BedParsingError('start and end must be integers')

@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('bed1', type=click.File())
@click.argument('bed2', type=click.File())
@click.option('--homology', type=click.File())
def main():
    pass

if __name__ == '__main__':
    main()
