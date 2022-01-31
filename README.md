# Pine-spruce synteny

This workflow aims to identify synteny between pine and spruce.

## Input

- Genome fasta files (indexed)
- GFF3 files with gene annotations
- CDS fasta files

## Configuration

### CDS files

If the names of the transcripts in the CDS file does not match with the names in the GFF3 file, then there are two options.
In case the correct name cannot be derived from the fasta headers, the headers would have to be manually edited.
If the names are contained in the headers, you can also define a regular expression that will be used for adjusting the headers.
These are defined in `config/config.yaml` under `clean_fasta_headers.name_regex`.
Set these to `null` or an empty string if you do not want to adjust the headers automatically.

The fasta files can contain either nucleotide or protein sequences.
In the config file, you can specify which it is under `cds.type`, and the options are `nucl` for nucleotide sequences, and `prot` for protein sequences.

## Output

The main output right now are plots that visualise the synteny between the two species.
These can be created using using snakemake with the command

```
snakemake --use-conda -c1 results/synteny/plots/<species1>_<seq>_vs_<species2>.png
```

where `<species1>` is the source species, `<seq>` is the sequence in the source species that should be visualised, and `<species2>` is the target species.
Every sequence in the target species that is represented in the synteny will be presented in the visualisation.
