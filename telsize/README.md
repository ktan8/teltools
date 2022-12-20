# TelSize
TelSize is a package for estimating telomere length from long-read sequencing data.

The TelSize package consists of two pipeline. The first pipeline extracts telomeric long-reads from sequencing data for deeper analysis. The second pipeline analyzes these telomeric reads in greater detail and provides telomeric length estimates.

An example of telomere length estimate from single long-reads is depicted in the image below:

![Examples](telsize/img/telomere_length_estimate_longreads.png)

## Requirements
- Perl 5
- Python 3
- numpy
- matplotlib
- minimap2

## Extract telomeric reads
This tool extracts candidate telomeric long-reads for deeper analysis

```
Usage: perl extractTelReads/extract_telomeric_reads.pl <input_fasta> <output_label>
```

A compressed fasta file with the name `<output_label>.telomeric.fasta.gz` will be generated at the end.


## Estimate length of extracted telomeric reads

Example command:
```
python estimateTelLength/estimateTelSize.py --noseq --nofig ./example/sample.telomeric_longreads.fasta.gz > test.telsize.txt
```

A full description of options and arguments are as follow:

```
usage: estimateTelSize.py [-h] [--cutoff cutoff] [--movave movave]
                          [--movmed movmed] [--motif motif] [--format format]
                          [--folder folder] [--noseq] [--nofig]
                          [--threads threads] [--penalty]
                          [--penaltyval movave]
                          fasta

Get length for repeats in fasta

positional arguments:
  fasta                fasta file to evaluate (unzipped|gzipped)

optional arguments:
  -h, --help           show this help message and exit
  --cutoff cutoff      cutoff of repeat signal (default: 0.35)
  --movave movave      Window size for moving average window (default: 50)
  --movmed movmed      Window size for moving median window (default: 501)
  --motif motif        motif to assess for repeat signal (e.g. TTAGGG)
                       (default: TTAGGG)
  --format format      output file type (e.g. pdf, png) (default: pdf)
  --folder folder      folder for images (default: ./img_1671479887/)
  --noseq              do not print sequence of read in output (default:
                       False)
  --nofig              do not plot figures for sequences (default: True)
  --threads threads    number of threads to use (default: 12)
  --penalty            calculate telomere length using penalty based method
                       (default: False)
  --penaltyval movave  penalty value to use for penalty method (default: 0.5)
```

