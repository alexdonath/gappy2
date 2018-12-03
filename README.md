# gappy2

gappy2 identifies splids (split-inducing indels) from multiple sequence alignments.

gappy2 is the successor of [gappy v1](https://github.com/alexdonath/gappy/). It has been completely reimplemented for various reasons: 1) The main algorithm has been redesigned using a new data structure (interval trees). 2) The functionality offered by BioPython saves a lot of additional work. 3) I wanted to learn Python ;-)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine.

### Prerequisites

gappy2 needs **Python 3.6.6** or later (earlier Python 3.x versions might work but haven't been tested), **BioPython version 1.72** or later (earlier versions >= 1.69 might work but haven't been tested), and the **intervaltree library** by Chaim-Leib Halbert and Konstantin Tretyakov.

The easiest way to set up Python is to install the Anaconda Python distribution:

* [Anaconda](https://www.anaconda.com/download/)

   Download the [Python 3 installer from anaconda.com](https://www.anaconda.com/download/), then open it and follow the prompts. Anaconda Python distribution 5.3.0 or higher should work.

   If you don't want to install Anaconda, you can install Python 3 using your preferred method. In this case, make sure you install the following two dependencies as well.

* [BioPython 1.72](https://www.biopython.org)

   The BioPython package can be obtained from [anaconda.com](http://docs.anaconda.com/anaconda/packages/pkg-docs/). Simply use the terminal or an anaconda prompt and type:

   ```bash
   conda install biopython
   ```

   This should install the latest BioPython version.

* [intervaltree 2.1.0](https://github.com/chaimleib/intervaltree)

   The intervaltree package is distributed by various organizations through the [anaconda cloud](https://anaconda.org/search?q=intervaltree). I chose to use the one distributed by [conda-forge](https://conda-forge.org/) but others should work as well. Again, use the terminal or an anaconda prompt and type:

   ```bash
   conda install -c conda-forge intervaltree
   ```

   This should install the latest version of the intervaltree package.

### Installing

Download the latest version of gappy2 from [github.com/alexdonath/gappy2](https://github.com/alexdonath/gappy2) and move it to wherever you want to install gappy2.

## Usage

Go to the location where you have saved gappy2 and execute

```bash
python3 ./gappy2.py -h
```

Or use the full/relative path to where it's located, e.g.

```bash
python3 ~/bin/gappy2.py -h
```

This will show all available options and parameters:

```txt
usage: gappy2.py [-h] [-f FASTA] [-m MAF] [-o OUTPUT_PREFIX] [-l MINLENGTH]
                 [-u MAXLENGTH] [-s LENGTH] [-c CHAR] [-z] [-v]

gappy2 extracts splids (split-inducing indels) from multiple sequence alignments.

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        input alignment file in FASTA format
  -m MAF, --maf MAF     input alignment file in MAF format
  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        prefix for the output files [GAPPY2]
  -l MINLENGTH, --minlength MINLENGTH
                        minimum indel size [2]
  -u MAXLENGTH, --maxlength MAXLENGTH
                        maximum indel size [inf]
  -s LENGTH, --length LENGTH
                        search for indels of exactly size [not set]
  -c CHAR, --unknownchar CHAR
                        character of unknown sites [?]
  -z, --fuzzy           use fuzzy search [False]
  -v, --version         output version information and exit
```

### Input

The only required input is an alignment in `FASTA` (`-f|--fasta`) or `MAF` (`-m|--maf`) format. The latter can contain a series of multiple alignments in a single file.

### Parameter selection

By default, gappy extracts all splids >= 2 bp and <= the length of the alignment. This is a reasonable choice. Smaller splids have a higher chance to appear by chance. You can adjust these settings by choosing a larger minimum (`-l|--minlength`) and/or a lower upper boundary for the splid length (`-u|--maxlength`).

Alternatively, only splids that have a fixed size (`-s|--length`) can be identified. For example,

```bash
python3 ./gappy2.py -l 10 -u 10 -f input.fasta
```

is equivalent to

```bash
python3 ./gappy2.py -s 10 -f input.fasta
```

Splids are defined as insertions/deletions (indels) that define a bipartition of the sequence set. According to this definition, overlap of a splid with another indel is not allowed. However, sometimes it might be desirable to ignore single residue indels because they appear more often by chance. gappy2 offers a fuzzy search mode (`-z|--fuzzy`) for this.

During a multiple sequence alignment, gaps can be incorporated at the beginning and end of a sequence. This is the case, _e.g._, if incomplete sequenced genes are aligned. In splid analyses it might be necessary to mask these gaps because they provide no valid information about the presence/absence of residues in these sequence(s). It is up to the user to mask these gaps using a specific character (`-c|--unknownchar`) that indicates the difference between a true missing residue (a deletion, `-`) and residues for which no information is available (_e.g._, `?`; default in gappy2). Masking needs to be done by the user before running gappy2.

### Output

Two output files will be created. The first file contains the binary presence/absence coding of the identified splids in `relaxed phylip` format. The second one is a tab-separated file with a detailed overview of the location and length of the splids.

The former can store more than one alignment, which is desired if your input is, _e.g._, a whole genome alignment in `MAF` format.

The output of the latter is in `BED` style, _i.e._, counting starts at 0, the starting position is included, and the ending position is not included.

The output files will be named according to the input file with the prefix `GAPPY2`. Furthermore, a suffix will be appended that indicates the minimum and maximum splid size and the output file type. For example, if your input alignment is called `alignment.fasta` and splids >= 2 bp and <= 10 bp are extracted from the alignment, the result will be written to files called `GAPPY2_alignment.fasta_2-10.tsv` and `GAPPY2_alignment.fasta_2-10.phy`.

Application of the fuzzy search mode will be indicated by an additional `_z`.

### Example

An example alignment is provided in the `example/` directory which looks like this:

```txt
>sequence_1
AAAAA--AAAAAAA-AAAAAAAA-----A
>sequence_2
AAAAAAAA-AAAAAAAA---AAAAAAAAA
>sequence_3
AAAAAAAA-----A-AAAAAAAA-----A
>sequence_4
AAA--AAA-----AAAA---AAAAAAAAA
>sequence_5
AAAAA--AAA-AAA-AAA----AA-AAAA
>sequence_6
A----AAAAAAAAAAAAA----AAAA-AA
>sequence_7
AAAAAAAAAAAAAAAA?????????????
>sequence_8
A-------------AAAAAAAAAAAA-AA
```

Running

```bash
python3 ./gappy2.py -f example.fasta
```

will create a file called `GAPPY2_example.fasta_2-inf.phy` that contains:

```bash
 8 2
sequence_1  10
sequence_2  00
sequence_3  01
sequence_4  01
sequence_5  10
sequence_6  00
sequence_7  00
sequence_8  00

```

and a file `GAPPY2_example.fasta_2-inf.tsv` with the following detailed information about the location of each splid:

```bash
#aln_no start   stop    sequences       length
1       5       7       sequence_1,sequence_5   2
1       8       13      sequence_3,sequence_4   5
```

Likewise,

```bash
./gappy2.py -f example.fasta -z
```

will create a file called `GAPPY2_example.fasta_2-inf_z.phy` that contains:

```bash
 8 3
sequence_1  101
sequence_2  000
sequence_3  011
sequence_4  010
sequence_5  100
sequence_6  000
sequence_7  00?
sequence_8  000
```

and a file `GAPPY2_example.fasta_2-inf_z.tsv`:

```bash
#aln_no start   stop    sequences       length
1       5       7       sequence_1,sequence_5   2
1       8       13      sequence_3,sequence_4   5
1       23      28      sequence_1,sequence_3   5
```

Note that the `?` appears in the third site of `sequence_7`. This is because the sequence contains unknown characters at the sites where the splid in `sequences_1` and `sequence_3` has been identified.

## Benchmarking/Resources

Splid identification in a single multiple alignment with 36 sequences and ca. 2.5 million sites took 1 hour, 15 minutes and needed ca. 1.2 GB of RAM.
Checking 3,846 alignments in a single `MAF` file took 4.5 minutes and 90 MB of RAM.
Both tests were run on a Intel(R) Core(TM) i7-4770 CPU @ 3.40GHz with 8 GB RAM.

## Release History

* **Version 0.1.0**

  Initial release of gappy2.

## Citation

If you find gappy2 useful for your research, please consider citing it:

Donath A., Stadler P.F. 2018. Split-inducing indels in phylogenomic analysis. *Algorithms for Molecular Biology*. 13:12. doi:[10.1186/s13015-018-0130-7](https://doi.org/10.1186/s13015-018-0130-7).

## License

This project is licensed under the GPLv3 License - see the [LICENSE.txt](LICENSE.txt) file for details.