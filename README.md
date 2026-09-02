# AKHAL: Assembly Graph Analysis Tool

## Overview

`akhal` is a toolkit for processing and analyzing r/GFA (Graphical Fragment Assembly) files. 
It validates graphs, reports statistics, and converts GAF (Graph Alignment Format) to SAM (Sequence Alignment Map).

It is structured like `htslib` + `samtools`: a reusable library (`libakhal`) that owns the file formats and data structures, and a thin command-line tool (`akhal`) built on top of it. 

This README documents the command-line tool. 
The library has its own reference under [`docs/`](docs/README.md) - one page per module, one section per public function, each with a short snippet showing how to call it and the data preparation it needs. 
Start at the [library index](docs/README.md), or go straight to a module: [`gfa`](docs/gfa.md), [`call`](docs/call.md), [`annot`](docs/annot.md), [`vg`](docs/vg.md), [`fasta`](docs/fasta.md), [`gaf`](docs/gaf.md), [`sam`](docs/sam.md), [`vcf`](docs/vcf.md), [`io`](docs/io.md), [`kstr`](docs/kstr.md), [`util`](docs/util.md), [`error`](docs/error.md).

## Installation

Build the library and the CLI with:

```sh
make
```

This produces `libakhal.a` and the `akhal` executable.

Requirements: a C compiler and `zlib` (`-lz`, used by `vg2gfa` to decompress `.vg` files). No protobuf or other external libraries are needed.

## Usage

`akhal` provides the following commands:

```sh
./akhal <PROGRAM> [...ARGS]
```

### Commands

#### 1. `parse`
Validates an r/GFA file and ensures its correctness. It checks segments and links and makes sure that everything is consistent. It also checks the overlapings, if it is presented.

**Usage:**
```sh
./akhal parse <r/GFA file>
```

#### 2. `stats`
Computes and outputs statistics about an r/GFA file.

**Usage:**
```sh
./akhal stats <r/GFA file>
```

The statistics include:
- **Segment count**: Number of segments in the graph.
- **Segment avg length**: Average segment length.
- **Segment std length**: Standard deviation of segment lengths.
- **Segment min length**: Minimum segment length.
- **Segment max length**: Maximum segment length.
- **Link count**: Number of links between segments.
- **Link overlapping avg length**: Average length of overlapping links.
- **Link overlapping std length**: Standard deviation of link overlap lengths.
- **Minimum in degree**: Minimum number of incoming links.
- **Maximum in degree**: Maximum number of incoming links.
- **Minimum out degree**: Minimum number of outgoing links.
- **Maximum out degree**: Maximum number of outgoing links.

#### 3. `extract`
Extract information from the r/GFA file. The first argument picks what to pull out.

**Usage:**
```sh
./akhal extract fa   <r/GFA file> <OUTPUT .fa/.fasta file>
./akhal extract path <r/GFA file> <OUTPUT .fa/.fasta file> [path name]
./akhal extract vcf  <r/GFA file> <OUTPUT .vcf file> [--ref <name>] [--fasta <FASTA file>]
```

##### `fa`
Writes one FASTA record per `P` line, exactly as the graph stores them.

##### `path`
Same output, but fragmented paths are stitched back together first.
Builders like `vg` split one reference across consecutive `P` lines (`chr22[0]`, `chr22[21]`, ...), so asking for `chr22` collects every fragment whose name reduces to it - vg's `name[start]` and `name:start-end` decorations are stripped, and a PanSN name such as `GRCh38#0#chr22[0]` is found by its bare contig name too.
The fragments are then ordered by the offset in their name and chained through the `L` lines: one fragment follows another when a link joins its last segment to the other's first with matching orientations, the other has no predecessor yet, and the join does not close a cycle.
Chaining only follows the forward strand, so a fragment stored reverse-complemented relative to its neighbours stays on its own.

Everything is written, merged or not: a chain of several fragments takes their shared base name (`chr22`, or `chr22_1`, `chr22_2`, ... when one base yields more than one chain), and a fragment that nothing joined keeps the name it already had. 
With no path name given, every path in the graph is grouped by its own base name.

##### `vcf`
Reads the graph as a reference plus its differences.
One walk is declared the backbone - a path by default (`--ref <name>`, else the graph's first path; fragmented `P` lines are stitched together exactly as `extract path` does), or with `--fasta` an external sequence traced through the nodes, `--ref` then naming the record to trace.
Every node on that walk gets a reference coordinate, and everything the backbone does not cover is by definition a variation: a detour that leaves the backbone at one node and rejoins it at a later one spells an alternate allele over the reference span between them, and a link that skips forward along the backbone itself spells a deletion of the bases it jumps over.
Detours that leave and rejoin at the same nodes share a `REF` span, so they come out as one multi-allelic row.

The result is VCF 4.2, with the backbone name as `CHROM`, `INFO` carrying `END` (1-based inclusive end of the `REF` allele) and `TYPE` (one of `SNP`/`MNP`/`INS`/`DEL`/`COMPLEX` per `ALT`).
Alleles with an empty side carry the preceding reference base on both sides, as VCF requires.
Coordinates assume blunt joins (0 bp link overlaps), which is what graph builders produce; a backbone carrying real overlaps is reported.

Example:
```sh
$ ./akhal extract vcf graph.gfa out.vcf --ref chr22
$ grep -v '^##' out.vcf
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr22	11	.	C	G,T	.	.	END=11;TYPE=SNP,SNP
chr22	21	.	GTTTT	G	.	.	END=25;TYPE=DEL
chr22	35	.	A	AAAA	.	.	END=35;TYPE=INS
```

#### 4. `sort`
Topologically sorts a graph and renumbers its nodes `1..N` in the sorted order. Ordering uses Kahn's algorithm; ties in the ready set (the "hops") are broken alphabetically by node **sequence content**, so the result is independent of the input's node numbering — two graphs identical in topology and content sort the same way. Nodes in cycles, if any, are appended after the acyclic prefix. The sorted graph is re-emitted with the new ids (S/L/P lines remapped, link orientation and overlap preserved).

**Usage:**
```sh
./akhal sort <input .gfa file> [output .gfa file]
```

Note: If no output file is given, the sorted GFA is written to standard output.

#### 5. `rank`
Rewrites a graph's `SR:i:` ranks against its backbone. 
In rGFA, rank 0 is the reference and anything higher came from a sample; `rank` decides which segments are which and re-emits the graph.

By default the backbone is the graph's own `P` lines: every segment any path visits becomes rank 0, everything else rank 1. 
Fragmented paths are consolidated first - a reference that arrived as `chr22[0]`, `chr22[21]`, ... leaves as a single `P chr22`, chained through the links exactly as `extract path` does. 
A graph with no `P` lines has no backbone to rank against, so the command stops rather than silently flattening an rGFA's own tags; use `--fasta` to give it one.

With `--fasta`, an external reference becomes the backbone instead: the sequence is traced through the graph, the nodes it walks become rank 0 and every other node rank 1, and the old `P` lines are replaced by one named after the FASTA record spelling that walk. 
This is how you re-root a graph on a different reference. 
`--ref <name>` picks which record to trace when the FASTA holds several.

**Usage:**
```sh
./akhal rank <input .gfa file> [output .gfa file] [--fasta <FASTA file>] [--ref <sequence name>]
```

Note: If no output file is given, the graph is written to standard output.

Note on reading: a file's own `SR:i:` tags are always authoritative and are never overwritten on load. 
The reader only derives ranks when the file carried none, so a plain GFA with `P` lines comes back with a rank-0 backbone while a real rGFA keeps exactly the ranks it shipped with. 
`rank` is the explicit way to overwrite them.

Example:
```sh
$ ./akhal rank graph.gfa --fasta GRCh38.chr22.fa
[INFO] (call) backbone 'chr22': 50818468 bp over 41233 node(s)
[INFO] backbone taken from GRCh38.chr22.fa: 41233 node(s) at rank 0, 18904 at rank 1
```

#### 6. `vg2gfa`
Converts vg's native `.vg` format (gzip/BGZF-compressed Protobuf) to GFA. The `.vg` parser is hand-written pure C - no protobuf or libvgio needed, only zlib.

**Usage:**
```sh
./akhal vg2gfa <input .vg file> [output .gfa file]
```

Note: If no output file is given, the GFA is written to standard output.

#### 7. `gaf2sam`
Converts a GAF file to a SAM file.

**Usage:**
```sh
./akhal gaf2sam <r/GFA file> <GAF file>  <FASTA file> <OUTPUT file>  [--simple]
```

Note: The reads should be stored in FASTA format and provided to the program. 
GAF file does not store sequences, hence, reads are needed when converting to SAM.

Note: `simple` option is optional. 
If provided, CIGAR string matches `(=)` and mismatches `(X)` will be replaced with sequence match `(M)`.

#### 8. `sampoke`
Validate SAM file (converted from gaf). It takes reference file and SAM to process CIGAR strings. 
Optionally, it can print the filtered SAM file that contains valid lines.

**Usage:**
```sh
./akhal sampoke <FASTA file> <GAF file> <OUTPUT file>
```

Note: Output file here is optional.

#### 9. `annotate`
Traces the origin of every graph node and saves the result as a binary `.annot` file. 
The backbone reference path (a `P` line) is identified first - by `--ref <name>`, or defaulting to the graph's first path - and its nodes are marked `backbone` with their reference coordinates. 
With `--vcf`, each variant is matched to the alternative side of its bubble: the shared REF/ALT prefix is stripped, the branch point is located on the backbone, and the nodes spelling the alternate allele are annotated like `SNP chr1:12345 A>G rs99`. 
Pure deletions produce only an edge, so they have nothing to annotate. 
With `--fasta`, each sequence is walked through the graph (link overlaps honoured) and every visited non-backbone node is annotated `SEQ <name> <offset>` with its origin and position in that sequence. 
Both inputs are optional; nodes that no input explains are simply left out and later report as `unknown`. 
A node explained more than once accumulates its annotations separated by `; `.

**Usage:**
```sh
./akhal annotate <r/GFA file> <OUTPUT .annot file> [--vcf <VCF file>] [--fasta <FASTA file>] [--ref <path name>]
```

#### 10. `annotget`
Looks up node annotations in a `.annot` file - the graph itself is not needed. 
Each queried node prints one line, `<id> <status> [<info>]`, where status is `backbone`, `annot`, or `unknown`. 
With no node ids, every record in the file is dumped.

**Usage:**
```sh
./akhal annotget <.annot file> [node id ...]
```

Example:
```sh
$ ./akhal annotget graph.annot 4 6 2 999
4	annot	SNP chr1:5 A>C rs1; SEQ sample1 4
6	annot	INS chr1:9 C>CTT
2	backbone	REF chr1 4-5
999	unknown
```

## License
is released under the BSD 3-Clause License, which allows for redistribution and use in source and binary forms, with or without modification, under certain conditions. 
For more detailed terms, please refer to the [license file](https://github.com/akmami/akhal/blob/main/LICENCE).

## Author
Developed by Akmuhammet Ashyralyyev.

**Note from the author**:

This tool is named after one of the most elegant horses, the [Akhal-Teke](https://en.wikipedia.org/wiki/Akhal-Teke). 
This breed is one of the oldest domesticated animals and is considered one of the most beautiful and intelligent horses in the world.
