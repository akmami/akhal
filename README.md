# AKHAL: Assembly Graph Analysis Tool

## Overview

`akhal` is a toolkit for processing and analyzing r/GFA (Graphical Fragment Assembly) files. 
It validates graphs, reports statistics, and converts GAF (Graph Alignment Format) to SAM (Sequence Alignment Map).

It is structured like `htslib` + `samtools`: a reusable library (`libakhal`) that owns the file formats and data structures, and a thin command-line tool (`akhal`) built on top of it. 

This README documents the command-line tool. 
The library has its own reference under [`docs/`](docs/README.md) - one page per module, one section per public function, each with a short snippet showing how to call it and the data preparation it needs. 
Start at the [library index](docs/README.md), or go straight to a module: [`gfa`](docs/gfa.md), [`call`](docs/call.md), [`diff`](docs/diff.md), [`rgfa`](docs/rgfa.md), [`compact`](docs/compact.md), [`dot`](docs/dot.md), [`annot`](docs/annot.md), [`vg`](docs/vg.md), [`fasta`](docs/fasta.md), [`gaf`](docs/gaf.md), [`sam`](docs/sam.md), [`vcf`](docs/vcf.md), [`io`](docs/io.md), [`kstr`](docs/kstr.md), [`util`](docs/util.md), [`error`](docs/error.md).

## Installation

Build the library and the CLI with:

```sh
make
```

This produces `libakhal.a` and the `akhal` executable.

Requirements: a C compiler and `zlib` (`-lz`, used by `vg2gfa` to decompress `.vg` files). No protobuf or other external libraries are needed.

To install the tool, the library and its headers, and tab completion for bash and zsh:

```sh
sudo make install
```

`PREFIX` defaults to `/usr/local`; override it for a home install that needs no root:

```sh
make install PREFIX=~/.local
```

`PREFIX=.` installs into the build tree itself, which is handy for trying things out - the library and the headers are already where they would be copied, so those two steps are skipped. 
`make uninstall` takes it all back out (and leaves the source tree's own `lib/` and `include/akhal/` alone), and `DESTDIR` is honoured for staged installs.

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
Computes and outputs statistics about an r/GFA graph or a GAF alignment file. 
The extension picks which: `.gfa`/`.rgfa` report on the graph, `.gaf` on the alignments.

**Usage:**
```sh
./akhal stats <r/GFA file> [--no-degrees] [--ranks]
./akhal stats <GAF file> [--cigar]
```

##### Graph statistics
The graph is never loaded: the file is read once and the counts and length
distributions are running totals. Two figures are properties of the whole
file rather than of any line and cost memory to answer, so each has a switch.
Degrees are computed by default: the segment and link ids are kept as flat
arrays, about 4 bytes per S line and 8 per L line, and sorted in place at the
end so the degrees can be read off as run lengths - roughly 8 GB on a graph of
700 million segments and links, whatever the ids look like; `--no-degrees`
skips that and the pass holds nothing. `--ranks` (off by default) also keeps
every P-line step, which on a file that carries all its haplotypes as P lines
dwarfs the graph itself; a file with its own `SR` tags answers the rank
question for free either way.

- **Segment count**: Number of segments in the graph.
- **Rank 0 segment count**: Segments at rank 0, from the file's `SR` tags when it has them, otherwise derived from the P lines under `--ranks` (reported as n/a without it).
- **Segment avg length**: Average segment length.
- **Segment std length**: Standard deviation of segment lengths.
- **Segment min length**: Minimum segment length.
- **Segment max length**: Maximum segment length.
- **Link count**: Number of links between segments.
- **Link overlapping avg length**: Average length of overlapping links.
- **Link overlapping std length**: Standard deviation of link overlap lengths.
- **In degree avg / std**, **Minimum / Maximum in degree**: Incoming links per segment, over segments that have any (omitted under `--no-degrees`).
- **Out degree avg / std**, **Minimum / Maximum out degree**: The same for outgoing links (omitted under `--no-degrees`).

Ids named by an L line (or, under `--ranks`, a P line) that no S line defines are counted and reported as a warning; under `--no-degrees` without `--ranks` nothing is checked.

##### Alignment statistics
Everything below the alignment counts describes *primary* alignments with a mapping quality above 0, which is the population `gaftools stat` reports on as well. 
An alignment is primary unless it carries a `tp:A` tag that is not `P`, so a file without the tag is entirely primary.

- **Total alignments**, with the **Primary** / **Secondary** split.
- **Reads with at least one alignment**: distinct query names among the primary alignments.
- **Total aligned bases**: the sum of column 10 (residue matches).
- **Average highest sequence identity**: per read, the best `matches / block_len` it achieved, averaged over reads - so a read aligned in ten places contributes once.
- **Average highest map ratio**: the same, for `(qend - qstart) / qlen`.
- **Mapping quality**, **Sequence identity**, **Map ratio**, **Query length**, **Alignment block length** and **Path node count**, each as avg / std / min. / max. over the alignments themselves. Path node count is the number of oriented nodes in the path string; a stable path name counts as one.

`--cigar` adds a pass over the `cg:Z` difference CIGARs, reporting the number of deletion, insertion, substitution and match runs, how many of each span 50 bp or more, and how many alignments are a single `=` run end to end. 
Primary alignments with no `cg:Z` tag are counted and reported as a warning rather than skipped silently.

The figures line up with `gaftools stat` on the same file, with two deliberate differences. 
gaftools divides its "Average mapping quality" by the total alignment count while summing only over primaries; `Mapping quality avg` here divides by the primary count, so the two disagree whenever a file has secondaries. 
And gaftools calls an alignment perfect when its CIGAR is a single run of any operation; this counts only a single `=` run.

Example:
```sh
$ ./akhal stats aln.gaf
Total alignments: 501
  Primary: 312
  Secondary: 189
Reads with at least one alignment: 196
Total aligned bases: 1295243
Average highest sequence identity: 0.917408
Average highest map ratio: 0.509446
Mapping quality avg: 69.294872
...
Path node count max.: 12

* Figures cover primary alignments with a mapping quality above 0
```

The whole file is read in one streaming pass and the distributions are accumulated as they arrive, so memory scales with the number of distinct read names rather than with the number of alignments.

Graph statistics are read the same way. 
Nothing the summary reports needs a graph - the counts are running totals, the two distributions are accumulated as they arrive, and the only per-segment state is its two degrees and whether a path visits it - so the file is read once and no graph is built. 
On chr22 that is 248 MB and 5 s instead of 1.5 GB and 27 s. 
Segment ids scattered too far apart to index that way are noticed, and the file is summarized through the graph instead: the same numbers, at the old price.

#### 3. `extract`
Extract information from the r/GFA file. The first argument picks what to pull out.

**Usage:**
```sh
./akhal extract fa   <r/GFA file> <OUTPUT .fa/.fasta file> [wrap length]
./akhal extract path <r/GFA file> <OUTPUT .fa/.fasta file> <path name> [path name ...] [wrap length]
./akhal extract vcf  <r/GFA file> <OUTPUT .vcf file> [--ref <name>] [--fasta <FASTA file>]
```

`wrap length` sets how many bases each output FASTA line carries; it defaults to 80 when omitted. 
For `path`, a trailing numeric argument is read as the wrap length rather than a path name - but only when a name precedes it, so a lone `60` is still the path called `60`.

##### `fa`
Writes every `P` line in the graph as one FASTA record, in file order, named after the path. 
One `P` line is one record: a reference the file splits over several `P` lines comes out as several records, which is what the file actually says.

##### `path`
The same records, but only for the paths you name - one name at least, as many as you like. 
`fa` is how to take every path in the graph; this is how to take a few of them.

Names match the `P` line exactly, except that a bare contig name also finds a PanSN path - `chr22` selects `GRCh38#0#chr22` without your having to spell the whole thing. 
Every path the name selects is written, so if the file carries the name more than once you get a record for each. 
Several names are written to the one file, in the order given, and a name that matches nothing stops the command rather than leaving a half-written file behind.

##### `vcf`
Reads the graph as a reference plus its differences.
One walk is declared the backbone - a `P` line by default (`--ref <name>`, else the graph's first path), or with `--fasta` an external sequence traced through the nodes, `--ref` then naming the record to trace.
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

#### 4. `compare`
Compares two files of the same kind and reports how much of one is present in the other. The first argument picks the kind.

**Usage:**
```sh
./akhal compare gfa <A .gfa file> <B .gfa file> [--verbose]
./akhal compare gaf <A .gaf file> <B .gaf file> [--verbose]
```

The exit status is 0 when the two files agree, 1 when they differ and 2 when the comparison could not be made, so either target can be used as a test.

##### `gfa`
Compares two graphs that need not agree on segment ids.

Nothing is compared by id: builders number nodes as they emit them, and `sort` renumbers them again. 
Segments are matched on their **sequence content** instead - both files are sorted by sequence and walked once, side by side, and every distinct sequence is put on one shared label. 
Links are then relabelled and compared as (from, orientation, to, orientation, overlap) tuples, canonicalized so that the two spellings of one edge - `L a + b +` here and `L b - a -` there - count as one link rather than two differences.

A label covers a whole class of equal sequences, not one segment. 
Variation graphs carry `A` as a node dozens of times over, and which copy in one file is which copy in the other is not something the two files agree on; pairing them off one by one would be a coin toss that every link touching them then inherits. 
Counting is still multiset-style - with three copies of a sequence in one file and two in the other, two match and the third is reported as unmatched - but a link between repeated sequences asks whether the other graph has that edge *between those sequences*, not between those particular nodes.

Segment ids are still listed for the unmatched surplus, but it is *how many* there are that is meaningful, not which ones file order happened to leave over.

Paths are compared by what they spell, not by what they walk over. 
`P` lines are paired by name across the two files and each pair's sequence is compared - so two graphs that chop one reference into different nodes still agree on it. 
Lengths are checked before the bases, so a path that differs in length costs nothing to reject. 
A file with no `P` lines at all is fine; it simply has no paths to report on.

Ids, node numbering, line order and `SR` ranks are not compared: a graph and its own `sort` output are identical to this command.

Two caveats on paths. 
Pairing is on the `P` line's own name, so a reference carried whole in one file (`chr1`) but split in the other (`chr1:0-1000`, `chr1:1000-2000`, ...) shares no name at all and reports as unmatched rather than as one difference. 
And the sequence comparison is byte-for-byte, hence case-sensitive.

Note: `--verbose` additionally lists the ids of the segments and links that only one of the two graphs carries, in that graph's own numbering.

Example:
```sh
$ ./akhal compare gfa graph.gfa sorted.gfa
Segments A: 12
Segments B: 12
Segments shared: 12
Segments only in A: 0
Segments only in B: 0
Links A: 15
Links B: 15
Links shared: 15
Links only in A: 0
Links only in B: 0
Paths A: 1
Paths B: 1
Paths identical: 1
Paths differing: 0
Paths only in A: 0
Paths only in B: 0
Path chr22: identical (48 bp)
[INFO] the graphs are identical
```

##### `gaf`
Compares two sets of alignments over the same graph, and reports how many reads both files put along the same walk.

The question is whether two aligners - or two runs, or two versions of one aligner - agree on where the reads go, so an alignment is reduced to the pair (read name, path) and nothing else. 
Coordinates, scores, mapping quality and CIGARs are not part of it: an alignment that starts 3 bp further into the same walk is the same alignment here.

The walk is compared in a canonical spelling. 
`>1>2<3` and `>3>2<1` are one walk read from its two ends, and an aligner that hit the read's reverse complement writes the second where another writes the first; both spellings are built and the smaller one is kept, so the two count as one alignment rather than two differences. 
A path field naming a stable sequence instead of walking nodes (`chr1`, as an unplaced alignment carries) has no orientation to flip and is kept verbatim.

Both files are then sorted by (read name, first node id of the canonical walk, the whole walk) and walked once, side by side - so the comparison costs a sort per file and a single pass, not a lookup per alignment. 
Reads are handled a name at a time: within a name both files carry, alignments pair off one-to-one on their walk, so a read aligned three ways here and twice there reports two pairs and one leftover rather than simply "matching". 
That splits the shared names three ways - every alignment paired, some paired, none paired - and those three plus the one-sided names account for every read in either file.

Malformed lines are skipped with a warning, and an empty file is not an error; it simply pairs with nothing.

Note: `--verbose` additionally lists every alignment that only one file carries, marking the ones whose read the other file never aligned at all - a read aligned *elsewhere* is the more interesting of the two.

Example:
```sh
$ ./akhal compare gaf minigraph.gaf graphaligner.gaf
Alignments A: 6
Alignments B: 5
Alignments shared: 3
Alignments only in A: 3
Alignments only in B: 2
Reads A: 5
Reads B: 5
Reads shared: 4
Reads only in A: 1
Reads only in B: 1
Reads on the same path(s): 2
Reads partly on the same path(s): 1
Reads on no shared path: 1
[INFO] the alignments differ
```

#### 5. `compact`
Folds every run of non-branching segments into a single node, the way `vg mod -u` and `odgi unchop` do.
Builders leave graphs full of chains - node after node joined by one link, with nothing branching anywhere along the way - and nothing is expressed by keeping them apart.

A join `u -> v` is taken when the link leaves `u`'s right end and enters `v`'s left end on the same strand (`L u + v +`, or the same edge read from the other side, `L v - u -`), it is the only link on either of those two ends, and it is a blunt join with no overlap - bases the two share would have to be dropped, and this never rewrites sequence it was not given.
A join that flips strand is left alone, and so is a self-loop.

Nothing else changes: the merged node keeps the first segment's id, links between runs are repointed at the nodes that swallowed their endpoints, and a path that walked a run step by step comes out walking the one node.
Two things stop a run where it would otherwise continue.
A path that starts or ends in the middle of one cuts it there, since a `P` line cannot stop halfway through a node.
And on an rGFA, two segments merge only if they sit on the same stable sequence at the same rank with contiguous offsets, so the merged node inherits a coordinate that still means something; a plain GFA has no tags to preserve and comes back plain.

Note on tags: the reader keeps a file's `SR` but not its `SN`, and takes each segment's offset from the layout of the path that owns it, so the `SN`/`SO` written here are the owning path and the offset along it rather than whatever the input said. 
Run `compact` before `gfa2rgfa`, not after, and both come out saying what they mean.

**Usage:**
```sh
./akhal compact <input .gfa file> [output .gfa file]
```

Note: If no output file is given, the compacted graph is written to standard output.

Example:
```sh
$ ./akhal compact graph.gfa
[INFO] 12 segment(s) in, 5 out: 7 folded into a neighbour
```

#### 6. `sort`
Topologically sorts a graph and renumbers its nodes `1..N` in the sorted order. Ordering uses Kahn's algorithm; ties in the ready set (the "hops") are broken alphabetically by node **sequence content**, so the result is independent of the input's node numbering — two graphs identical in topology and content sort the same way. Nodes in cycles, if any, are appended after the acyclic prefix. The sorted graph is re-emitted with the new ids (S/L/P lines remapped, link orientation and overlap preserved).

**Usage:**
```sh
./akhal sort <input .gfa file> [output .gfa file]
```

Note: If no output file is given, the sorted GFA is written to standard output.

#### 7. `rank`
Rewrites a graph's `SR:i:` ranks against its backbone. 
In rGFA, rank 0 is the reference and anything higher came from a sample; `rank` decides which segments are which and re-emits the graph.

By default the backbone is the graph's own `P` lines: every segment any path visits becomes rank 0, everything else rank 1. 
The path block is left exactly as it was read; only the `SR` tags change. 
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

#### 8. `vg2gfa`
Converts vg's native `.vg` format (gzip/BGZF-compressed Protobuf) to GFA. The `.vg` parser is hand-written pure C - no protobuf or libvgio needed, only zlib.

A `.vg` file has nowhere to store a path: a `Path` message only exists nested inside a `Graph` chunk, so each chunk carries the mappings of a path whose nodes it happens to hold, and one path arrives as thousands of pieces sharing a name. 
Each mapping carries a `rank` - its position along the whole path - and that is what the reader joins them by, so every path leaves as a single `P` line no matter how the file was chunked. 
This matters for anything that went through `vg convert`: those chunks follow the source graph's iteration order rather than the reference, so a path's pieces interleave instead of abutting and pasting them together in file order would be quietly wrong. 
A file whose mappings carry no usable ranks falls back to file order, with a warning.

**Usage:**
```sh
./akhal vg2gfa <input .vg file> [output .gfa file]
```

Note: If no output file is given, the GFA is written to standard output.

#### 9. `gfa2rgfa`
Labels a plain GFA as an rGFA, giving every segment the three tags that say where it sits on a real sequence: `SN:Z:` the name of that sequence, `SO:i:` the offset on it, `SR:i:` how far it is from the linear reference. 
A GFA carries none of them, but its `P` lines say everything needed to work them out.

One `P` line is declared the backbone - `--ref <name>`, else the graph's first - and its segments become rank 0, named after it, with offsets running the length of the walk. 
One `P` line is one path here too, so a reference the file splits over several of them contributes only the backbone line at rank 0; join it into a single `P` line first if that is not what you want.

Every other path is walked in turn. 
While it runs over ground that is already labelled it only follows along; where it leaves that ground, the segments it visits alone are one rank deeper, named after that path, and offset onward from the point it left - so a bubble carries the coordinate of the reference stretch it detours around, and the counting stops as soon as the path merges back down to a lower rank. 
A detour hanging off a rank-1 stretch is rank 2, and so on.

Where two paths disagree, nothing is invented. 
A segment that one path places at one offset and another at a different one has no single answer, so it keeps its rank and loses `SN`/`SO`, as does everything the walk reaches after it - until the walk arrives somewhere an authoritative path, the backbone above all, has already pinned down. 
Segments no path visits at all are left untagged. 
Both are counted in the summary.

This reads a **variation graph**, where the paths are haplotypes over a shared backbone. 
An assembly graph, whose paths are unrelated contigs, will label badly - there is no meaningful backbone for the rest to be an offset from.

**Usage:**
```sh
./akhal gfa2rgfa <input .gfa file> [output .rgfa file] [--ref <path name>]
```

Note: If no output file is given, the rGFA is written to standard output.

Example:
```sh
$ ./akhal gfa2rgfa graph.gfa --ref ref
S	1	AAAA	SN:Z:ref	SO:i:0	SR:i:0
S	2	C	SN:Z:ref	SO:i:4	SR:i:0
S	3	G	SN:Z:samp	SO:i:4	SR:i:1
S	4	TTTT	SN:Z:ref	SO:i:5	SR:i:0
S	5	GG	SN:Z:samp	SO:i:9	SR:i:1
```

Node 3 is the alternate allele of a SNP: rank 1, named after the sample path that explains it, and offset 4 - the reference coordinate of the node it detours around. Node 5 hangs off the end of the reference and carries on from where it stopped.

#### 10. `gfa2dot`
Writes a graph as a Graphviz digraph, for looking at rather than computing on.

One box per segment and one arrow per link. 
A segment carries its id and, when it is short enough to read, its bases; anything longer shows a length instead, and `--ids` drops the sequences entirely for a graph that is getting crowded. 
GFA is bidirected and DOT is not, so a link that flips strand says so on the arrow (`+/-`) and one that does not is left plain, which keeps a clean variation graph free of clutter. 
An overlap is labelled where there is one.

On a file that ranks itself the backbone is shaded apart from what hangs off it, so an rGFA's structure reads at a glance - pipe `gfa2rgfa` into this and the reference stands out from the bubbles.

This is meant for graphs small enough to look at: a bubble, a test case, a region pulled out of something bigger. 
Graphviz will accept a graph of any size and spend a very long time on it, so a graph over ten thousand nodes is reported before it is written.

**Usage:**
```sh
./akhal gfa2dot <input r/GFA file> [output .dot/.gv file] [--ids]
```

Note: If no output file is given, the DOT is written to standard output.

Example:
```sh
$ ./akhal gfa2dot graph.gfa graph.dot
[INFO] wrote 5 node(s) and 5 edge(s) to graph.dot; render with: dot -Tpng graph.dot -o graph.png
```

#### 11. `gaf2sam`
Converts a GAF file to a SAM file.

**Usage:**
```sh
./akhal gaf2sam <r/GFA file> <GAF file>  <FASTA file> <OUTPUT file>  [--simple]
```

Note: The reads should be stored in FASTA format and provided to the program. 
GAF file does not store sequences, hence, reads are needed when converting to SAM.

Note: `simple` option is optional. 
If provided, CIGAR string matches `(=)` and mismatches `(X)` will be replaced with sequence match `(M)`.

#### 12. `sampoke`
Validate SAM file (converted from gaf). It takes reference file and SAM to process CIGAR strings. 
Optionally, it can print the filtered SAM file that contains valid lines.

**Usage:**
```sh
./akhal sampoke <FASTA file> <GAF file> <OUTPUT file>
```

Note: Output file here is optional.

#### 13. `annotate`
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

#### 14. `annotget`
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
