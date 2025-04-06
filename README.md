# AKHAL: Assembly Graph Analysis Tool

## Overview

`akhal` is a command-line tool designed to process and analyze r/GFA (Graphical Fragment Assembly) files. It provides functionality for validating, analyzing statistics, and converting GAF (Graphical Alignment Format) to SAM (Sequence Alignment Map).

## Installation 

You can build this tool by running:

```sh
make
```

## Usage

`akhal` provides three main commands:

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
Extract information from the r/GFA file.

**Usage:**
```sh
./akhal gaf2sam extract [OPTION] <r/GFA file> <OUTPUT file>
```

Options:
- **fa**: Reference genome. Output file should end with `.fa` or `.fasta`

#### 4. `gaf2sam`
Converts a GAF file to a SAM file.

**Usage:**
```sh
./akhal gaf2sam <FAI file> <r/GFA file> <GAF file> <OUTPUT file>
```

## License
is released under the BSD 3-Clause License, which allows for redistribution and use in source and binary forms, with or without modification, under certain conditions. For more detailed terms, please refer to the [license file](https://github.com/akmami/akhal/blob/main/LICENCE).

## Author
Developed by Akmuhammet Ashyralyyev.

**Note from the author**:

This tool is named after one of the most elegant horses, the [Akhal-Teke](https://en.wikipedia.org/wiki/Akhal-Teke). This breed is one of the oldest domesticated animals and is considered one of the most beautiful and intelligent horses in the world.
