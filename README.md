# HClink - The HierCC linking tool

## About

HClink is a powerful tool designed to bridge the gap between [Pathogenwatch](https://pathogen.watch/)
and [Enterobase](https://enterobase.warwick.ac.uk/) hierarchical clustering (HierCC) systems for bacterial genomics. It
serves the following key functions:

1. **Profile Matching**: HClink takes a [Pathogenwatch](https://pathogen.watch/) cgMLST (core genome Multi-Locus
   Sequence Typing) profile and matches it to the nearest [Enterobase](https://enterobase.warwick.ac.uk/) profile.

2. **Distance Calculation**: The tool calculates the genetic distance between profiles using the method described in
   the [Zhou et al. (2021) publication](https://academic.oup.com/bioinformatics/article/37/20/3645/6212647). This
   approach ensures consistency with established methodologies in the field.

3. **HierCC Code Inference**: Based on the calculated distance, HClink infers the HierCC (Hierarchical Clustering of
   cgMLST) code up to the allowed threshold. This provides valuable information about the hierarchical relationship of
   the input genome to known clusters.

4. **Flexible Input**: HClink can process both allele ST (Sequence Type) codes (numeric) and SHA1 checksum codes. The
   Pathogenwatch cgMLST code will use a SHA-1 checksum when there isn't a corresponding locus code in the database at
   the time it was built and run. hclink will map these to the corresponding allele if it has been added at Enterobase
   since. It can even be run entirely with checksum codes and these will be mapped to the current allele codes.

5. **Integration**: By linking Pathogenwatch and Enterobase clustering systems, HClink facilitates cross-platform
   comparisons and enhances the interoperability of different genomic databases.

This tool is particularly useful for researchers working on bacterial population genetics, epidemiology, and outbreak
investigations. It allows for rapid classification of new isolates within the context of existing genomic databases,
aiding in the understanding of bacterial evolution and spread.

For more information on the underlying concepts:

- [cgMLST and whole genome MLST](https://pubmed.ncbi.nlm.nih.gov/21143983/) - The BIGSdb paper which provided
  foundations for core genome MLST.
- [HierCC explanation](https://enterobase.readthedocs.io/en/latest/features/hierarchical-clustering.html)
- [Pathogenwatch cgMLST implementation](https://cgps.gitbook.io/pathogenwatch/technical-descriptions/cgmlst)

The software is designed with efficiency and ease of use in mind, making it a valuable addition to the toolkit of
microbial genomics researchers.

Hclink currently supports the _E. coli_ and _S. enterica_ schemes.

## Running HClink

### Via Docker

```
> cat my_genome.cgmlst.json | docker run --rm -i hiercc > assignment.json
```

### On the command line with uv

The easiest way to run hclink directly is using `uv`. This will automatically install all dependencies including the
correct version of python.

Use `--help` to get the complete list of options. In `assign` mode it can either read a file directly or take input
from STDIN by passing "-" as the filename.

```
# Help
> uv run hclink --help
```

```
# Build the latest E. coli scheme
> uv run hclink build my_ver ${API_KEY} ecoli
```

```
# Run against a PW cgMLST result file
> uv run hclink assign my_genome.cgmlst.json > assignment.json
> cat my_genome.cgmlst.json | uv run hclink - > assignment.json
```

## Building the pipeline

This software is intended to be distributed within Docker, so the recommended build path is to run the Dockerfile.
It will download all required files and dependencies, and build the latest database from fresh.

Replace `${VERSION}` with the version for reporting in the results.
Replace `${SPECIES}` with either:

* ecoli
* senterica

```
> docker build --build-arg VERSION=${VERSION} --build-arg API_KEY=${API_KEY} --rm hclink:latest .
```

## Input

A plain text file containing a JSON object with a field called `code` - extra fields will be ignored. This field should
be the allele codes joined by underscores. Note that the allele codes must be in the same order as given by Enterobase.

```
{
  "code": "1_1_50_44_2"
}
```

## Example output

```
{
  "st": "333640",
  "distance": 0,
  "hierCC_distance": 0,
  "gaps_both": 0,
  "gaps_a": 0,
  "gaps_b": 0,
  "hierCC": [
    [
      "d0",
      "333572"
    ],
    [
      "d2",
      "47842"
    ],
    [
      "d5",
      "47842"
    ],
    [
      "d10",
      "47842"
    ],
    [
      "d20",
      "47842"
    ],
    [
      "d50",
      "938"
    ],
    [
      "d100",
      "938"
    ],
    [
      "d150",
      "883"
    ],
    [
      "d200",
      "401"
    ],
    [
      "d400",
      "401"
    ],
    [
      "d900",
      "401"
    ],
    [
      "d2000",
      "44"
    ],
    [
      "d2600",
      "2"
    ],
    [
      "d2850",
      "2"
    ]
  ]
}
```

## About the checksum codes

Allele checksums are calculated using a SHA-1 digest of the _lower case_ DNA sequence of the allele. The code is
guaranteed unique for an individual allele sequence, though the schemes are not guaranteed to not have duplicate loci.

