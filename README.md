# About

HClink is a tool designed to bridge the gap between [Pathogenwatch](https://pathogen.watch/)
and [Enterobase](https://enterobase.warwick.ac.uk/) hierarchical clustering (HierCC) systems for bacterial genomics.
Genomes in Pathogenwatch are assigned the nearest possible HierCC code based on the cgMLST profile assigned using the
Pathogenwatch mlst tool.

While designed for integration into Pathogenwatch it also functions as a standalone tool. Once
the search database of Enterobase profiles is compiled, the tool is extremely fast, typically less than 2s per genome
queries against >100,000 profiles. To achieve this it uses the [usearch](https://unum-cloud.github.io/usearch/)
similarity search engine.

## Core functionality

1. **Profile Matching**: HClink takes a [Pathogenwatch](https://pathogen.watch/) cgMLST (core genome Multi-Locus
   Sequence Typing) profile and matches it to the nearest [Enterobase](https://enterobase.warwick.ac.uk/) profile.

2. **Distance Calculation**: The tool calculates the genetic distance between profiles using the method described in
   the [Zhou et al. (2021) publication](https://academic.oup.com/bioinformatics/article/37/20/3645/6212647) to ensure
   consistency with Enterobase.

3. **HierCC Code Inference**: Based on the calculated distance, HClink infers the HierCC (Hierarchical Clustering of
   cgMLST) code up to the threshold matched by the distance.

4. **Flexible Input**: HClink can process both allele ST (Sequence Type) codes (numeric) and SHA1 checksum codes. The
   Pathogenwatch cgMLST software will use a SHA-1 checksum when there isn't a corresponding locus code in the database
   at the time it was built and run. HClink will map these to the corresponding allele if it has been added at
   Enterobase since. It can even be run entirely with checksum codes and these will be mapped to the current allele
   codes.

5. **Integration**: This tools allows local searching of the Enterobase genomes with a query profile without installing
   a complex bioinformatics tool chain.

For more information on the underlying concepts:

- [cgMLST and whole genome MLST](https://pubmed.ncbi.nlm.nih.gov/21143983/) - The BIGSdb paper which provided
  foundations for core genome MLST.
- [HierCC explanation](https://enterobase.readthedocs.io/en/latest/features/hierarchical-clustering.html)
- [Pathogenwatch cgMLST implementation](https://cgps.gitbook.io/pathogenwatch/technical-descriptions/cgmlst)

Hclink currently supports the _E. coli_ and _S. enterica_ schemes.

## Author

This software was written by Corin Yeats for
the [Centre for Genomic Pathogen Surveillance](https://www.pathogensurveillance.net/) with help from Enterobase.

# Running HClink

## Important note

The HCLink database must be the same or more recent than the database used by the cgMLST tool, unless the cgMLST profile
is purely composed of checksums.

## Via Docker

```
> cat my_ecoli.cgmlst.json | docker run --rm -i hiercc:2025_03_06-2-ecoli > assignment.json
```

## On the command line with uv

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

## Building the software

### With Docker

This software is intended to be distributed within Docker, so the recommended build path is to run the Dockerfile.
It will download all required files and dependencies, and build the latest database from fresh.

Replace `${VERSION}` with the version for reporting in the results.
Replace `${SPECIES}` with either:

* ecoli
* senterica

```
> docker build --build-arg VERSION=${VERSION} --build-arg API_KEY=${API_KEY} --rm hclink:latest .
```

## With uv

```
uv build
uv pip install dist/*.whl
hclink --help
```

## Input

A plain text file containing a JSON object with a field called `code` - extra fields will be ignored. This field should
be the allele codes joined by underscores. Note that the allele codes must be in the same order as given by Enterobase.

```
{
  "code": "1_1_50_44_2_43..."
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
  "hierCC": [ [ "d0", "333572" ], [ "d2", "47842" ], [ "d5", "47842" ], [ "d10", "47842" ], [ "d20", "47842" ], [ "d50", "938" ], [ "d100", "938" ], [ "d150", "883" ], [ "d200", "401" ], [ "d400", "401" ], [ "d900", "401" ], [ "d2000", "44" ], [ "d2600", "2" ], [ "d2850", "2" ] ]
}
```

# About the checksum codes

Allele checksums are calculated using a SHA-1 digest of the _lower case_ DNA sequence of the allele. The code is
guaranteed unique for an individual allele sequence, though the schemes are not guaranteed to not have duplicate loci.

# Citing

If you use the outputs from this tool in a scientific publication, please cite Enterobase. The tool itself is currently
unpublished, but we recommend you upload your data to [Microreact](https://microreact.org) and genomes
to [Pathogenwatch](https://pathogen.watch) for sharing with the community (and cite those resources instead).

# Citations

1. Enterobase:
   Zhou, Z., Alikhan, N. F., Mohamed, K., Fan, Y., Achtman, M., & Brown, D. (2020). The EnteroBase user's guide, with
   case studies on Salmonella transmissions, Yersinia pestis phylogeny, and Escherichia core genomic diversity. Genome
   Research, 30(1), 138-152. https://doi.org/10.1101/gr.251678.119

2. Pathogenwatch:
   Argimón, S., Yeats, C. A., Goater, R. J., Abudahab, K., Taylor, B., Underwood, A., Sánchez-Busó, L., David, S.,
   Mariani-Kurkdjian, P., Lefevre, M., Baquero, F., Campos, J., Chaudhry, R., Chilton, D., Ciesielczuk, H., Comas, I.,
   Crook, D. W., de Mendoza, C., de Pinna, E., ... Aanensen, D. M. (2021). A global resource for genomic predictions of
   antimicrobial resistance and surveillance of Salmonella Typhi at pathogenwatch. Nature Communications, 12(1)
   ,2879. https://doi.org/10.1038/s41467-021-23091-2

3. cgMLST and whole genome MLST:
   Jolley, K. A., & Maiden, M. C. (2010). BIGSdb: Scalable analysis of bacterial genome variation at the population
   level. BMC Bioinformatics, 11, 595. https://doi.org/10.1186/1471-2105-11-595

4. HierCC method:
   Zhou, Z., Charlesworth, J., & Achtman, M. (2021). HierCC: A multi-level hierarchical clustering scheme for population
   genomic analyses. Bioinformatics, 37(20), 3645-3646. https://doi.org/10.1093/bioinformatics/btab234

# License

See the [License](/LICENSE) file