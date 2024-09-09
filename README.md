# HClink - The HierCC linking tool

## About
This program enables linking between Pathogenwatch and Enterobase clusters. It takes a Pathogenwatch cgMLST profile, and
matches it to the nearest Enterobase profile, calculating the distance as according to 
[their publication](https://academic.oup.com/bioinformatics/article/37/20/3645/6212647?login=true). It then infers the
HierCC code to the threshold allowed by the calculated distance.

## Running HClink
### Via Docker
```
> cat my_genome.cgmlst.json | docker run --rm -i hiercc > assignment.json
```

### On the command line
Use `--help` to get the complete list of options. In `assign` mode it can either read a file directly or take input 
from STDIN.
```
> python3 assign_enterobase_codes.py --help
> python3 assign_enterobase_codes.py assign my_genome.cgmlst.json > assignment.json
> cat my_genome.cgmlst.json | python3 assign_enterobase_codes.py - > assignment.json
```

## Building the pipeline
This software is intended to be distributed within Docker, so the recommended build path is to run the Dockerfile. 
It will download all required files and dependencies.

Replace `${VERSION}` with the version for reporting in the results.
Replace `${API_KEY}` with the Enterobase API key.

```
> docker build --build-arg VERSION=${VERSION} --build-arg API_KEY=${API_KEY} --rm hclink:latest .
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


