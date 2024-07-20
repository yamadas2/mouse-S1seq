# Directory organization
Map results will be stored in the following format. Project_XXX and NAME01, NAME02, ... are specified in a mapping script. The fastqc directory contains FASTQC outputs of fastq files. The GENOME directory contains map results against a reference genome <GENOME>. The trimming directory contains trimming reports and FASTQC outputs of trimmed fastq files.


- Format
```
MAP/
└── Project_XXX/
    ├── fastqc/
    ├── GENOME/
    │   ├── NAME01/
    │   └── NAME02/
    └── trimming/
        └── fastqc/
```

- Example
```
MAP/
└── Project_XXX/
    ├── fastqc/
    ├── mm10/
    │   ├── SAMPLE01/
    │   └── SAMPLE02/
    └── trimming/
        └── fastqc/
```
