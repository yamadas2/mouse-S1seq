# Directory organization
Fastq files need to be stored in the following format.
NAME will be the name of the output files. NAME can be any letters except space and underscore symbols. XXX, YYY and ZZZ can be any letters except space. Project_XXX and NAME will be specified in a mapping script.

- Format
```
FASTQ/
└── Project_XXX/
    └── YYY/
        └── Sample_NAME_ZZZ/
            ├── NAME_ZZZ_R1_001.fastq.gz
            └── NAME_ZZZ_R2_001.fastq.gz
```

- Example
```
FASTQ/
└── Project_example/
    └── data/
        ├── Sample_SAMPLE01_S01/
        │   ├── SAMPLE01_S01_R1_001.fastq.gz
        │   └── SAMPLE01_S01_R2_001.fastq.gz
        └── Sample_S1seq02_S01/
            ├── SAMPLE01_S01_R1_001.fastq.gz
            └── SAMPLE01_S01_R2_001.fastq.gz
```
