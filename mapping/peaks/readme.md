# Directory organization
Peak center files need to be stored in the following format. NAME can be any letters except space. 


- Format
```
peaks/
└── GENOME/
    └── primary/
        └── NAME.bed
```

- Example
```
peaks/
└── mm10/
    └── primary/
        └── hotspot.center.SPO11.B6.Cell.mmc2.bed
```


- Peak center file format
```
<chromosome name>   <center position - 1>   <center position>   <peak height>
```

- Peak center file example
```
chr1	100	101	10.5
chr1	105	106	100.5
chrY	135	136	50
```

The list of SPO11-oligo hotspot centers can be obtained from Lange et al., 2016, PMID: 27745971, doi: 10.1016/j.cell.2016.09.035.
