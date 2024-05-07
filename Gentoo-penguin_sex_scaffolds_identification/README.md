
## 1. Samtools idxstats to obtain report of aligment
#!/bin/bash
source activate preSNPcalling

file="13A2
15A2
18A1
19A1
22A1
24A1
CRO05
CRO07
CRO08
G4
G5
G6
G7
GEN04
GEN05
GEN08
GEN09
GEN10
GSA10
GSA12
GSA14
GSA15
GSA7
GSA8
GSA9
GTP06
GTP07
GTP12
MQ9209
MQ9312
P0048
P0051
P0052
P0056
P0058
P0065
P0167
P0169
P0170
P0173
P0196
P0198
P0200
P0402
P0404
P0405
P0407
P0408
P0409
P0412
P1884
P1885
P1887
P1889
P1893
P1896
P2134
P2136
P2143
PAM10
PAM1
PAM3
PAM5
PAM7"

```for i in $file
do
samtools idxstats ${i}.bam > ${i}.idxstats
done
```


###### The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped read-segments and # unmapped read-segments. It is written to stdout. Note this may count reads multiple times if they are mapped more than once or in multiple fragments. (https://www.htslib.org/doc/samtools-idxstats.html)

