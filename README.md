## Synopsis
```sh
# install gffio
git clone https://github.com/lh3/gffio
cd gffio && make

# examples with test files
gffio view test/hs38-gc42-part.gtf.gz > out.gff       # output GFF3
gffio view -t test/hs38-gc42-part.gtf.gz > out.gtf    # output GTF
gffio view -g test/hs38-gc42-part.gtf.gz > out.gff    # group by gene/mRNA
gffio view -L test/hs38-gc42-part.gtf.gz > out.gff    # select the longest mRNA
gffio gff2bed test/hs38-gc42-part.gtf.gz > out.bed    # output BED12
gffio gff2bed -i test/hs38-gc42-part.gtf.gz > out.bed # output introns in BED6

# examples without test files
gffio gff2fa hg38.fa test/hs38-gc42-part.gtf.gz > out.fa    # extract mRNA
gffio gff2fa -p hg38.fa test/hs38-gc42-part.gtf.gz > out.fa # extract proteins
gffio gff2fa -c hg38.fa test/hs38-gc42-part.gtf.gz > out.fa # extract CDS
```

## Introduction

gffio is a software tool to process GFF3 and GTF files. It can convert between
GFF3 and GTF, generate 12-column BED, extract CDS/transcript/protein sequences,
reorder features and select the longest CDS/transcript. Many gffio features are
also available in [gffread][gffread]. I implemented gffio mainly for a few use
cases I needed for my work.

gffio is largely a prototype. Let me know if see bugs or want more features.

## Limitations

* Only recognize one Parent field only.

* Requiring universally unique ID across all features.

* Not checking many possible errors in GFF.

[gffread]: https://github.com/gpertea/gffread
