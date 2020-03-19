# gene_conversion_detect
## Usage
```
perl wtf_conversion_search.pl  -fsi flank_snp.identity.txt  -fa wtf.fasta -db dbname  -o output  -cv copy_number.var.txt
-h|--help:    print manual;
-fsi:   file recording pair-wise identies of flanking SNPs, which should involve 4 columns:strain1,strain2,identity(upstream),identity(downstream)
Format example:
JB938   JB22    0.99    1
Notice: when strainA vs strainB has occurred,strainB vs strainA is accepeted,but not recommended;
-fa:    file containing all wtf sequnces in all strains in one wtf locus;
-db:    database_name for blast to search donors;
-o:     output
-cv:    copy number varations
```
## schematic diagram
![image](https://github.com/xyhcelia/Readme_images/blob/master/gene_conversion_detect/wtf_gene_conversion_search_sketchmap.png)
