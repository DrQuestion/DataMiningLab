# Bash commands used during project

We downloaded 537 unique cases(551 totali) from GDC database.
```
gunzip ./*/* -r to extract the data from the zipped files.
```
to see how many and what sample sites (we have = 34)
```
cut -f7 gdc_sample_sheet.2019-10-10.tsv | tail -n537 | cut -d '-' -f2 | sort | uniq | wc -l 
```
to check if alla the files in the folders have the same number of genes (=60483)
```
for i in $(ls); do wc -l ./$i/*.FPKM.txt ; done | cut -f1 -d ' ' | sort | uniq 
```
to count how many samples per patient
```
cut -f7 gdc_sample_sheet.2019-10-11.tsv | cut -f3 -d '-' | sort | head -n551 | uniq -c | sort -n 
```
Generate list of duplicates Cases ID (9 of them):
```
cut -f6 gdc_sample_sheet.2019-10-10.tsv | grep -v Case | sort | uniq  | sort -n | tail -n 9 > duplicates_cases_id.txt
```
