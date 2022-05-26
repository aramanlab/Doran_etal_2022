download Uniprot Reference dataset from https://www.uniprot.org/proteomes/?query=reference%3Ayes+taxonomy%3A%22Bacteria+%5B2%5D%22+excluded%3Ano

UP000538659 and UP000581646 were removed based on having 0 proteins measured

write the "Organism ID" column to text file `UP8881TaxID.txt`

download taxonkit

```bash
brew install brewsci/bio/taxonkit
```

Downloaded NCBI taxonomy data from NCBI taxonomy FTP

downloaded: taxdump.tar.gz version 2022-05-26 1:26:00 PM
into ~/.taxonkit
extract ``tar vxzf taxdump.tar.gz``


```bash
taxonkit lineage UP8881TaxID.txt | awk '$2!=""' > UP8881_Lineage.txt
13:20:14.973 [WARN] taxid 1264596 was merged into 2878388
13:20:14.981 [WARN] taxid 2478744 was merged into 1076744
13:20:14.986 [WARN] taxid 2029116 was merged into 2029117
13:20:15.040 [WARN] taxid 68210 was merged into 47760
```
these 4 warnings were found

reformat to machine readable format

```bash
cat UP8881_Lineage.txt \                                            
    | taxonkit reformat \
    | cut -f 1,3 \
    | tee UP8881_CleanLineage.txt
```
