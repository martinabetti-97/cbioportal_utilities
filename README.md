# cbioportal_utilities

## Data Mutations file 

```{r}
library(tidyr) 

files = list.files(path = ".",pattern = 'maf') 

test = read.csv2('file.maf',sep='\t',skip = 1) 

df = data.frame(matrix(ncol=ncol(test))) 

colnames(df) = colnames(test) 

for (file in files){ 

  df2 = read.csv2(file,sep='\t',skip = 1) 

  df = rbind(df,df2)} 

df = df %>% drop_na(Hugo_Symbol) %>% drop_na(Tumor_Sample_Barcode) 

df$n_ref_count <- NULL 

df$n_alt_count <- NULL 

df$n_depth <- NULL 

df <- df[,colSums(is.na(df))<nrow(df)] 

df = df[,c('Hugo_Symbol','Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','Tumor_Sample_Barcode','t_alt_count','t_ref_count','Variant_Classification','HGVSp_Short','Protein_position','SWISSPROT')] 

write.table(df,'data_mutations_gersom.maf',row.names = F, quote = F,sep= '\t') 

```
### Macchina AWS

#### Access 
```{bash}
ssh -i $key.pem ec2-user@15.160.44.209
```
Study folders in ```{bash}/home/cbioftp/study```

#### Generate dataset 
Compile dataset generation. This step both validate and update the dataset.


```{bash}
cd /home/ec2-user/cbioportal-docker-compose
docker-compose run cbioportal metaImport.py -u http://cbioportal:8080 -s /study/name --override_warning
```


### Genome-nexus 

If the annotation of the mutations file is not compliant, re-annotate with the following command:

```{bash}
docker run -v ${PWD}:/wd genomenexus/gn-annotation-pipeline:master --filename /wd/data_mutations_backup.txt  --output-filename /wd/output.txt 
```


## GITSIC  

To transform copy number variation values 

Gistic_format = round( log2(copy number + 2) - 1 )
