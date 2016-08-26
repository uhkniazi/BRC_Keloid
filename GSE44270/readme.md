# Title
Keratinocyte and fibroblast gene expression in skin and keloid scar tissue

**DB ID** Project.id=1, Data.id=4  

[Evernote Link - Private](https://www.evernote.com/shard/s288/nl/38698211/8b9eb5c2-8222-41dd-badc-5b7ad1375eef?title=Keloid%20dataset%204)
## Comments
GSE7890. Author comments:  In order to improve understanding of the molecular basis of keloid scarring, we have assessed the genomic profiles of keloid fibroblasts and keratinocytes. 9 Keloid scar and 3 normal skin samples. 
The dataset has very small numbers of normal skin cell lines.  
**My Comments:** The Keratinocyte data results are similar to the original paper. The fibroblast results are also comparable to the paper with some interesting additions i.e. the role of glycosylation and cell motility. Other pathways like GPCR signalling, extra cellular matrix organisation ECM, are similar to the previous study results as well.  

## Scripts
1. download_data.R
  * Downloads the data in raw format from geo database.  
2. de_analysis.R  
  * Does a DE analysis and clustering for the data.  
3. de_analysis_simulated.R
  * Due to small sample sizes, repeated 2 group comparison instead of 3 and using simulated posterior means to calculate p-values.  
  
