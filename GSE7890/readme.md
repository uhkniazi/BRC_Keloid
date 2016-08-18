# Title
Gene profiling of keloid fibroblasts shows altered expression in multiple fibrosis-associated pathways.  

**DB ID** Project.id=1, Data.id=1  

[Evernote Link - Private](https://www.evernote.com/shard/s288/nl/38698211/d7134d10-ba53-4756-8473-b6b76fe95496?title=Keloid%20dataset%201)
## Comments
GSE7890. Author comments: 500 genes seen as DE. Groupings used normal scars, keloid scars, male, female; samples grown with and without hydrocortisone. My Comments: When comparing the Normal cells vs Keloid after adjusting for treatment with or without hydrocortisone, the number of DE genes at 1% FDR is 852. These genes belog to cell cycle, hemostasis, signal transduction, extra-cellular matrix organisation pathways when doing an over represenation analysis with the reactome database.  

## Scripts
1. download_data.R
  * Downloads the data in raw format from geo database.  
2. de_analysis.R  
  * Does a DE analysis and clustering for the data.  

