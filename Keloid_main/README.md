# Title
Identifying what goes wrong during wound healing to cause excessive scars known as keloids.

**DB ID** Project.id=1, Data.id=2  

[Evernote Link - Private](https://www.evernote.com/shard/s288/nl/38698211/0a8bc5dc-e07d-40e4-9077-725c94e97bcd?title=00%20Keloid%20RNA-Seq)
## Comments
300 fastq files processed for quality checks. All seems fine, some trimming on the 3' side will be required on the first library.  

## Scripts
1. download_data.R
  * Downloads the data in raw format from geo database.  
2. de_analysis.R  
  * Does a DE analysis and clustering for the data.  