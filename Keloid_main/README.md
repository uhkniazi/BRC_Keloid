# Title
Identifying what goes wrong during wound healing to cause excessive scars known as keloids.

**DB ID** Project.id=1, Data.id=2  

[Evernote Link - Private](https://www.evernote.com/shard/s288/nl/38698211/0a8bc5dc-e07d-40e4-9077-725c94e97bcd?title=00%20Keloid%20RNA-Seq)
## Comments
300 fastq files processed for quality checks. All seems fine, some trimming on the 3' side will be required on the first library.  

## Scripts
1. rna_seq_qa.R
  * Process a batch of FASTQ files and produces html output reports. Adds results to Project.File table.  
2. rna_seq_qa_2.R 
  * Similar to previous script but creates one large object with all ~ 260 FASTQ files from the first experiment, and creates one pdf document with read quality plots. **NOTE** memory and time intensive.
  
