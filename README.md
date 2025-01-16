# Saureus_6850_methylation

All scripts relevent to the project for course BIO253, funded by micro_innovation grant from UZH.

# Working Plan and things to do
* [JG] use protein fasta w/ StringDB to get functional classes (https://string-db.org/)
* [JG] use protein fasta directly with eggNogg to get GO categories (http://eggnog6.embl.de/) -> eggNog-mapper is down today?
* [JG] calculate all possible contrasts for proteomics data (log2FC, T, FDR) for multiple contrasts
* [JG] add proteomics data (contrast, log2FC, T, FDR) to mtx data frame (all combinations of contrasts?)





# Brainstorming 16/01/2025
* For methylations in the coding regions, what codon (position and AA) are affected?
* Add proteomics data (contrast, log2FC, T, FDR) to mtx data frame (all combinations of contrasts?)
** What contrasts shall be calculated? (all clones vs ancestor? irrespective of growing condition?)
* Find very consistent methylations (found in all the 3 replicates) and indicate them in the mtx data frame
** Here NZ reduced dataset to have it represented if and only if found in 2 out of 3 replicates from the same groupingConditions
* Regions (from all strains) in the genome.. how "uniform" is the genome methylated? (white or black spots?)
* PCA from this heatmap?
* try to bring in categories for the genes based on eggNOGG or StringDB?
* we have evidence that we can filter for IPD ratio > 2.5 for m6A (there is a bimodal distribution)

## Questions to potentially adress
* Intergene distances?
* Do we see Methylation in gene bodies or in intergenic regions? (should also be in context of how much is coding from the full genome vs not coding)
* Methylation in promoter regions? (how to define promoter regions?) -> according to google: -35bp before start codon and -10bp after TSS (transcription start site)
* For methylations in the coding regions, what codon (position and AA) are affected?
* Is there a PA-SN effect (global analysis allSN vs allTSB)
    --> How does this look for each individual clone?
* Do we find methylations specific to a group of clones? (evolution effect)
* Comparing the ancestor but also the clones (PASN vs TSB) -> each clone individually -> Lukas
* 6mA (and 4mC) density per gene (no intergenics) -> Natalia



# Resources
https://www.nature.com/articles/s41467-024-54033-3;
https://www.nature.com/articles/srep29390
https://journals.asm.org/doi/full/10.1128/mbio.01773-24 #conserved 50S ribosomal protein methylation 
