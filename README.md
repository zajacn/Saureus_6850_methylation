# Saureus_6850_methylation

All scripts relevent to the project for course BIO253, funded by micro_innovation grant from UZH.

# Working Plan


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
* Methylation in gene bodies or in intergenic regions?
* Methylation in promoter regions?
* For methylations in the coding regions, what codon (position and AA) are affected?
* Is there a PA-SN effect (global analysis allSN vs allTSB)
    --> How does this look for each individual clone?
* Do we find methylations specific to a group of clones? (evolution effect)
* Comparing the ancestor but also the clones (PASN vs TSB) -> each clone individually
* 6mA (and 4mC) density per gene (no intergenics)




# Resources
https://www.nature.com/articles/s41467-024-54033-3;
https://www.nature.com/articles/srep29390
https://journals.asm.org/doi/full/10.1128/mbio.01773-24 #conserved 50S ribosomal protein methylation 
