# Saureus_6850_methylation

All scripts relevent to the project for course BIO253, funded by micro_innovation grant from UZH.

# Summary of the methods 

HiFi reads were converted to add kinetic information using ccs-kinetics-bystrandify script. The reads were then mapped to an NCBI reference genome using pbmm2, an official wrapper software for minimap2, to calculate interpulse duration ratios. Modification detection and motif prediction were performed using ipdSummary from the kineticsTools with default parameters. The workflow was performed as part of the SMRTLink v13.1.

# Working Plan and things to do
* ~~[JG] use protein fasta w/ StringDB to get functional classes (https://string-db.org/)~~ Results from StringDB not really useful since less structured
* ~~[JG] use protein fasta directly with eggNogg to get GO categories (http://eggnog6.embl.de/)~~ -> eggNog-mapper done, this one is very useful (2025-02-05: new rds file in resources)
* ~~[JG] calculate all possible contrasts for proteomics data (log2FC, T, FDR) for multiple contrasts~~ -> done, outside of this repos using prolfquapp
* ~~[JG] add proteomics data (contrast, log2FC, T, FDR) to mtx data frame (all combinations of contrasts?)~~ -> done, new rds file in resources
* [JG] join proteomics data w/ mtx in both directions (mtx to proteomics and proteomics to mtx) using locus_tag



# Brainstorming 16/01/2025
* For methylations in the coding regions, what codon (position and AA) are affected?
* Add proteomics data (contrast, log2FC, T, FDR) to mtx data frame (all combinations of contrasts?)
** What contrasts shall be calculated? (all clones vs ancestor? irrespective of growing condition?)
* Find very consistent methylations (found in all the 3 replicates) and indicate them in the mtx data frame
** ~~Here NZ reduced dataset to have it represented if and only if found in 2 out of 3 replicates from the same groupingConditions~~ -> done
* Regions (from all strains) in the genome.. how "uniform" is the genome methylated? (white or black spots?)
* PCA from this heatmap?
* try to bring in categories for the genes based on eggNOGG or StringDB?
* we have evidence that we can filter for IPD ratio > 2.5 for m6A (there is a bimodal distribution) ..

## Questions to potentially adress
* Intergene distances? (from gff -> NZ looks into it)
* methylation density: -> done by NZ: Do we see Methylation in gene bodies or in intergenic regions? (should also be in context of how much is coding from the full genome vs not coding)
* Methylation in promoter regions? (how to define promoter regions?) -> according to google: -35bp before start codon and -10bp after TSS (transcription start site)
* For methylations in the coding regions, what codon (position and AA) are affected?
* Is there a PA-SN effect (global analysis allSN vs allTSB)
    --> How does this look for each individual clone?
* Do we find methylations specific to a group of clones? (evolution effect)
* Comparing the ancestor but also the clones (PASN vs TSB) -> each clone individually
* 6mA (and 4mC) density per gene (no intergenics) # no difference with respect to conditions

# Brainstorming 11/04/2025
* LS: Think about how we would structure a potential manuscript
    * Idea is that in our previous proteomics paper we saw: Proteom mainly driven by environment rather then genotype --> How is this regulated? Our hypothesis = Methylation
    * I propose a Methylation centric approach
        * Focus on 6850_anc: 1. General Methylation pattern (TSB), Intergene distances, methylation density, Methylation in promoter regions 2. Compare 6850 (TSB) vs 6850 (SN): Where do we find differences 3. If we find differences, can we match these with proteomic data
* JG:

# Resources
https://www.nature.com/articles/s41467-024-54033-3; (methylation in SA from Group X)
https://www.nature.com/articles/srep29390 (methylation in SA from Group X)
https://journals.asm.org/doi/full/10.1128/mbio.01773-24 #conserved 50S ribosomal protein methylation 
https://www.nature.com/articles/s41467-023-38291-1 #conserved pattern of 6mA across 28 strains of L.paracasei
https://pmc.ncbi.nlm.nih.gov/articles/PMC11544029/ #genomic methylation promotes SA persistance, they found 6mA but needed to adjust their seq for 4mC
