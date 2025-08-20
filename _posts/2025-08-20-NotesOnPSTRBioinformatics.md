---
layout: post
title: Notes on Pstrigosa bioinformatics data analysis
date: '2025-08-20'
categories: Analysis
tags: [reciprocal transplant, coding, Pegasus, rRNA, de novo transcriptome assembly]
---

Dr. Brad Weiler's new publication on diel transcriptomics of *Pseudodiploria strigosa* was just published on bioRxiv [here](https://www.biorxiv.org/content/10.1101/2025.08.06.668741v1), and his pipeline is available now on GitHub [here](https://github.com/delCampoLab/diel_pstr/blob/main/pstr_diel_rnaseq_submit.Rmd). 

I want to go through his steps and follow his code since he constructed a *de novo*  transcriptome for *P. strigosa*.

The sequencing he performed was "Poly-A capture/enrichment methodology by Novogene (Beijing, China) on the
158 Illumina NovaSeq PE 150bp. The methodology follows a standard workflow of fragmentation, reverse
159 transcription, cDNA synthesis, end-repair and A-tailing."

3' RNA-Seq doesn't necessarily equate to Poly-A selection/capture/enrichment. 3' RNA-Seq can be done with no selection, meaning that it would just sequence the 3'ends of RNA molecules but isn't targeting molecules with polyA tails. Poly-A selection/capture/enrichment is performed to enrich mRNA molecules over other types or RNA (typically rRNA) that dominate total RNA content. 

The question is, did my reciprocal transplant reads get any sort of selection? Looking back through my notes, I don't see anything about Poly-A selection.

## His pipeline: 
1. **FastQC**
2. Adapter trimming using **CutAdapt** wrapper script TrimGalore
3. Ribosomal RNA partitioned from total RNA using **SortMeRNA** and default fasta reference database
4. Remaining RNA mapped to index of genomes (Ostreobium quekettii,
Symbiodinium sp. 57, Breviolum sp. 58, Cladocopium sp. 57, Durusdinium sp. 59 183 , Eimeria tenella, Emitis
184 houghton, Sarcocystis neurona, Toxoplasma gondii, Malassezia restricta, Nitzschia putrida, Gracilaria
gracilis, and Uronemita sp.) to isolate non-coral RNA reads using **BBsplit**
5. Unmapped reads concaenated for *de novo* transcriptome assembly using **Trinity**
6. Decontamination using Standard PlusPF reference database in **Kraken2**
7. Reads aligned to transcriptome using **Salmon**
8. **STAR** used to align *Breviolum* binned reads
9. Gene transcripts were first predicted into longest open reading frames (ORFs) using **transDecoder**
10. Annotations using **EggNOG-mapper**
11. Gene ontology (GO) annotations supplemented using GoPredSim ProtT5 protein language model 65 195 derived topGO terms using **FANTASIA**

### Using SortMeRNA for removing rRNA
- Looking online, there seems to be some thoughts that removing rRNA reads is not required when Poly-A selection is performed during sequencing, as preferentially selecting polyA tails automatically creates a bias towards mRNA (see [this article](https://www.biostars.org/p/419845/) and [this article]()).
- However, some suggest that removing rRNA would save computation time and lead to a cleaner assembly ([see here](https://www.researchgate.net/post/Should_I_remove_rRNAs_from_transcriptome_data_while_moving_on_to_Denovo_assembly)).
- This [tutorial/paper](https://trhvidsten.com/docs/Delhomme-EpiGeneSysProtocol2015.pdf) is a good resource - it says that if GC content is abnormal, it could be a sign of high rRNA content (GC content of over 50% is typical for rRNA) and then SortMeRNA should be used. They also said "if there is any doubt, this step should be performed."

Based on the multiQC report of my trimmed reads, it looks like the average GC content is not normally distributed and is higher for some samples. To be safe, I'll follow Brad's code and run **SortMeRNA** on it.

<img width="959" height="528" alt="Screenshot 2025-08-20 at 1 54 00 PM" src="https://github.com/user-attachments/assets/82bfdc4b-c67a-42c8-bafb-a946812b8e2f" />

