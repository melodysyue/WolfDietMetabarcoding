Blocking Primer Design
================
Yue Shi, PhD candidate, University of Washington,
05/07/2018

-   [Background](#background)
-   [DPO (dual priming oligonucleotide) primers](#dpo-dual-priming-oligonucleotide-primers)
    -   [Principle](#principle)
    -   [Design strategy](#design-strategy)
    -   [Examples of DPO blocking primers](#examples-of-dpo-blocking-primers)
-   [Wolf blocking primer design](#wolf-blocking-primer-design)
    -   [Workflow](#workflow)
    -   [Result](#result)
-   [References](#references)

Background
----------

DNA metabarcoding has been widely used with environmental DNA (eDNA), DNA extracted from environmental samples including water, soil, scat etc. The diversity of sequences within a sample that can be detected by universal primers is often compromised by high concentrations of host DNA templates. If the DNA within the sample contains a small number of sequences in relatively high concentrations, then less concentrated sequences are often not amplified because the PCR favours the dominant DNA types. This is a particular problem in molecular diet studies, where predator DNA is often present in great excess of food/prey-derived DNA. A simple method is to use predator specific blocking primer (Vestheim and Jarman 2008).

DPO (dual priming oligonucleotide) primers
------------------------------------------

### Principle

Generally, there are two approaches to design blocking primers (Vestheim and Jarman 2008). One is to design a blocking primer which overlaps the 3' end of one universal primer, but extend into the predator-specific sequences and is modified with a C3 spacer at the 3'end. C3 spacer (3 hydrocarbons) CPG is a standard primer modification avaialble from most suppliers of custom oligonucleotides. Adding this C3 modification prevents elongation during PCR without noticeably influencing its annealing properties. In most cases, finding a predator-specific binding site next to a binding site of a universal primer is often difficult and very long primers often do not work (longer than 25 bp). Another approach is to design a dual priming oligonucleotide (DPO) containing two seprate priming regions joined by 5 consecutive deoxyinosin bases, called a poly(I) linker. The linker itself is not involved in priming. The 5' segment is the longest (18-25 bp) and crucial for positioning and stable annealing of the primer. The 3' segment is short (6-12 bp) and will only bind if there is already stable annealing of the 5' end. The short length and low annealing temperature of the 3' segment result in a low tolerance for mismatches and are reported to ensure target-specific binding. DPOs does not suffer from the limitations of a high Tm since the linker assumes a bubble-like structure resulting in two primer segments with distinct annealing properties. Furthermore, the bubble-like structure of linker efficiently prevents primer-dimer and hairpin structure formation.

### Design strategy

The position of the 3'end segment is determined first, at a site where 6–12 bases had a 40–80% GC content and the TM was not considered. Five deoxyinosines are designated for the poly(I) linker since they had generated the best result when 3–8 deoxyinosines were tested to determine the optimum length of the linker. The 5'end segment of the DPO is automatically determined by the sequence of the bases upstream of the 3'end segment and extended 18–25 bases, until the TM was &gt; 65 degC. In order to generate predator-specific DPOs, the length of the 5' segment of the DPO can be increased so that it had an even higher Tm (70-90 degC). The secondary structure and dimers are not considered in the design of the DPO because the 3' segment alone, which is physically separated by the linker from the 50-segment, is too short to form such structures stably (Chun et al. 2007).

### Examples of DPO blocking primers

**Short28SF-DPO-blkKrill** (Vestheim and Jarman 2008)

CCTGCCG (overlapping with universal forward primer, size=7bp)

CCTGCCGAAGCAACTAGCCCTGAAAATGGATGGCGCTCAAGCGTCCTC (GC% = 58%, size=48bp)

ACTCGACCGTTG (GC% = 58%, size=12bp)

**UrsusV5B2 (bear blocking)** (De Barba et al. 2013)

CCACTATGC (overlapping with universal forward primer, size=9bp)

CCACTATGCTTAGCCTTAAACAT (GC% = 39%, size=23bp)

AAATAATTTATTAAAC (GC% = 6%, size=16bp)

**HomoB (human blocking)** (De Barba et al. 2013)

CTATGC (overlapping with universal forward primer, size=6bp)

CTATGCTTAGCCCTAAACCTCAA (GC% = 43%, size=23bp)

CAGTTAAATCAACAAAACTGCT (GC% = 32%, size=22bp)

Wolf blocking primer design
---------------------------

### Workflow

-   Get 12S alignment from wolf and potential prey species (keep only one sequence per species). See `12S.fas`;
-   Extract the potentially amplication products with universal primers and align them. See `12S_amplified.fas`;
-   Extract the consensus sequence and compare it with that portion of wolf. See `wolf.vs.consensus.fas`;
-   Find out the wolf specific region and blast it to check. See `12S_wolf_specific.fas`;
-   Select 3' segment (6-12 bp, GC 40-80%), ideally 3'end should be wolf-specific;
-   Go upstream, skip 5 bps for the linker;
-   Go upstream and select 5' segment so that its Tm is greater than 65 degC.

### Result

**5' segment:** AGA TAC CCC ACT ATG CTT AGC CCT AAA CAT AGA TAA TTT TAC AAC (Size:45bp, GC: 36%, Tm: 65 degC)

**linker:** 5 deoxyinosines IIIII (AAAAT)

**3' segment:** AAT TCG CCA GAG G (Size:13bp, GC:54%)

References
----------

Chun, J Y, K J Kim, I T Hwang, Y J Kim, D H Lee, I K Lee, and J K Kim. 2007. “Dual priming oligonucleotide system for the multiplex detection of respiratory viruses and SNP genotyping of CYP2C19 gene.” *Nucleic Acids Research* 35 (6): e40–e40.

De Barba, M, C Miquel, F Boyer, C Mercier, D Rioux, E Coissac, and P Taberlet. 2013. “DNA metabarcoding multiplexing and validation of data accuracy for diet assessment: application to omnivorous diet.” *Molecular Ecology Resources* 14 (2): 306–23.

Vestheim, Hege, and Simon N Jarman. 2008. “Blocking primers to enhance PCR amplification of rare sequences in mixed samples a case study on prey DNA in Antarctic krill stomachs.” *Frontiers in Zoology* 5 (1): 12–11.
