# B. subtilis Protein Secretion Toolkit

The efficient secretion of enzymes out of bacterial cells can greatly reduce the cost and complexity of enzyme production and purification.
Bacillus subtilis is widely used in industrial protein manufacturing (e.g., enzymes) due to its Generally Regarded As Safe (GRAS) status and its ability
to secrete enzymes tagged with the appropriate N-terminal peptide signal (SecTag).  We have designed a collection of 148 Bacillus subtilis SecTags that
can be used to screen for efficient secretion of any given protein.

### Genetic Design and Construction

Each plasmid contains an array of many SecTag variants, each preceded by a B. subtilis ribosome binding site.
RBS-SecTag pairs are bounded by BsaI restriction enzyme cut sites, and overhangs define each tag as an ‘RBS/Localization tag’
part in the extended Modular Cloning ‘MoClo’ assembly standard used for FreeGenes’s E. coli Protein Expression Toolkit.
Parts have overhang syntax B and T1 (5'--TACT---part---CCAT--3') from the Protein Expression Toolkit. Between each variant are 60 randomly generated base pairs.
In the middle of each random DNA spacer is a 8-base, blunt-cutting PmeI restriction site, which can be cleaved during assembly reactions to reduce the
chance of erroneous incorporation of two or more tandem RBS-SecTags into any genetic construct built from these plasmids.

### How to Use the Library

The B. subtilis SecTag library plasmids are designed for use in Golden Gate assembly reactions that attach both the RBS/SecTag pairs and a
fluorescent protein (FP) tag to an protein of interest, assembling a library of FP-tagged enzymes with different SecTags in a single Golden Gate reaction.
The assembled library can then be transformed into B. subtilis, and secretion can be screened for either by imaging a fluorescent ‘halo’ of secreted,
diffusing protein around colonies on agar plates, or by high-throughput culturing of colonies in 96-deep-well plates, followed by centrifugation and measurement
of secreted fluorescence via a plate reader.  The collection enables low-cost screening for efficient SecTag/enzyme pairs through a 1-pot Golden Gate reaction
starting from only six or fewer SecTag variant plasmids.

---

### Versions and Platemaps

|Version|Name|Platemap|Distribution Period|
|---|---|---|---|
|1.0|Freegenes Current Distribution|[Plate 1](https://github.com/Reclone-org/Open-DNA-Collections/tree/main/Bacillus%20subtilis%20Protein%20Secretion%20Toolkit/Platemaps/BPST-v1_0.csv)||

---

### Plasmids

|Name|ID|Freegenes ID|
|---|---|---|
| BsubSecTagLP_1 | [ODC_0331](https://github.com/Reclone-org/Open-DNA-Collections/blob/main/Bacillus%20subtilis%20Protein%20Secretion%20Toolkit/Plasmids_Genbank/ODC_0331.gb) | BBF10K_000485 |
| BsubSecTagLP_2 | [ODC_0332](https://github.com/Reclone-org/Open-DNA-Collections/blob/main/Bacillus%20subtilis%20Protein%20Secretion%20Toolkit/Plasmids_Genbank/ODC_0332.gb) | BBF10K_000486 |
| BsubSecTagLP_3 | [ODC_0333](https://github.com/Reclone-org/Open-DNA-Collections/blob/main/Bacillus%20subtilis%20Protein%20Secretion%20Toolkit/Plasmids_Genbank/ODC_0333.gb) | BBF10K_000487 |
| BsubSecTagLP_4 | [ODC_0334](https://github.com/Reclone-org/Open-DNA-Collections/blob/main/Bacillus%20subtilis%20Protein%20Secretion%20Toolkit/Plasmids_Genbank/ODC_0334.gb) | BBF10K_000488 |
| BsubSecTagLP_5 | [ODC_0335](https://github.com/Reclone-org/Open-DNA-Collections/blob/main/Bacillus%20subtilis%20Protein%20Secretion%20Toolkit/Plasmids_Genbank/ODC_0335.gb) | BBF10K_000489 |
| BsubSecTagLP_6 | [ODC_0336](https://github.com/Reclone-org/Open-DNA-Collections/blob/main/Bacillus%20subtilis%20Protein%20Secretion%20Toolkit/Plasmids_Genbank/ODC_0336.gb) | BBF10K_000490 |
