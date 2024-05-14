---
name: "[Collection ID] Issue Report"
about: 'Create a report to help fix issues within the Open DNA Collections. '
title: "[ODC] Short description of issue "
labels: ''
assignees: ''

---

Thank you for taking the time to file an issue report. To help with addressing this issue, please provide as much of the following information as possible. 

### **Naming the Issue**
To help with pinpointing which of the Open DNA Collections that the issue is about, please use these codes in the "Issue Title" as follows, before adding a short description of the issue: 
- [OEC] = Open Enzyme Collection (parts in vector pOpen_v3). 
- [ORC] = Open Reporters Collection (parts in vector pOpen_v3). 
- [MDT] = Molecular Diagnostics Toolkit (ready-to-express constructs in vector pOpen_v3). 
- [RiDC] = Research in Diagnostics Collection (as MDT, but in vector pTI).
- [ePET] = E. coli Protein Expression Toolkit (parts in vector pOpen_v3).  
- [bPST] = Bacillus subtilis Protein Secretion Toolkit (parts in vector pOpen_v3). 
- [OYC] = Open Yeast Collection (parts in vector pOpen_v3). 
- [PURE] = PUREiodic Table Construction Kit (parts in pET28a). 
- [OP] = Open Plasmids. 
- [ODC] = Open DNA Collections - if the issue is related to all the above collections in general. 
e.g. [OP] pOpen backbones not available 

### **Describe the Problem**
A clear and concise description of what the issue is. 
Please include the [Collection ID], the associated plate wells, and the part names. 

### **Screenshots/Images**
If applicable, add screenshots to help explain the problem.

### **Sequence Data**
If applicable, add any sequence (DNA/protein) information associated.
e.g. 
>avGFP: 
>MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFS<sup>**H**</sup>GVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYN<sup>**F**</sup>NSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK
- Expected Y where there is a ***H***, which makes this BFP instead!
- In avGFP, the ***F*** should also be a Y. 
- Looks to be that a BFP sequence has been associated with this part when it should be avGFP. 

### **Expected State/Behaviour**
A clear and concise description of what you expected a part to have/step to happen.

### **Steps To Reproduce**
If applicable, add steps to reproduce the behaviour. Please indicate at which step the problem happened, and how. 
1. When ordering the collections(s) - e.g. no confirmation that my request was received. 
2. Upon receipt of the collection(s) - e.g. the foil covering the plates was loose. 
3. When recovering the collection(s) - e.g. no DNA could be recovered from Plate [ORC], wells A03 and H11; few/no colonies grew after transformation into E.coli DH5a (bought from NEB#C2987H with 10^9 cfu competency), nor with homemade batch made using NEB stock (with 10^6 cfu competency). 
4. Upon analysing sequencing results - e.g. the sequencing data does not match the sequence associated with the part (see above sections on "Sequence Data" and "Expected State/Behaviour"). 
5. When expressing/purifying/using protein - e.g. able to express all RiDC constructs and can visualise proteins on 10% SDS-PAGE (see "Screenshots/Images"), but low activity from the HIV-RT construct despite following the [expression and purification protocol](https://www.protocols.io/view/recombinant-expression-and-purification-of-hiv-1-r-ck8vuzw6). 
6. Documentation - e.g. could not find any information regarding the parts from https://github.com/Reclone-org/Open-DNA-Collections/ or https://stanford.freegenes.org/ . Nor were there any information of best practices for protein expression on https://www.protocols.io/workspaces/recloneorg-the-reagent-collaboration-network/
7. Generating issues - e.g. posted issue to discuss on https://forum.reclone.org/c/open-dna-collections/26 but had no feedback as of yet. 

### **Additional Context**
Add any other context about the problem here.

### **Follow-Up Actions (Assignees, Labels)**
Each issue will also be given a label to indicate what it refers to. For now, we are mainly using labels: sequence, part function, plasmid recovery, part proposal, enhancement, documentation. 

We welcome contributors to assign themselves to help with addressing these issues. Please let us know if you wish to do so. We're a small team, so any issues raised and support you can provide is greatly appreciated! 

Thanks for taking the time to fill out this Issue Report and for contributing to improving the Open DNA Collections!
