# Pfmi3Dsc
Pfmi3Dsc: Protein functional mutation Identification by 3D Structure Comparison of protein Families
Selective pressures that trigger cancer formation and progression, shape the mutational landscape of somatic mutations in cancer. Given the limits within which cells are regulated, a growing tumor has access to only a finite number of pathways that it can alter. As a result, tumors arising from different cells of origin often harbor identical genetic alterations. Recent expansive sequencing efforts have identified recurrent hotspot mutated residues in individual genes. Here, we introduce Pfmi3Dsc, a novel statistical method developed based on the hypothesis that, functional mutations in a recurrently aberrant gene family can guide the identification of mutated residues in the family’s individual genes, with potential functional relevance. Pfmi3Dsc combines 3D structural alignment of related proteins with recurrence data for their mutated residues, to calculate the probability of randomness of the proposed mutation. The application of this approach to the RAS and RHO protein families returned known mutational hotspots as well as previously unrecognized mutated residues with potentially altering effect on protein stability and function. These mutations were located in, or in proximity to, active domains and domains related to protein-protein interactions.


the mutation porfile of all proteins must be downloaded from :
https://hive.biochemistry.gwu.edu/prd//biomuta/content/biomuta-master.csv
in order to run the program.

run main.py file to start processing the results.

The software outputs an HTML file with aligned residues and probabilities after receiving a UniProt Protein name as input.
