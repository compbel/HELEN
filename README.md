# HELEN (Heralding Emerging Lineages in Epistatic Network)

A computational framework for identifying epistatically linked sets of SAV alleles and merging them into haplotypes using statistical inference, population genetics, and graph theory.

## Algorithm:
![alt text](/img/algFlow.png)

HELEN (Heralding Emerging Lineages in Epistatic Networks): a computational framework for inference of viral variants as dense communities in epistatic networks


## Pre-requisites:
   - Matlab
   - Gurobi

## Instructions:

The main script is 
``[vocInfer, vocSupport] = HELEN_run(fastafile,fastaref,k_min,k_max,nSol,timeLimit,outdirRes,outToken)``

input:  
- ``fastafile`` -      fasta file with aligned viral sequences
- ``fastaref`` -       fastafile with the reference
- ``k_min`` -          minimal size of a densest candidate subgraph to be generated by HELEN
- ``k_max`` -          maximal size of a densest candidate subgraph to be generated by HELEN 
- ``nSol`` -           the number of densest subgraphs of a given size k (k_min <= k <= k_max) that are not contained in densest subgraphs of higher sizes
                         generated by HELEN
- ``timeLimit`` -      time limit (in seconds) for each ILP solver execution by HELEN
- ``outdirRes`` -      output directory, where intermediate HELEN results are saved. Set outdirRes = [], if saving of intermediate results is not required
- ``outToken`` -       an identifier to be attached to each saved output file (can be set if outdirRes ~= [])

 output:  
 - ``vocInfer`` -      a cell array of inferred genomic variants. Each variant is represented by an array of genomic positions defining this variant  
 - ``vocSupport`` -    vector of support values for inferred variants

 Example: [vocInfer, vocSupport] = HELEN_run('myData.fas','ref.fas',7,27,100,20000,'HELEN_results','myData')

For more information, see "Mohebbi, Zelikovsky, Mangul, Chowell,
Skums, Community structure and temporal dynamics of SARS-CoV-2 epistatic network allows for early detection of emerging variants with altered phenotypes"

## Additional notes
In addition, the repository contains the following data and scripts performing specific subroutines:

1) ``E = constructEpisNetwork(M,rho)`` : a script that construct epistatic networks from genomic data

   input:  
   - ``M`` - 	mutation matrix
   - ``rho`` -   p-value for detection of epistatically linked pairs
   
   output: 
   - ``E`` - 	list of edges of the constructed epistatic network
		
2) ``[varInfer, vocSupport] = HELEN_infer(G,k_min,k_max,nSol,timeLimit,outdirRes,outToken)``: a script that infers viral variants from the epistatic network G

   input:  
   - ``G`` -	epistatic network (as a graph object)
   - ``k_min``, ``k_max``, ``nSol``, ``timeLimit``, ``outdirRes``, ``outToken`` - see above
   
   output: see above

3) ``run_analysis`` : a script that generates the data and analysis results used in the paper "Mohebbi, Zelikovsky, Mangul, Chowell,
		Skums, Community structure and temporal dynamics of SARS-CoV-2 epistatic network allows for early detection of emerging variants with altered phenotypes"

4) ``collectResDrawPlotsSamp``, ``collectResDrawPlotsDensest``, ``collectResDrawPlotsInfer``: scripts that generate plots with the analysis results for sampling-based p-values, densest subnetworks, and haplotype inference used in the paper

5) HELEN data : secondary data generated for the paper

## Data:
Genomic data and associated metadata analyzed in this study were obtained from GISAID[^1].
SARS-CoV-2 epistatic networks and data derived from them can be downloaded from the following links, the links contain data from the complete dataset, the first truncated dataset, and the second truncated dataset respectively:

https://uconn-my.sharepoint.com/personal/pavel_skums_uconn_edu/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fpavel%5Fskums%5Fuconn%5Fedu%2FDocuments%2FMatlab%2FHELEN%20release%2FHELEN%5Fdata1&ga=1

https://uconn-my.sharepoint.com/personal/pavel_skums_uconn_edu/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fpavel%5Fskums%5Fuconn%5Fedu%2FDocuments%2FMatlab%2FHELEN%20release%2FHELEN%5Fdata2&ga=1

https://uconn-my.sharepoint.com/personal/pavel_skums_uconn_edu/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fpavel%5Fskums%5Fuconn%5Fedu%2FDocuments%2FMatlab%2FHELEN%20release%2FHELEN%5Fdata3&ga=1

## Citation:
Mohebbi F, Zelikovsky A, Mangul S, Chowell G, Skums P. Early detection of emerging viral variants through analysis of community structure of coordinated substitution networks. Nat Commun. 2024 Apr 2;15(1):2838. doi: 10.1038/s41467-024-47304-6. PMID: 38565543; PMCID: PMC10987511.

https://pubmed.ncbi.nlm.nih.gov/38565543/

## Acknowledgements:
We gratefully acknowledge all data contributors, i.e., the Authors and their Originating laboratories responsible for obtaining the specimens, and their Submitting laboratories for generating the genetic sequence and metadata and sharing via the GISAID Initiative, on which this research is based. The provided GISAID supplemental table includes a DOI where you can find all the associated authors and their originating laboratories.

[^1]: Khare, S., et al (2021) GISAID’s Role in Pandemic Response. China CDC Weekly, 3(49): 1049-1051. doi: 10.46234/ccdcw2021.255 PMCID: 8668406
