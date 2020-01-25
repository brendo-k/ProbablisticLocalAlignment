#Probablistic Local Alignment

Algorithm designed to align query sequence against a probablistic genome. The algorithm functions simlarly to BOWTIE using burrows wheeler transform to compress the genome. The compressed genome can then be used to find seeds using Forward Mapping indexing (FM-indexing). Seeds are then expanded using Smith-Waterman algorithm. More information about the algorithm can be found in the report.
