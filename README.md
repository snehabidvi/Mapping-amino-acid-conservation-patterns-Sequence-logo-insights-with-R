# Mapping amino acid conservation patterns : Sequene logo Insights with R

In the field of bioinformatics, a sequence logo is a visual representation illustrating the conservation of nucleotide sequences (within DNA/RNA strands) or amino acid sequences (within protein molecules). The construction of a sequence logo is derived from aligned sequences, effectively illustrating both the consensus sequence and the variability present among the sequences. Sequence logos are commonly employed to represent sequence features, such as protein-binding regions within DNA or functional domains within proteins.

The overall height of a logo position is contingent upon the level of conservation within the respective multiple sequence alignment position. More conserved alignment positions produce high logo positions. The height of each character within a logo position is directly correlated to the frequency with which the corresponding amino acid is observed in the alignment column.

This R script is to perform multiple sequence alignment and visualize the conservation patterns of amino acids in peptide sequences using a **sequence logo** plot. It starts with reading a FASTA file containing peptide sequences (eg. plant_def.fasta) and converting them into an amino acid string set using the AAStringSet() function from the Biostrings package. Then, the msa() function from the msa package is used to align the amino acid sequences. 

The ggplot2 and ggseqlogo libraries are used to create a SeqLogo plot that displays the conservation patterns across the aligned sequences. Additionally, the consensus sequence is calculated using the consensusString() function with a threshold of 0.6, which determines the most frequent amino acids at each position in the alignment. This consensus string is annotated beneath the sequence logo plot. 


## Installing and loading necessary libraries 

```{r}


install.packages("Biostrings", "ampir", "msa", "ggplot2", "ggseqlogo")

library(Biostrings)
library(ampir)
library(msa)
library(ggplot2)
library(ggseqlogo)
```


## Input data 

The function read_faa() reads a fasta file to import data in rstudio. The input sequences can be given as a character vector.

```{r}
df <- read_faa("plant_def.fasta")   # replace with your file name 

```

    Arguments -
    1. file :- file path to the FASTA format file containing the protein sequences
  

The sequences are then converted to an amino acid string set for further analysis by readAAStringSet() function.

```{r}
aa <- AAStringSet(df$seq_aa)

```

    Arguments -
    1. filepath	:- A character vector (of arbitrary length when reading, of length 1 when writing) containing the path(s) to the file(s) to read or write. Reading files in gzip format (which usually have the '.gz' extension) is supported

    2. format	:- Either "fasta" (the default) or "fastq"

    3. nrec	:- Single integer. The maximum of number of records to read in. Negative values are ignored

    4. skip	:- Single non-negative integer. The number of records of the data file(s) to skip before beginning to read in records

    5. seek.first.rec	:- TRUE or FALSE (the default). If TRUE, then the reading function starts by setting the file position indicator at the beginning of the first line in the file that looks like the beginning of a FASTA (if format is "fasta") or FASTQ (if format is "fastq") record. More precisely this is the first line in the file that starts with a '>' (for FASTA) or a '@' (for FASTQ). An error is raised if no such line is found.
    Normal parsing then starts from there, and everything happens like if the file actually started there. In particular it will be an error if this first record is not a valid FASTA or FASTQ record
    Using seek.first.rec=TRUE is useful for example to parse GFF3 files with embedded FASTA data

    6. use.names :- TRUE (the default) or FALSE. If TRUE, then the returned vector is named. For FASTA the names are taken from the record description lines. For FASTQ they are taken from the record sequence ids. Dropping the names with use.names=FALSE can help reduce memory footprint e.g. for a FASTQ file containing millions of reads

    7. with.qualities	:- TRUE or FALSE (the default). This argument is only supported when reading a FASTQ file. If TRUE, then the quality strings are also read and returned in the qualities metadata column of the returned DNAStringSet object. Note that by default the quality strings are ignored. This helps reduce memory footprint if the FASTQ file contains millions of reads
    Using readQualityScaledDNAStringSet() is the preferred way to load a set of DNA sequences and their qualities from a FASTQ file into Bioconductor. Its main advantage is that it will return a QualityScaledDNAStringSet object instead of a DNAStringSet object, which makes handling of the qualities more convenient and less error prone 

    8. seqtype :- A single string specifying the type of sequences contained in the FASTA file(s)
    Supported sequence types:
    "B" for anything i.e. any letter is a valid one-letter sequence code.
    "DNA" for DNA sequences i.e. only letters in DNA_ALPHABET (case     ignored) are valid one-letter sequence codes
    "RNA" for RNA sequences i.e. only letters in RNA_ALPHABET (case ignored) are valid one-letter sequence codes
    "AA" for Amino Acid sequences i.e. only letters in AA_ALPHABET (case ignored) are valid one-letter sequence codes
    Invalid one-letter sequence codes are ignored with a warning

    9. x :- For writeXStringSet, the object to write to file.
    For saveXStringSet, the object to serialize

    10. append :-	TRUE or FALSE. If TRUE output will be appended to file; otherwise, it will overwrite the contents of file

    11. compress :- Like for the save function in base R, must be TRUE or FALSE (the default), or a single string specifying whether writing to the file is to use compression. The only type of compression supported at the moment is "gzip"
    Passing TRUE is equivalent to passing "gzip"

    12. ... :- Further format-specific arguments.

    If format="fasta", the width argument can be used to specify the maximum number of letters per line of sequence. width must be a single integer

    If format="fastq", the qualities argument can be used to specify the quality strings. qualities must be a BStringSet object. If the argument is omitted, then the quality strings are taken from the qualities metadata column of x (i.e. from mcols(x)$qualities). If x has no qualities metadata column and the qualities argument is omitted, then the fake quality ';' is assigned to each letter in x and written to the FASTQ file. If x is a QualityScaledXStringSet and qualities is not defined, the qualities contained in x are used automatically

    13. objname	:- The name of the serialized object

    14. dirpath	:- The path to the directory where to save the serialized object

    15. save.dups	:- TRUE or FALSE. If TRUE then the Dups object describing how duplicated elements in x are related to each other is saved too. For advanced users only

    16. verbose	:- TRUE or FALSE
    
    
## Performing multiple sequence alignment 

The msa() function provides a unified interface to the three multiple sequence alignment algorithms in this package: ‘ClustalW’, ‘ClustalOmega’, and ‘MUSCLE’.

Multiple sequence alignment (MSA) is performed on the peptide sequences using the msa package. The aligned sequences are then converted into a character format, which is required for plotting the sequence logo.

```{r}
maln <- msa(aa)

maln_chr <- as.character(maln)

```

Arguments -

    1. inputSeqs :- input sequences; this argument can be a character vector, an object of class XStringSet (includes the classes AAStringSet, DNAStringSet, and RNAStringSet), or a single character string with a file name. In the latter case, the file name is required to have the suffix ‘.fa’ or ‘.fasta’, and the file must be in FASTA format

    2. method	:- specifies the multiple sequence alignment to be used; currently, "ClustalW", "ClustalOmega", and "Muscle" are supported

    3. cluster :-	parameter related to sequence clustering; its interpretation and default value depends on the method; see msaClustalW, msaClustalOmega, or msaMuscle for algorithm-specific information

    4. gapOpening	:- gap opening penalty; the defaults are specific to the algorithm (see msaClustalW, and msaMuscle). Note that the sign of this parameter is ignored. The sign is automatically adjusted such that the called algorithm penalizes gaps instead of rewarding them

    5. gapExtension	:- gap extension penalty; the defaults are specific to the algorithm (see msaClustalW, and msaMuscle). Note that the sign of this parameter is ignored. The sign is automatically adjusted such that the called algorithm penalizes gaps instead of rewarding them

    6. maxiters	:- maximum number of iterations; its interpretation and default value depends on the method; see msaClustalW, msaClustalOmega, or msaMuscle for algorithm-specific information

    7. substitutionMatrix	:- substitution matrix for scoring matches and mismatches; format and defaults depend on the algorithm; see msaClustalW, msaClustalOmega, or msaMuscle for algorithm-specific information

    8. type	:- type of the input sequences inputSeqs; possible values are "dna", "rna", or "protein". In the original ClustalW implementation, this parameter is also called -type; "auto" is also possible in the original ClustalW, but, in this package, "auto" is deactivated. The type argument is mandatory if inputSeqs is a character vector or the file name of a FASTA file. If inputSeqs is an object of class AAStringSet, DNAStringSet, or RNAStringSet, the type of sequences is determined by the class of inputSeqs and the type parameter is not necessary. If it is nevertheless specified and the type does not match the class of inputSeqs, the function stops with an error

    9. order :- how the sequences should be ordered in the output object; if "aligned" is chosen, the sequences are ordered in the way the multiple sequence alignment algorithm orders them. If "input" is chosen, the sequences in the output object are ordered in the same way as the input sequences. For MUSCLE, the choice "input" is not available for sequence data that is read directly from a FASTA file. Even if sequences are supplied directly via R, the sequences must have unique names, otherwise the input order cannot be recovered. If the sequences do not have names or if the names are not unique, the msaMuscle function assignes generic unique names "Seq1"-Seqn to the sequences and issues a warning

    10. verbose	:- if TRUE, the algorithm displays detailed information and progress messages

    11. help :-	if TRUE, information about algorithm-specific parameters is displayed. In this case, no multiple sequence alignment is performed and the function quits after displaying the additional help information

    12. ...	:- all other parameters are passed on to the multiple sequence algorithm, i.e. to one of the functions msaClustalW, msaClustalOmega, or msaMuscle. An overview of parameters that are available for the chosen method is shown when calling msa with help=TRUE
    
  
## Generating the sequence logo Plot

geom_logo() function plots sequence logo as a layer on ggplot.

```{r}
p <- ggplot() +
  geom_logo(maln_chr, seq_type = "aa", stack_width = 1) +
  theme(axis.text.x = element_text(angle = 90, size = 5, vjust = 0.5)) +
  labs(title = "Amino acid conservation patterns in some defensin peptides")

```

Arguments -

    1. data	:- Character vector of sequences or named list of sequences. All sequences must have same width

    2. method	:- Height method, can be one of "bits" or "probability" (default: "bits")

    3. seq_type	:- Sequence type, can be one of "auto", "aa", "dna", "rna" or "other" (default: "auto", sequence type is automatically guessed)

    4. namespace :- Character vector of single letters to be used for custom namespaces. Can be alphanumeric, including Greek characters

    5. font	:- Name of font

    6. stack_width :-	Width of letter stack between 0 and 1 (default: 0.95)

    7. rev_stack_order :-	If TRUE, order of letter stack is reversed (default: FALSE)

    8. col_scheme	:- Color scheme applied to the sequence logo. See list_col_schemes for available fonts. (default: "auto", color scheme is automatically picked based on seq_type). One can also pass custom color scheme objects created with the make_col_scheme function

    9. low_col, high_col :-	Colors for low and high ends of the gradient if a quantitative color scheme is used (default: "black" and "yellow").

    10. na_col :- Color for letters missing in color scheme (default: "grey20")

    11. plot :- If FALSE, plotting data is returned

    12. ...	:- Additional arguments passed to layer params
    
    
## Adding Consensus Sequence to the Plot
   
The consensusString() function creates the consensus sequence from the consensus matrix based upon specified criteria.

```{r}
consensus_seq <- consensusString(maln, ambiguityMap = "X", threshold = 0.6)

p <- p +
  #annotate("text", x = seq(1, nchar(consensus_seq)), y = -0.3, 
           #label = strsplit(consensus_seq, "")[[1]], size = 2, hjust = 0.5)

p

```

Arguments -

    1. x :- An XStringSet or XStringViews object for consensusString
    
    2. ambiguityMap	:- Either a single character to use when agreement is not reached or a named character vector where the names are the ambiguity characters and the values are the combinations of letters that comprise the ambiguity (e.g. link{IUPAC_CODE_MAP}). When ambiguityMap is a named character vector, occurrences of ambiguous letters in x are replaced with their base alphabet letters that have been equally weighted to sum to 1 
    
    3. threshold :- The minimum probability threshold for an agreement to be declared. When ambiguityMap is a single character, threshold is a single number in (0, 1]. When ambiguityMap is a named character vector (e.g. link{IUPAC_CODE_MAP}), threshold is a single number in (0, 1/sum(nchar(ambiguityMap) == 1)]
    
    4. shift :- An integer vector (recycled to the length of x) specifying how each sequence in x should be (horizontally) shifted with respect to the first column of the consensus matrix to be returned. By default (shift=0), each sequence in x has its first letter aligned with the first column of the matrix. A positive shift value means that the corresponding sequence must be shifted to the right, and a negative shift value that it must be shifted to the left. For example, a shift of 5 means that it must be shifted 5 positions to the right (i.e. the first letter in the sequence must be aligned with the 6th column of the matrix), and a shift of -3 means that it must be shifted 3 positions to the left (i.e. the 4th letter in the sequence must be aligned with the first column of the matrix)

    5. width :- The number of columns of the returned matrix for the consensusMatrix method for XStringSet objects. When width=NULL (the default), then this method returns a matrix that has just enough columns to have its last column aligned with the rightmost letter of all the sequences in x after those sequences have been shifted (see the shift argument above). This ensures that any wider consensus matrix would be a "padded with zeros" version of the matrix returned when width=NULL
    
    
## Saving the Plot as a PNG File 

ggsave() is a convenient function for saving a plot. It defaults to saving the last plot that you displayed, using the size of the current graphics device. It also guesses the type of graphics device from the extension.


```{r}
ggsave("seqlogo.png", plot = p, height = 2, width = 7, units = "in")

```


Arguments -

    1. filename	:- File name to create on disk

    2. plot	:- Plot to save, defaults to last plot displayed

    3. device	:- Device to use. Can either be a device function (e.g. png), or one of "eps", "ps", "tex" (pictex), "pdf", "jpeg", "tiff", "png", "bmp", "svg" or "wmf" (windows only). If NULL (default), the device is guessed based on the filename extension

    4. path	:- Path of the directory to save plot to: path and filename are combined to create the fully qualified file name. Defaults to the working directory

    5. scale :- Multiplicative scaling factor

    6. width, height :- Plot size in units expressed by the units argument. If not supplied, uses the size of the current graphics device

    7. units :-	One of the following units in which the width and height arguments are expressed: "in", "cm", "mm" or "px"

    8. dpi :- Plot resolution. Also accepts a string input: "retina" (320), "print" (300), or "screen" (72). Applies only to raster output types

    9. limitsize :- When TRUE (the default), ggsave() will not save images larger than 50x50 inches, to prevent the common error of specifying dimensions in pixels

    10. bg :- Background colour. If NULL, uses the plot.background fill value from the plot theme

    11. create.dir :-	 Whether to create new directories if a non-existing directory is specified in the filename or path (TRUE) or return an error (FALSE, default). If FALSE and run in an interactive session, a prompt will appear asking to create a new directory when necessary

    12. ... :- Other arguments passed on to the graphics device function, as specified by device


