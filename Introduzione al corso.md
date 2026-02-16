## Obiettivi dell'insegnamento
Bioinformatics represents an application domain of increasing interest and importance for computer science. This discipline arose from the increasing need to develop appropriate computational methods for different problems of molecular biology and genetics, including the analysis of biological sequences (DNA, RNA, proteins). More specifically, recent next-generation sequencing technologies allow to generate massive amounts of biological sequence data that pose complex problems and require specific algorithms and computational analyses for their efficient and effective processing.

The goal of the course is to introduce students to some of the most important algorithms (and their underlying concepts) used in computational biology/bioinformatics, especially regarding those used for aligning biological sequences and identifying their locations in the reference genome. Also the assembly of a genome from short sequences, the modeling of functional biological sequences and the identification of mutations in a sequenced genome will be discussed. Each topic is treated both theoretically and practically. By discussing also alternative approaches for several of the algorithms, including their advantages and diadvantages, the students can learn how effective and efficient algorithms can be designed and when approximate solutions are preferable.

## Risultati di apprendimento attesi
| Dublin Descriptors | Expected learning outcomes |
|-------------------|----------------------------|
| Knowledge and understanding (DD1) | Students will learn how to align two or more sequences to each other including global and local sequence alignment and multiple sequence alignment; align short sequences produced by modern sequencing technology to a reference genome through read mapping; assemble a genome from short sequences produced by modern sequencing technology using genome assembly; model functional biological sequences by means of Hidden Markov Models; identify differences in a sequenced genome with respect to a reference genome through variant or mutation calling |
| Applying knowledge and understanding (DD2) | Given specific biological sequences such as DNA sequences including next generation sequencing reads and reference genomes students will be able to perform global and local sequence alignments; apply Burrows Wheeler transformations on genome sequences and use the transforms to perform exact and inexact matching of short sequences; understand the complexity of sequence analysis algorithms and discard inefficient algorithms; understand how functional sequences can be modeled by statistical means using a Hidden Markov Model |
| Making judgements (DD3) | Given the discussion of multiple alternative approaches and underlying concepts students will learn to evaluate the appropriateness of an algorithm based on its complexity as well as the complexity of the task; appreciate the importance of algorithms that approximate an optimal solution |
| Communication (DD4) | Not applicable |
| Lifelong learning skills (DD5) | Given the discussion of multiple alternative approaches and underlying concepts students can learn to approach the same algorithmic problem in different ways and from different sometimes opposing starting points; learn which aspects to pay attention to when studying a specific approach such as time and space complexity and inexact versus exact solutions |


## Argomenti trattati
1. (Very brief) Introduction to the biological background
    - Cells, biological molecules (DNA, RNA, proteins) and cellular processes
    - DNA sequencing techniques

2. Sequence alignment
    - Defining and measuring distances and similarities between sequences
    - Pairwise alignment of sequences based on similarity measures
    - Dynamic programming for constructing global and local sequence alignments (Needleman-Wunsch and Smith-Waterman algorithms)
    - Longest common substring problem
    - Multiple sequence alignment (more than two sequences)
        - Approximate versus optimal solution
        - Center star algorithm
        - Tree-based, progressive multiple sequence alignment
        - Weighted dynamic programming based on alignment/sequence profiles
        - Local multiple sequence alignment and motif scoring
    - Computational complexity of sequence alignment algorithms

3. Read alignment of next-generation sequencing reads
    - Suffix trie, suffix tree and suffix array
    - Burrows-Wheeler transform
    - FM index for searching substrings
    - Exact matching and inexact matching of sequencing reads
    - Computational complexity of the different approaches

4. De-novo genome assembly from next-generation sequencing reads
    - Modeling sequence data and computing expected results
    - Graph-based genome assembly (de Bruijn graphs, Eulerian path problem and Hierholzer’s algorithm)

5. Notions of variant calling
    - Single nucleotide variants/mutations
    - Copy number variants (deletions and duplications)
    - Germline variants in hereditary disorders versus somatic variants in cancer

6. Introduction to Hidden Markov Models for sequence modeling
    - Markov chains and Hidden Markov Model (HMMs)
    - Probability of a sequence being generated by an HMM
    - Viterbi decoding and posterior decoding
    - HMM training: learning model parameters


## Prerequisiti
Basic knowledge of algorithms and algorithmic complexity

## Modalità di valutazione
The assessment will be based on a written exam to be taken in the exam sessions defined by the school and covering all aspects in the syllabus. The exam will assign up to 60 points (+ 6 possible bonus points). 30 cum laude will be assigned when the total score exceeds 64 points.