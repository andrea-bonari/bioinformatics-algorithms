>[!note]
>Variant calling is one of the key challenges in many areas of genomics research and diagnostics. Having aligned the fragments from DNA sequencing of one or more individuals to a reference genome, "single-nucleotide variant (SNV) calling" identifies variable sites, whereas "genotype calling" determines the genotype for each individual at each site. Also structural variants can be identified.

We assume we have data from WGS. Genomic variations include small mutations such as SNPs/SNVs, indels or structural variants such as duplication, deletions, inversions or translocation.

A Single Nucleotide Polymorphism (SNP) is a variant that is polymorphic within a population, while a Single Nucleotide Variant (SNV) is a variant identified in an individual genome.

Structural variants can either be balanced (no change in the amount of DNA) or unbalanced (amount of DNA increased or decreased).

>[!tip] PHRED quality scores
>Two measures are fundamental when evaluating how confident we are with the information obtained from high-throughput sequencing:
>- Mapping quality: quality/confidence with a given read maps to a unique genomic location
>- Base quality: quality/confidence with which an individual base was sequenced or "called"
>
>Both of these measures are usually reported as PHRED quality scores, defined as: $$Q_\text{PHRED}=-10\log_{10} p$$
>Where $p$ is the probability that the corresponding base call (or read mapping) is wrong.

>[!tip] Pileups
>The pileup format facilitates SNP/indel calling and bried alignment viewing by eyes, where each line consists of chromosome, $1$-based coordinate, reference base, the number of reads covering the site, read bases and base qualities.
>
>![[Pasted image 20260605130711.png|center]]
>
>Where read bases are indicated as:
>- A dot: sequence match to the reference base on the forward strand
>- A comma: sequence match on the reverse strand
>- "ACGTN": a mismatch on the forward strand
>- "acgtn": a mismatch on the reverse strand

For germline variants we can do the following:
- Variant calling: identifying variants and their genomic locations
- Genotype calling/genotyping: determining the "genotypes" of individual genomic positions

Possible genotypes in the germline (considering two copies) include:
- Homozygous reference: both copies have the base of the reference genome
- Heterozygous variant: one copy has an alternate base
- Homozygous variant/alternate: both copies have the same alternate base

### Germline variants
>[!note]
>We need to characterise each column of the sequence alignment, in other words, the alignment of all reads which map to the corresponding position in the reference genome. Lets assume, of $n$ reads which map to a given genomic position, we observe $k$ reference bases and $n-k$ alternate bases, we can make these assumptions:
>- If the true genotype is homozygous reference, then the $n-k$ observed alternate bases must be sequencing errors
>- If the true genotype is homozygous alternate, then the $k$ observed reference bases must be sequencing errors
>- If the true genotype is heterozygous, then we can approximate the probability of observing $k$ reference bases (out of $n$) as: $$\text{dbinom}(n,k,\underbrace{p}_{0.5})= \binom{n}{k} \frac{1}{2^{n}}$$

Early NGS studies basically filtered base calls according to quality and then used a frequency filter:
1. Typically, a base quality filter of PHRED score $Q\geq 20$ was used
2. Then, the following thresholds were used for the alternate base $b$ frequency

$$f(b)= \frac{n-k}{n}$$

| $f(b)$      | Genotype call        |
| ----------- | -------------------- |
| $[0,0.2)$   | Homozygous reference |
| $[0.2,0.8]$ | Heterozygous         |
| $(0.8,1]$   | Homozygous variant   |
The frequency $f(b)$ is often also referred as Variant Allele Frequency (VAF) or Mutant Allele Frequency (MAF).

### MAQ
>[!note]
>To address problems such as undercalling of heterozygous genotypes, loss of information and no measure of confidence in the previous approach, we have a number of probabilistic methods, we will discuss MAQ, which is based on the MAP estimation procedure.
>

We consider the probability that a read $z$ comes from position $u$ of a reference sequence/genome $\mathcal{R}$. Let $\set{\mathcal{M}\space\mathcal{M}}$ refer to the set of mismatched positions in the read, we have that: $$p(z|\mathcal{R}, u)=\prod\limits_{i\in\set{\mathcal{M}\space\mathcal{M}}}10^{- \frac{Q_{i}}{10}}= 10^{-{\frac{\sum\limits_{i}Q_{i}}{10}}}$$
We now calculate the posterior probability $p_{s}(u|\mathcal{R},z)$ that the read $z$ actually maps to position $u$, using Bayes' law: $$p_{s}(u|\mathcal{R},z)= \frac{p(z|\mathcal{R},u)\cdot p(u|\mathcal{R})}{\sum\limits_{v}p(z|\mathcal{R},v)\cdot p(v|\mathcal{R})}$$
MAQ actually uses various heuristics to calculate the mapping error probability, and thus the mapping quality of the read $z$, resulting in a PHRED scaled score $Q_{z}$ for the probability that the read is wrongly mapped. MAQ then redefines base qualities as the minimum of base quality $Q_{i}$ and mapping quality $Q_{z}$: $$Q_{i}=\min(Q_{i}, Q_{z})$$
MAQ then uses MAP methodology to identify the genotype that maximizes the probabilities: $$\begin{align*}
&p(\langle a,a\rangle|D)\qquad\text{homozygous reference}\\
&p(\langle a,b\rangle|D)\qquad\text{heterozygous}\\
&p(\langle b,b\rangle|D)\qquad\text{homozygous alternate}
\end{align*}$$
Where $a$ refers to the reference base, $b$ to the alternate base and $D$ to the data. MAQ assumes the priors for the genotypes are: $$p(\langle a,a\rangle)=p(\langle b,b\rangle)= \frac{1-r}{2}\qquad p(\langle a,b\rangle)= r$$
Where $r$ is the probability of observing a heterozygous genotype. MAQ typically uses $r=0.001$ for new SNPs, and $0.2$ for known SNPs. MAQ thus calls the genotypes as: $$\widehat{g}=\arg\max_{g}p(g|D)$$
The quality of this genotype call can then be calculated as: $$Q_{g}=-10\log_{10}(1-P(\widehat{g}|D))$$
We now need a way of calculating the probability of observing $k$ errors in $n$ nucleotides in the alignment $\alpha_{n,k}$, if we assume that errors arise independently, and error rates are identical for all bases, then we can use a binomial distribution: $$\text{dbinom}(n,k,p=\varepsilon)= \binom{n}{k}\varepsilon^{k}(1-\varepsilon)^{n-k}$$
In practice, MAQ doesn't use a binomial distribution, but a heuristic that reflects the probabilities of observing an alignment with the given pattern of per base error probabilities: $$\alpha_{n,k}=c'_{n,k}\prod\limits_{i=0}^{k-1}\varepsilon_{i+1}^{\theta^{i}}$$
Here, $\varepsilon_{i}$ is the $i$-th smallest base error probability for the $k$ observed error, $c'_{n,k}$ is a constant and the $\theta$ is a parameter that controls the frequency of errors.

### VarScan2
>[!note]
>VarScan2 is a heuristic, pileup-based approach towards somatic variant calling algorithm. For each position in the genome, the following steps are performed in parallel for the tumour sample and the matched normal sample:
>1. Check the minimum coverage requirement for both samples, by default $n\geq 3$ reads with base quality $Q\geq 20$ at that position
>2. Based upon the observed read bases, determine a genotype for each sample individually. By default, a variant allele must be supported by at least $n-k\geq 2$ independent reads and at least $8\%$ of reads
>3. Variants are called homozygous if supported by at least $75\%$ of all reads at a position

If the genotypes of the two samples do not match, then their read counts are evaluated by a one-tailed Fisher's exact test in a two-by-two contingency table.

|                      | ALT base           | REF base           | total            |
| -------------------- | ------------------ | ------------------ | ---------------- |
| Tumour sample reads  | $N_{T,\text{alt}}$ | $N_{T,\text{ref}}$ | $N_{T}$          |
| Control sample reads | $N_{C,\text{alt}}$ | $N_{C,\text{ref}}$ | $N_{C}$          |
| Total                | $N_{\text{alt}}$   | $N_{\text{ref}}$   | $N_{\text{tot}}$ |
And the probability to have exactly $N_{T,\text{alt}}$ reads in the tumour sample is: $$P(X=N_{T,\text{alt}})=\frac{\binom{N_{T}}{N_{T,\text{alt}}}\binom{N_{C}}{N_{C,\text{alt}}}}{\binom{N_\text{tot}}{N_\text{alt}}}$$
In practice, we actually need the probability to have at least $N_{T,\text{alt}}$ ALT reads in the tumour sample: $$P(X\geq N_{T,\text{alt}})= \sum\limits_{i=N_{T,\text{alt}}}^{N_{T}} \frac{\binom{N_{T}}{i}\binom{N_{C}}{N_{C,\text{alt}}}}{\binom{N_\text{tot}}{i+ N_\text{alt}}}$$
If $P(X\geq N_{T,\text{alt}})<0.05$ (significant $p$-value) then the variant is called somatic. Another possibility is loss of heterozygosity (LOH), where there's a heterozygous variant in the control, but it's homozygous in the tumour.

### Read depth analysis
>[!note] 
>Analysis of read depth can identify deletions/duplications. Read coverage is often plotted as the $\log_{2}$ of the number of reads which start within windows of fixed size. It's best done with WGS data.
>
>For tumour samples, we often compare the coverage in the tumour sample (center) with the coverage in the control sample (top) to evidence somatic differences (bottom).
>
>![[Pasted image 20260605142312.png|center]]
>
>We use the $z$-score to describe the read count measurements, which represent how many standard deviations $\sigma$ the read count $x$ of a given window is above or below the mean read count $\mu$: $$z= \frac{x-\mu}{\sigma}$$

First of all we convert the $z$-score to its upper and lower tail probability: $$P_{i}^{\text{Upper}}=p(Z> z_{i})\qquad P_{i}^{\text{Lower}}=p(Z< z_{i})$$
For an interval $\mathcal{A}$ of $\mathcal{l}$ consecutive windows, we call it an unusual event if (for duplications): $$\max_{i}\left\{P_{i}^{\text{Upper}}|i\in\mathcal{A}\right\}< \left( \frac{\mathcal{l}}{L}\cdot \text{FPR} \right)^{\frac{1}{\mathcal{l}}}$$
Where $L$ is the length of the entire chromosome and $\text{FPR}$ is the desired false-positive rate for the entire chromosome. The same procedure is now done separately for deletions, using the formula: $$\max_{i}\left\{P_{i}^{\text{Lower}}|i\in\mathcal{A}\right\}< \left( \frac{\mathcal{l}}{L}\cdot \text{FPR} \right)^{\frac{1}{\mathcal{l}}}$$