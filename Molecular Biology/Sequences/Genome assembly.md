>[!note]
>Genome assembly consists in reconstructing the genome from short sequencing reads. It's inverse process consists into mapping the reads to a known reference genome, and it's called read mapping.

There are two major classes of assembly algorithms, OLC and DBG.

### Model of a genome assembly
>[!note]
>Let's consider an idealised genome, that represents a long random sequence of four bases and that doesn't contain repeats or other complex structures. Let's also consider a simple and error free sequencing strategy. We sample equal length fragments with starting points randomly distributed across the genome.
>
>The "shotgun" process can be compared to a process that samples bases from all genomes positions at random. The chance that any particular base is sampled is very low in a single sampling process. However we perform the sampling process a very large number of times.
>
>The Poisson distribution expresses the probability of a given number of events occurring in a fixed interval of time if these events are independent and identically distributed.
>$$f(k, \lambda)=P(X=x)= \frac{e^{-\lambda}\lambda^{k}}{k!}$$
>Where $k$ refers to the number of reads that overlap a certain genomic position, and $\lambda$ refers to the mean sequencing depth.
>
>Let's consider $G$ as the genome size, $L$ as the read length, $N$ as the number of reads and $n_{b}$ as the total number of sequenced bases. We know that:
>$$n_{b}=NL\qquad \lambda= \frac{n_{b}}{G}$$
>

We define a $k$-mer as a subsequence with $k$ nucleotides, in general there are $L-(k-1)$ $k$-mer subsequence in a sequence of length $L$ with $k<L$. In general the total number of $k$-mers win out whole genome sequencing (WGS) is: $$n_{k}=N\cdot(L-k+1)$$
While the coverage depth in term of $k$-mers is: $$d_{k}= \frac{n_{k}}{G}$$Therefore the ration between the coverage depth for bases ($\lambda$) and that for $k$-mers is then: $$\frac{\lambda}{d_{k}}= \frac{L}{L-k+1}$$
It's possible to estimate the genome size and base coverage as: $$G\approx \frac{n_{k}}{d_{k}}\qquad \lambda\approx \frac{L\cdot d_{k}}{L-k+1}$$
We can the estimate the probability that a given base will not be covered: $$P(X=0)= \frac{e^{-\lambda}\lambda^{0}}{0!}=e^{-\lambda}$$
And therefore, the probability of seeing at leas one read at a given position: $$P(X>0)=1-P(X=0)=1-e^{-\lambda}$$
>[!tip] Contigs
>Contigs are combinations of overlapping reads that represent contiguous sequence. The result of an assembly is a a set of contigs with gaps.
>
>![[Pasted image 20260604134100.png|center]]
>
>It's possible to piece together contigs, with paired-end sequencing we can assemble some contigs into "scaffolds".
>
>Let's consider an interval $l$ that is as long as a read ($L$ nucleotides), the probability that at least one read starts in $l$ is: $$P(X>0)=1-e^{-\lambda}$$
>Now lets consider a nucleotide at genomic position $i$, this nucleotide is in a gap between contigs if no read starts in the interval: $$[i-(L-1),i]$$
>This interval has length $L$ and therefore the probability that no read starts in it is $e^{-\lambda}$, we can then estimate the number of nucleotides in gaps across the entire assembly as: $$G\cdot e^{-\lambda}$$
>
>Each contig has a unique rightmost read $R$, and the probability that a given read is the rightmost read is the same as the probability that no other read starts within that read. Also, the number of contigs must be equal to the number of rightmost reads, thus the expected number of contigs is: $$N\cdot e^{-\lambda}$$
>The expected number of reads per contig is then: $$\frac{N}{Ne^{-\lambda}}= \frac{1}{e^{-\lambda}}$$
>We can also get the expected size of a contig: $$\frac{(1-e^{-\lambda})G}{Ne^{-\lambda}}$$
>Let $\theta$ refer to the minimum portion of $L$ that is required to detect an overlap, a group of reads combine to a contig if they are connected by overlaps of length greater or equal to $\theta L$. Given this we can estimate the expected number of contigs, given that we demand an overlap of at least $\theta L$ nucleotides: $$E[\# \text{contigs}]=Ne^{-(1-\theta) \frac{LN}{G}}$$
> 

>[!tip] N50
>The N50 is the measure used to estimate the quality of a genome assembly, to find it we arrange the contigs from largest to smallest, we find the position where the contigs cover $50\%$ of the total genome size, and the length of the contig in this position is defined as the N50.
>
>![[Pasted image 20260604135927.png]]
>
>The longer the N50 is, the better the assembly.

We need to find an algorithm that will allow us to take a collection of short NGS sequence reads, and to output a longer string representing the Genome that was sequenced.

### De Bruijn graph
>[!note]
>Imagine we have a function called $\text{composition}_{k}$ that takes a DNA sequence and returns a set of all $k$-mers contained in it, however since we don't know the original order of the $k$-mers in the genome we show the lexicographically. Let us now put each of the $k$-mers into the node of a graph and connect the graph by edges.
>![[Pasted image 20260604141324.png]]
>To find the sequence based only on a collection of $k$-mers subsequences we search for overlaps between them: $$\text{suffix}(k\text{-mer}_{i})=\text{prefix}(k\text{-mer}_{j})\qquad i\neq j$$
>In our graph the solution to our problem was a path that visited every node exactly once, that is an Hamiltonian path.
>
>We can instead label the nodes with the $k$-mer subsequences, and label the edges with these $k$-mers.
>![[Pasted image 20260604142224.png]]
>
>We can now merge identically labeled nodes in this graph, whilst retaining the edges.
>![[Pasted image 20260604142424.png|center]]
>This is the de Bruijn graph of the string, and to find the sequence we need a path that visits every edge exactly once, that is the Eulerian path problem, which fortunately has more efficient algorithms to solve compared to Hamiltonian paths.

>[!tip] Hamiltonian paths
>An Hamiltonian path is a path in an undirected or directed graph that visits each vertex exactly once (no need to use all edges).
>
>An Hamiltonian cycle is a Hamiltonian path that is a cycle.
>
>Determining whether Hamiltonian paths and cycles exist in a given graph is NP-complete, as $T(n)=\mathcal{O}\left(n^{k}\right)$.

>[!tip] Eulerian cycles
>An Eulerian cycle is a path that traverses every edge exactly once and returns at the end of the traversal to the start node.
>
>A traversable graph is one that can be drawn without taking a pen from the paper and without retracing the same edge. In such a case the graph is said to have an Eulerian path.
>
>Given that the degree of a vertex is the number of edges starting/ending at the vertex, when the degree of all the vertices is even, the graph is traversable (we have an Eulerian cycle), and we can draw it starting at any vertex (in directed graphs all vertices have $\text{degree}_\text{in}-\text{degree}_\text{out}=0$).
>
>If there are exactly two vertices of odd degree and all other vertices are of even degree, there is an Eulerian path (starting at one of the odd vertices), but no cycle (in directed graphs one vertex has $\text{degree}_\text{in}-\text{degree}_\text{out}=1$, one vertex has $\text{degree}_\text{in}-\text{degree}_\text{out}=-1$, while all others have $\text{degree}_\text{in}-\text{degree}_\text{out}=0$).
>
>If there are more than two odd vertices the graph cannot be traversed without repeating an edge.

The minimum read overlap $\theta L$ must be $k-1$ bases, such that the corresponding nodes of size $k-1$ can be merged.

To find an Eulerian path in a De Bruijn we use the Hierholzer's algorithm, with $T(E)=\mathcal{O}(E)$:
1. Make sure the Graph satisfies the degree requirements for an Eulerian path (or cycle) to exist
2. Initialise two stacks, a temporary `tpath` and a `epath` for the final solution
3. Choose a suitable starting vertex `v` (the one with one more outgoing edge for an Eulerian path, or any for a cycle)
4. Push `v` to `tpath`
5. Let `u=tpath.TOP`
6. Check all outgoing edges from `u`, if they've all been visited pop `u` from `tpath` and push it to `epath`, if not select a random outgoing edge $(u,x)$, push `x` to `tpath` and delete the edge $(u,x)$ from the set of available edges $E$
7. Repeat from step 5 until `tpath` is empty

The Eulerian path will be in `epath`.