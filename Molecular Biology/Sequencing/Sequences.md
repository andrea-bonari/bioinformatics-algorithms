>[!note]
>Information is usually encoded in living cells in the form of sequences, such as DNA, RNA and protein sequences.
>
>Biological function largely depends on the three dimensional structure of large molecules, but any characteristics or functions can also be inferred from the one dimensional structure (the sequences).

>[!tip] Alphabets and strings
>Let $\Sigma$ be a finite set of characters, called the alphabet. A string over $\Sigma$ is any finite sequence of symbols from $\Sigma$.
>
>The length of a string $S$ is the number of symbols in $S$, and is denoted as $|S|$.
>
>The empty string is the unique string over $\Sigma$ of length $0$, and is denoted as $\varepsilon$.
>
>The set of all possible strings over $\Sigma$ of length $n$ si denoted as $\Sigma^{n}$. We can also define the Kleene closure of $\Sigma$ as:
>$$\Sigma^{*}=\bigcup_{n\in\mathbb{N}\cup\set{0}}\Sigma^{n}$$
>For any two strings $s$ and $t$ in $\Sigma^{*}$, their concatenation is defined as the sequence of symbols in $s$ followed by the sequence of characters in $t$, denoted as $st$.
>
>A string $s$ is a substring of $t$ if there exist strings $u$ and $v$ such that $t=usv$. The number of substrings on a string with length $n$ is:
>$$1+ \frac{n(n+1)}{2}$$
>The reverse of a string is a string with the same symbols but in reverse order. We can also define a string that's the reverse of itself as a palindrome.
>
>A string $s=uv\quad u,v\in\Sigma^{*}$ is said to be a rotation of $t$ if $t=vu$.


Protein sequences are usually represented by string over a 20-letter alphabet:
$$\Sigma=\set{A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y}$$

### Sequences comparison
>[!note]
>An important goal is to define suitable methods for the comparison of sequences, as such, we need to define suitable measures of distance or similarity, as well as efficient algorithms to compute them.

>[!tip] Sequencing error
>Sometimes the sequencing machine identifies a wrong base/nucleotide in a sequence. 
>
>Living cells can replicate their DNA and replication errors are naturally part of this process, such as base substitutions where the original pair is erroneously replaced by another one, or insertions/deletions of small DNA fragments (indels).
>
>Such errors cannot be ignored by algorithms employed in the data analysis.

>[!tip] Hamming distance
>The Hamming distance between two strings of equal length is the number of positions at which the corresponding symbols are different. However this criteria isn't good as this doesn't account for rotation, and is not a good measure for similarity in general.

>[!tip] Edit distance
>The edit distance is a way of quantifying how dissimilar two strings are to one another by counting the minimum number of operations required to transform one string into the other. We mostly use the Levenshtein definition of edit distance, where the allowed operations are removal, insertion and substitution. A possible mathematical definition is:
>$$\text{lev}(a,b)=\begin{cases}
>|a|\qquad&\text{if }|b|=0 \\
>|b|\qquad&\text{if }|a|=0 \\
>\text{lev}(\text{tail}(a),\text{tail(b)})\qquad&\text{if }a[0]=b[0] \\
>1+\min\begin{cases}\text{lev}(\text{tail}(a),b)\\ \text{lev}(a,\text{tail}(b))\\\text{lev}(\text{tail}(a),\text{tail}(b))\end{cases}\qquad&\text{otherwise}
>\end{cases}$$
>Where $\text{tail}(x)$ is the $x$ string without its first character, and $x[n]$ is the $n$-th character of a string $x$.
>
>From this mathematical definition, it's evident that computing the edit distance between two strings implies solving a minimization problem.

If two sequences are similar, we can often infer the function of one's biological function to the other.

### Sequence reading
>[!note]
>The DNA fragments inserted between the adapters are usually longer than the maximum read sequence length. We can either sequence only one end of the fragment (single-end sequencing) or both ends (paired-end sequencing).
>![[Pasted image 20260506173931.png]]

Advantages of paired-end sequencing include:
- Knowing an approximate distance between the two reads
- Helping to identify structural variants

