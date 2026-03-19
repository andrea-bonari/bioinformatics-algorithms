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

### Sequence alignment
>[!note]
>A sequence alignment is an intuitive overview of differences between two strings. It's done by writing the two strings one below the other and highlighting where edit operation took place.
>$$\begin{align*}
>\textcolor{blue}{\text{I}}\textcolor{red}{\text{NT}}\text{E}\textcolor{red}{\text{N}}\textcolor{blue}{\text{- }}\text{TION}\\
>\textcolor{blue}{\text{-}}\textcolor{red}{\text{EX}}\text{E}\textcolor{red}{\text{C}}\textcolor{blue}{\text{U}}\text{TION}
>\end{align*}$$

### Alignment for minimum edit distance
>[!note]
>Let's say that $E(a,b)$ is the edit distance between two strings $a$ and $b$, and that $x(i)$ is the prefix on length $i\leq |x|$ of an $x$ string. Given that $E(a,\varepsilon)=|a|$ we can give the recursive definition:
>$$E(a(i),b(j))=\min\begin{cases}
>1+E(a(i),b(j-1)) \\
>1+E(a(i-1),b(j)) \\
>\begin{cases}
>E(a(i-1),b(j-1))\qquad&\text{if }a_{i}=b_{j} \\
>1+E(a(i-1),b(j-1))\qquad&\text{if } a_{i}\neq b_{j}
>\end{cases}
>\end{cases}$$
>If we set initially $i=|a|$ and $j=|b|$ we can compute $E(a,b)$.

Practically speaking, we can make a matrix of $|a|+1$ rows and $|b|+1$ columns, where each row and column is labelled with one character of $a$ and $b$ respectively. The first row and column are labelled with $-$ as an exception. We can use this table to compute, for each $i$ and $j$, the edit distance for the prefixes $a(i)$ and $b(j)$.

Once we've computed the edit distance, we need to get what's the optimal series of operations to transform one string into another.
To do this we start from the cell $(|a|,|b|)$ and perform backtracking, where we choose from the top, left and diagonal adjacent cells and select the minimum value (keep in mind that the left and top cells add $1$), once we select the cell we mark it and proceed iteratively until we reach the cell $(0,0)$.

![[Pasted image 20260319142601.png|center]]

>[!tip] Efficient backtracking
>For some cells there are multiple neighbouring cells which give the same minimum value, the choice between them is usually arbitrary.

The time complexity of the backtracking algorithm is $\mathcal{O}(|a|+|b|)$, while the overall complexity of the algorithm is $\mathcal{O}(|a|\cdot|b|)$.

>[!tip] Weighted edit distance
>By assigning a different weight to each operation ($w_{\text{id}}$ for insertions and deletions, and $w(a_{i},b_{j})$ for the symbols) we can define the weighted edit distance as:
>$$E(a(i),b(j))=\min\begin{cases}
>w_{\text{id}}+E(a(i),b(j-1)) \\
>w_{\text{id}}+E(a(i-1),b(j)) \\
>\begin{cases}
>E(a(i-1),b(j-1))\qquad&\text{if }a_{i}=b_{j} \\
>w(a_{i},b_{j})+E(a(i-1),b(j-1))\qquad&\text{if } a_{i}\neq b_{j}
>\end{cases}
>\end{cases}$$

### Alignment based on sequence similarity
>[!note]
>We can evaluate similarity instead of distance between two strings. To do that we need to define and computer measures of similarity, such as:
>- Percent identity: after the alignment we determine the percentage of the alignment which does not indicate an edit operation.
>- Longest common subsequence (LCS): derived by deleting some, all or no character, the order of the remaining elements remains unchanged.
>
>Let $a$ and $b$ be two strings, and let $a(i)$ and $b(j)$ be prefixes of length $0\leq i\leq |a|$ and $0\leq j\leq |b|$. Let $\text{LCS}(a,b)$ be the length of the longest common subsequence. We can compute it using:
>$$\text{LCS}(a(i),b(j))=\max\begin{cases}
>\text{LCS}(a(i),b(j-1)) \\
>\text{LCS}(a(i-1),b(j)) \\
>\begin{cases}
>1+\text{LCS}(a(i-1),b(j-1))\qquad&\text{if }a_{i}=b_{j} \\
>\text{LCS}(a(i-1),b(j-1))\qquad&\text{if }a_{i}\neq b_{j}
>\end{cases}
>\end{cases}$$

Practically we apply the same criteria as the edit distance table, but the backtracking algorithm adds one for diagonal cells and adds nothing for left and top cells.

### Global sequence alignment based on distance and similarity
>[!note]
>By combining the previous two approaches into a single one we define a similarity measure for two strings $a$ and $b$ such that:
>- differences have a negative effect: there's a negative weight for gap and mismatch, often called penalty.
>- conserved letters have a positive effect: there's a positive weight for match.
>  
>We can do that by using a combined alignment score $A(a(i),b(j))$ as follows:
>$$A(a(i),b(j))=\max\begin{cases}
>w_{d}+A(a(i),b(j-1)) \\
>w_{d}+A(a(i-1),b(j)) \\
>\begin{cases}
>w_{m}+A(a(i-1),b(j-1))\qquad&\text{if }a_{i}=b_{j} \\
>w_{s}+A(a(i-1),b(j-1))\qquad&\text{if }a_{i}\neq b_{j}
>\end{cases}
>\end{cases}$$
>Where $w_{d}<0$ is the negative gap penalty for insertions/deletions, $w_{s}<0$ is the negative mismatch penalty for substitutions and $w_{m}>0$ positive match score for the same symbol.

If we know the weight used when applying the algorithm, we can determine the final score of an alignment directly from the alignment itself by associating with each column its score and sum over all individual scores.

>[!tip] Gap penalties
>There's actually different models of gap penalties:
>- Linear gap penalty: gap penalty for a single symbol is multiplied by the total length of a gap.
>- Constant gap penalty: fixed negative score for the entire gap, regardless of its length.
>- Affine gap penalty: combines the constant and linear approaches by defining a constant gap opening penalty and a linear gap extension penalty.

