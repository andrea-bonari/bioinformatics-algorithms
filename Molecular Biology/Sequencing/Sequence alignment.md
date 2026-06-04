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

>[!tip] Substitution matrix
>It's possible to simplify the process using a substitution matrix. Given an alphabet $\Sigma$, the substitution matrix $\sigma\in\mathbb{N}^{|\Sigma|\times|\Sigma|}$, which gives for every pair of symbols $s_{i},s_{j}\in\Sigma$ the weight of the substituting $s_{i}$ by $s_{j}$ in an alignment. The matrix is symmetrical.
>
>For $i=j$ the weight represents the positive score assigned by a match.
>
>The usage of a substitution matrix changes the recursive computation to:
>$$E(a(i),b(j))=\max\begin{cases}
>w_{d}+E(a(i),b(j-1)) \\
>w_{d}+E(a(i-1),b(j)) \\
>\sigma(a_{i},b_{j})+E(a(i-1),b(j-1))
>\end{cases}$$

### Local sequence alignment
>[!note]
>Let $a=a_{1}a_{2}\cdots a_{n}$ and $b=b_{1}b_{2}\cdots b_{m}$ be two strings over the alphabet $\Sigma$. Let also $\sigma(s_{i},s_{j})$ be a substitution matrix with scores for all $s_{i},s_{j}\in\Sigma$ and $w_{d}$ be a gap penalty.
>
>To find a two substrings that produce the alignment with the maximum possible score we use the Smith-Waterman algorithm. This algorithm uses the following computation:
>$$M(i,j)=\max\begin{cases}
>0 \\
>w_{d}+M(i,j-1) \\
>w_{d}+M(i-1,j) \\
>\sigma(a_{i},b_{i})+M(i-1,j-1)
>\end{cases}$$
>