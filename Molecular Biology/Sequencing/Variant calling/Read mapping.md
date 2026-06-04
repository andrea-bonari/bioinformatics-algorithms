>[!note]
>Considering DNA sequencing data (WGS or ES), referenced-base assembly often follows the goal of finding the differences between an individual's genome and the reference genome for the corresponding species, rather than characterising the genome of that species in the first place.
>
>Reads being aligned/mapped to the reference genome are a prerequisite for variant calling, the computational identification of high quality mutations from sequencing data.
>

Search algorithms like Smith-Waterman are quite slow, but faster search algorithms based on preprocessing of the text to build a substring index exist. A substring index is a data structure which gives substring search in a text or text collection in sublinear time.

>[!tip] Naive string search
>Given a genome with $G$ characters and read of $L$ characters, the simplest string algorithm would simply slide the pattern across the genome, extending it letter for letter as long as there is a match, with $T(G,L)=\mathcal{O}(G\cdot L)$.

### Trie
>[!note]
>A trie (from retrieval), is a multi-way tree structure useful for storing strings over an alphabet. It's defined as the smallest tree over an alphabet $\Sigma$ such that each edge of the trie is labelled with one character $c\in\Sigma$, a node has at most one outgoing edge labelled for each $c\in\Sigma$ (at most $|\Sigma|$), and each key is spelled out along some path starting at the root.
>
>![[Pasted image 20260604154206.png|center]]
>
>$ is a symbol that does not appear anywhere in our genome template $T$, we define it to be lexycografically less than our other characters, therefore the $ enforces a lexicographic rule that we know from dictionaries (for instance, "over" comes before "overture"). The $ also ensures that no suffix will be considered as a prefix of any other suffix.
>
>If $L'$ is the maximum length of any read, then the runtime of the trie algorithm is $T(L',G)=\mathcal{O}(L'\cdot G)$ for matching and $T(n_{b})=\mathcal{O}(n_{b})$ for trie construction, where $n_{b}$ is the combined length of our reads. This is quite time efficient, however the amount of memory required for the trie is in the worst case proportional to the total length of the reads, which can be enormous $S(n_{b})=\mathcal{O}(n_{b})$.

>[!tip] Suffix Tries
>Instead of using different words for the trie, let's use suffixes of $, hence, each path from the root to a leaf represents a suffix, and each suffix is represented by a path from the root to a leaf.
>
>Each substring $T$ is represented by a path from the root, because every substring is a prefix of some suffix of $T$. Thus to search for a substring $S$, start at the root and follow the edges labelled with the characters of $S$. If at some point there is no outgoing edge for the next character of $S$, then $S$ is not a substring of $T$.
>
>![[Pasted image 20260604155537.png|center]]

To construct a suffixTrie we use the following algorithm:
```
SuffixTrie(T):
	T += $
	root = {}
	for i <- 1 to length(T)
		n = root    # n is the current node
		for c in T[i:]    # for each char in the i-th suffix
			n[c] = {}    # add outgoing edge to n if needed
		
		n = n[c]    # switch current node to child node
	
	return root
```

The `followPath` algorithm returns the node at the end of the path or NULL if there is no path:
```
followPath(T, S):
	root = SuffixTrie(T)
	n = root    # n is the current node
	for i <- 1 to length(S)
		c = S[i]    # i-th char of S
		
		if c not in n then
			return NULL    # not found
		
		n = n[c]    # switch current node to child node
	
	return n
```

The `hasSubstring` algorithm checks if `followPath` does not "fall off" the tree and returns NULL:
```
hasSubstring(T, S):
	return followPath(T, S) != NULL
```

>[!tip] Suffix Trees
>Given $m=|T|$, the worst case for space complexity is $S(m)=\mathcal{O}\left(m^{2}\right)$, that is too large. To optimise this we can combine non branching paths into a single edge with a string label:
>
>![[Pasted image 20260604160914.png|center]]
>
>By doing this we make sure that our tree has at most as many internal nodes as a full binary tree, because if an internal node has more than $2$ children we have fewer parent nodes, thus there are less or equal to $2m$ total nodes, so $S(m)=\mathcal{O}(m)$. This works but the total length of the edge labels is still $\mathcal{O}\left(m^{2}\right)$. To reduce size complexity further we can store in each edge only the offset and length of the original labels in the original string, so there are only two integers for one edge ($S(m)=\mathcal{O}(1)$).
>
>![[Pasted image 20260604162331.png|center]]
>
>This is a Suffix Tree, and we can also build it directly using Ukkonen's linear time on-line suffix tree construction algorithm.
>
>To fine all matches of a sequencing read $P$ in a genome $T$, letting $k$ be the number of matches and $n$ be the length of $P$, the search is then $T(n,k)=\mathcal{O}(n+k)$.
>
>Suffix trees, although having a linear algorithm has a really high constant factor $c$ ($\mathcal{O}(n)$ means $c\cdot n+d$), therefore it can be quire impractical.

### Suffix arrays
>[!note]
>The suffix array, at leas in its simplest incarnation, requires only $4\text{ bytes}$ per character of the input sequence.
>
>Given a string $T$, the suffixes of this strings are: $$T[0\cdots N-1]\quad T[1\cdots N-1]\quad\cdots\quad T[N-1\cdots N-1]$$
>A naive implementation of the suffix array basically manipulates an array of pointers to the suffixes of $T$.

A naive approach to build a suffix array is:
1. Form all possible suffixes from the input string $T$, and with each suffix associate it's position/index in the original string
2. Sort lexicographically to bring repeated strings together, we only keep the index array

>[!tip] Manber and Myers algorithm
>The naive approach to build the suffix array is not very efficient, the Manber and Myers algorithm has $T(n)=\mathcal{O}(n\log n)$. It works by:
>1. Sorting only the first character of the suffixes (using key-indexed counting sort)
>2. Recursive phase ($i$): given an array of suffixes sorted on the first $2^{i-1}$ characters, create an array of suffixes sorted on the first $2^{i}$ characters

### Burrows Wheeler Transform
>[!note]
>The Burrows Wheeler Transform (BWT) applies a reversible transformation to a block of input text. The transformation itself does not compress the data, but reorders it to make it easy to compress with simple algorithms.
>
>To make the BWT of a input text $T$:
>1. Form all rotations of the input text $T$, appended with the $
>2. Sort the rotated strings lexicographically, the result is the Burrows Wheeler matrix
>3. The BWT is simply the last column of the Burrows Wheeler matrix
>
>We denote the BWT of an input string $T$ as: $$\text{BWT}(T)$$

>[!tip] BWT from Suffix Array
>The Burrows Wheeler matrix is nearly the same as the suffixes referred to by the suffix array of the same string.
>
>![[Pasted image 20260604165803.png|center]]
>
>We can now write an algorithm to create $\text{BWT}(T)$ from the suffix array $\text{SA}(T)$ of $T$ by noting that position $i$ of the BWT corresponds to the character that, in the original string, is just to the left of the $i$-th suffix in the SA.
>
>We can now create the BWT as follows: $$\text{BWT}(T)[i]=\begin{cases}
>T[\text{SA}[i]-1]\qquad&\text{if }\text{SA}[i]>0\\\$ &\text{if }\text{SA}[i]=0\end{cases}$$
>The naive algorithm for this:
>```
>bwtFromSuffixArray(T):
>	sa = constructSuffixArray(T$)
>	L = length(sa)
>	bwt = new string[L]
>	
>	for i <- 0 to L - 1
>		if sa[i] = 0
>			bwt[i] = $
>		else
>			bwt[i] = T[sa[i] -1]
>	
>	return bwt
>```

>[!tip] LF Mapping property
>For any character, the $T$-ranking of characters in the first column $F$ of the BW matrix is the same as order of characters in the last column $L$.
>
>The $T$-ranking of the character at any given position is the number of times that an identical character has preceded it in $T$. 
>
>![[Pasted image 20260604170651.png|center]]
>
>The $B$-ranking of a character at a specific position is the number of times the character has occured in the column above the current position.

The reversibility of the BWT depends on the LF Mapping property. We can reverse the BWT by following this visual algorithm:

![[Pasted image 20260604171459.png|center]]

Note that we can do this process starting only from the $\text{BWT}(T)$, which corresponds to the $L$ column, if we count the number of each character in $\text{BWT}(T)$, we can easily reconstruct the chunks of characters in the $F$ column of the BWM, using the cumulative index property.
### FM Index
>[!note]
>The Full-text index in Minute space (FM index) uses the BWT and some other auxiliary data structures to generate a fast and efficient index for small patterns within a larger string $S$.
>
>The main data structures are the $F$ and $L$ column from the BWM. Note that the $F$ column itself is not stored because it can be represented as an array of integers.
>
>Given a string $P$ to be searched in our genome $T$, we look for all the rows that start with the last letter of $P$ in $F$, we then look in $L$ to identify those rows whose last letter corresponds to the second to last letter in $P$, we now use the LF mapping to find the rows in $F$ that begin with these characters.
>
>![[Pasted image 20260604172735.png|center]]
>
>If no letters can be found in $L$, the search pattern $P$ doesn't exist in $T$.

The naive implementation of the algorithm has many problems:
- We need to find the preceding characters efficiently, at the worst case this has $T(|T|)=\mathcal{O}(|T|)$
- We still need a way to get the $B$-ranks of the characters in $L$
- We still need a way of figuring out at what positions matches occur in $T$

To fix these issues we construct a tally table that precalculates the number of each specific character in $L$ up to every row.

![[Pasted image 20260604173347.png|center]]

After we found all rows beginning with the last character of $T$, we need to find rows with the second to last character in the $L$ column. Say the range of rows is $[i,j]$, we look in the tally table in row $i-1$ and in row $j$. We now know how many characters $r$ occur in $L$, in that range, therefore there are only that many lookups instead of $\mathcal{O}(|T|)$. 

This approach however needs to store $\mathcal{O}(|T|\cdot |\Sigma|)$ integers, therefore we store only at every $k$-th row.

![[Pasted image 20260604173830.png|center]]

This also provides us with the $B$-ranks, if we subtract $1$ from the tally.
Finally, to find the position of the match we can use a similar approach, taking advantage of the LF mapping.

![[Pasted image 20260604174204.png|center]]

>[!tip] BWT/FM Index algorithms for read mapping
>There are lots of published read aligners for genomic resequencing. Perhaps the best known amongst them use the BWT/FM Index plus lots of bells and whistles.

### Burrows Wheeler Aligner
>[!note]
>BWA is an algorithm based on BWT, that takes into account inexact matching.
>
>Let $W,T$ be two strings over the alphabet $\Sigma$, let $|W|=m$ and $|T|=n$, let $d(a,b)$ be the distance measure between two strings $a$ and $b$, let $d_\text{max}$ be a threshold for the maximum allowed distance, and finally let $T^{p}=T[p\cdots p+(m-1)]$ be the substring of $T$ having length $m$ starting at position $p$.
>

>[!tip] Prefix trie
>Let's consider the prefix trie for the string $T$, that is equivalent to a suffix trie but a path from a leaf to the root gives a unique prefix of $T$.
>
>Like for BWT, we start searching for substrings with the last character of the substring using DFS, and for inexact matching we bruteforce it to allow for up to $d_\text{max}$ mismatches.
>
>![[Pasted image 20260604183240.png|center]]

>[!tip] Suffix array interval
>All occurrences of a substring $W$ in the original string $T$ appear next to each other in the suffix array, because a substring is equivalent to the prefix of a suffix of $T$ and we have lexicographically sorted all suffixes, therefore a substring $W$ can be represented as a SA interval: $$\left[\underline{R}(W),\overline{R}(W)\right]$$
>![[Pasted image 20260604183502.png|center]]
>
>The BWA paper presents the method for calculating the SA interval of the query word $W=W[0,\cdots m-1]$ that we've already seen.
>

To inexact search over a suffix array we use the algorithms:
```
InexactSearch(W,z):
	CalculateD(W)    # Precalculates lower bounds of the number of mismatches
	# returns the SA intervals of substrings in T that match W with no more than 
	# z differences
	return InexRecur(W, |W| - 1, z, 1, |T| - 1)

CalculateD(W):
	z <- 0
	j <- 0
	for i = 0 to |W| - 1
		if W[j ... i] is not a substring of T
			z <- z + 1
			j <- i + 1
		
		# sets the lower bound of the number of differences in W[0 ... i] to
		# the best match T
		D(i) <- z
	
	return D

InexRecur(W, i, z, k, l):
	if i < 0
		return { k, l }    # for instance a SA interval
	
	if z < D(i)
		return { }
	
	l <- { }
	
	for b in { a, c, g, t }
		k <- C(b) + O(b, k - 1) + 1
		l <- C(b) + O(b, l)
		
		if k <= l
			if b = W[i]
				l <- l union InexRecur(W, i - 1, z, k, l)    # match
			else
				l <- l union InexRecur(W, i - 1, z - 1, k, l)    # mismatch
	
	return l
```
