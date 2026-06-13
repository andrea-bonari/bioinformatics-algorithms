>[!note]
>Many applications (especially in the field of evolutionary genetics/genomics) it's useful to align multiple sequences to each other at the same time.
>
>Formally, a Multiple Sequence Alignment (MSA) of $k$ strings $S_{1},\cdots, S_{k}$ over an alphabet $\Sigma$ is a set of $k$ strings with gaps "$-$" $S_{1}',\cdots,S_{k}'$ over the alphabet $\Sigma'=\Sigma\cup\set{-}$ where: $$|S_{1}'|=\cdots=|S_{k}'|\qquad\text{and}\qquad \forall i\in(1,\cdots,k)\space S_{i}'\text{ is obtained from }S_{i}\text{ with the insertion of gaps}$$
>The possible ways to align $k$ strings depend on where we place gaps.
>
>Given this definition it's always possible to find a "consensus", which is exactly the common ancestor word.
>
>![[Pasted image 20260612085553.png|center]]

Multiple sequence alignment can give us much more information, and a more general view of what the relationships among all our sequences are.

>[!tip] Sum of pairs
>The possible ways to align $k$ strings depend on where we place gaps. Since multiple alignment implicitly induces a pairwise alignment for every pair $S_{i},S_{j}$, we can use the substitution matrix and compute the pairwise alignment for every pair. We can then compute the sums of pair: $$\text{SP}=\sum\limits_{i=1}^{k-1}\sum\limits_{j=i+1}^{k} a(S_{i},S_{j})$$
>The SP score is essentially a generalisation to $k$ strings of the principle we applied to $2$ strings.

### Neighbouring cells
>[!note]
 The optimal solution can also be found by dynamic programming, using a $k$-dimensional hypercube, where the last cell in the structure will contain the score of the optimal multiple alignment.
>
>![[Pasted image 20260612090312.png|center]]
> 
>A cell has $2^{k}-1$.
>For $k=3$, we can then use: $$s_{i,j,k}=\max\begin{cases}s_{i-1,j-1,k-1}+ \delta(v_{i},w_{j}, u_{k}) \\s_{i-1,j-1,k}+\delta(v_{i},w_{j},\_) \\s_{i-1,j,k-1}+\delta(v_{i},\_,w_{j}) \\s_{i,j-1,k-1}+ \delta(\_,w_{j},u_{k}) \\s_{i-1,j,k}+ \delta(v_{i},\_,\_) \\s_{i,j-1,k}+ \delta(\_, w_{j},\_) \\s_{i,j,k-1}+ \delta(\_,\_,u_{k})\end{cases}$$
>
>Where $\delta(x,y,z)$ is an entry in the 3D scoring matrix. Each step is one column in the alignment, so we can compute the corresponding $\delta(x,y,z)$: $$S(A_{i})=\sum\limits_{j=1}^{N}\sum\limits_{k>j} s(a_{i}^{(j)},a_{i}^{(k)})$$
>
>We then need a 2D scoring matrix for alignment matches, to compute $\delta(x,y,z)$ and a gap penalty. Also we set $s(-,-)=0$ because aligning two gaps has little meaning from a biological point of view.

This algorithm has a space complexity of $\mathcal{O}(n^{k})$, and a time complexity of $\mathcal{O}(2^{k},n^{k})$, therefore this problem is NP hard.

>[!tip] NP Hard problems
>When we face a NP hard problem, we have to find a way to solve it anyway. There are different possibilities:
>- Changing the formalisation
>- Explore the space of possible solutions, considering only a reasonable polynomial subset of them, hoping that it contains the optimal solution
>- Design a polynomial time/space algorithm, based on some heuristic hoping to find the optimal solution
>- Design an approximation algorithm
>
>Suppose that we design a polynomial time algorithm for an NP hard maximisation problem. For instance $I$ of the problem, let $A(I)$ be the value of the solution found, while let $M(I)$ be the optimal maximum value for instance $I$.
>If we can prove that: $$\exists\varepsilon\in(0,1]\qquad \text{such that}\qquad\forall I\quad \frac{A(I)}{M(I)}\geq \varepsilon$$
>Then we have a performance guaranteed approximation algorithm.
>
>Of course we want $\varepsilon$ to be as close as possible to $1$.

### Center star algorithm
>[!note]
>Given an input of $k$ strings $S_{1},\cdots,S_{k}$, we compute the pairwise alignment of every pair of input strings. Let $a(S_{i},S_{j})$ be the score of the alignment of strings $S_{i}$ and $S_{j}$, and let $\text{SP}(S_{i})=\sum\limits_{i\neq j}a(S_{i},S_{j})$ be the total pairwise alignment score of $S_{i}$ with all other strings. We then choose $S_{C}$ to be: $$S_{C}=\arg\max_{S_{i}}\text{SP}(S_{i})$$
>We call this the center of the star. The multiple alignment is built by starting with $S_{C}$ and iteratively adding the other strings as they were aligned to $S_{C}$ in step $1$.

Overall, the space/time complexity equals to $\mathcal{O}(k^{2},n^{2})$.

>[!tip] Approximation
>To approximate the Center star algorithm, we just compare all the strings to the center, hoping that the pairwise alignment of each of them would also be found in the optimal multiple alignment.
>
>By minimising the overall distance instead of the maximising the overall similarity, if the triangle inequality holds (for edit distance it does), we can prove that it's an approximation algorithm with $\varepsilon=2$.

### Alignment profile
>[!note]
>One useful way to describe a multiple alignment is to represent it as an alignment profile. For each column of the alignment, the profile represents the frequency with which we find each symbol of the alphabet (including the gap) in that column, where frequencies in each column sum up to $1$.
>
>![[Pasted image 20260612100038.png|center]]

We can add another sequence to this profile by using a dynamic programming matrix, where we replace one of the strings by an entire profile, where the rows represent the string to be added, while the columns have all letters of the alphabet and their respective frequencies.

We then fill the matrix with a "weighted" version of the global alignment rule. 
![[Pasted image 20260612100455.png|center]]

To compute the individual scores for the alignment matrix $A$, given the profile $R$ with frequencies $R[c,j]$ for character $c$ in column $j$, and given a string $X=x_{1}\cdots x_{n}$ to be aligned, we compute the matrix components by using the formula: $$\begin{align*}
A[i,j]&=\max\begin{cases}
A[i-1,j-1]+P(x_{i},j)\qquad&\text{align }x_{i}\text{ to column }j \\
A[i-1,j]+\text{gap}&\text{introduce a gap into the profile} \\
A[i,j-1]+P('\_',j)&\text{introduce a gap into }x
\end{cases}\\

P(x,j)&= \sum\limits_{c\in \Sigma}\text{sim}(x,c)\cdot R[c,j]
\end{align*}$$
Where $\text{sim}(x,c)$ is the cost of aligning $x$ with $c$. Note that $\text{sim}(-,-)=0$.
We then update the alignment profile by adding the new string we just aligned and recomputing the new frequencies.

The same principle can be applied to the alignment of two profiles. Also, for global alignment the update of the cells can be performed with a "double weighted" version of the rule.

### Guide trees
>[!note]
>The idea of progressive alignment is to keep aligning pairs of sequences following their evolutionary history backwards. That is, we start by aligning the strings that are closest relatives, then we move backwards to the more distant history. They key is having the "guide-tree" telling us the order of operations in which the alignment has to be built.
>
>However we need to build the guide tree ourselves.
>
>To do that we start by comparing all sequences and compute their pairwise distances, we choose the pair with minimal distance and build their profile. The alignment represents their common ancestor, and it will be part of the final multiple alignment. We repeat this process iteratively until we get a common ancestor to every string, that way, we constructed the guide tree and the final multiple alignment at the same time. This is actually a hierarchical clustering of the strings.
>![[Pasted image 20260612102615.png|center]]

>[!tip] Hierarchical clustering
>Hierarchical clustering is one of a set of related data mining techniques.
>
>![[Pasted image 20260612102717.png|center]]
>
>Where edge length is often depicted proportionally to the computed distance, showing how "close" the individual items and clusters are.
>
>Clustering is a very important branch for data mining, machine learning and statistics, it can be applied to different kinds of data, and each item to be clustered is seen as a point in a multi dimensional space. The data will then be measured with a different distance (or similarity) measures.
>

In the guide tree, the distance between two strings is the distance given by their alignment, while the distance between one string and an alignment can be measured by:
- Distance between the string and the closes string of the cluster
- Distance between the string and the furthest string of the cluster
- The average distance between the string and all the strings of the cluster
While, for the distance between two clusters:
- Distance of the closest pair of members
- Distance of the furthest pair of members
- Average distance between all pairs of members

For multiple sequence alignment we use neighbour joining, where the measure employed combines the average distance between two objects and the distance of the two objects from all the other objects.

### Local multiple alignment
>[!note]
>Given a set of input strings $S_{1},\cdots,S_{k}$, a fixed substring length $w$, and a function $f(x)$ for the evaluation of the multiple alignment, we need to find the set of $k$ substrings of length $w$, one per input string, that produces the alignment of maximum score according to function $f(x)$.
>
>This problem is also known as "motif finding".

For motif finding we consider, for each input string we take exactly one substring, all substrings which are part of the solution have the same length $k$, and there will be no gaps in the final solution. Therefore any candidate solution can be simply described by a vector of integers $(p_{1},\cdots,p_{k})$.

We can describe a solution by means of the corresponding profile, with one row per symbol of the alphabet $\Sigma$ and $f_{i,j}$ the frequency of symbol $i$ in column $j$. We can then calculate: $$S(P)=\sum\limits_{i=1}^{|\Sigma|}\sum\limits_{j=1}^{w}f_{i,j}\log_{2}\frac{f_{i,j}}{\frac{1}{|\Sigma|}}$$
We then align every substring of string $S_{i}$ to every substring of $S_{i+1}$, we obtain a $\mathcal{O}(n^{2})$ profiles, but we only keep the best $n$ profiles, according to their score. Repeating this iteratively gives us a $\mathcal{O}(n^{2})$.

### Combinatorial optimisation
>[!note]
>A solution to the problem is an assignment of a value between $1$ and $(n-w+1)$ to $k$ discrete variables $X_{1},\cdots, X_{k}$ (assuming all strings have length $n$).
>
>Hence, each candidate solution to the problem can be seen as a point in a $k$-dimensional "search space" of the problem.
>
>In out case there are $(n-w+1)^{k}$ points that can be visited in the search space.
>
>![[Pasted image 20260612131755.png|center]]

>[!tip] Local search
>We start at a point in the $k$-dimensional space, chosen randomly or with some heuristic. We then look around at the neighbouring points and compute their scores, and choose the point that will bring the best improvement. We keep moving like this until we find a point at which no further improvement of the objective function is possible. This point is a local maximum.
>
>"Looking around" in the search space at neighbouring points corresponds to changing the value of one or more of the variables in the current solution.
>One possible way is to select one variable $X_{i}$, keeping the other variables fixed and choose the value $X_{i}$ that brings the best improvement. We iterate this step, by cycling over the $X_{i}$ variables. This is equivalent to moving in one direction in the $k$-dimensional space, going to the "maximum" point we see in that direction.
>
>This technique is known as "hill climbing".
>
>![[Pasted image 20260612132420.png|center]]
>
>A good choice is to iterate the whole algorithm many times with different starting points and output the best solution found across all searches.
>
>![[Pasted image 20260612132518.png|center]]
>
>Using an approach similar to MAX-CUT we can:
>1. Select an $X_{i}$
>2. Remove the substring corresponding to $X_{i}$ from the current solution
>3. For each $j\in[1,n-w+1]$, we compute the score by setting $X_{i}=j$ and leaving all other variables unchanged
>4. Choose $j$ such that $X_{i}=j$ corresponds to the profile with the highest score
>5. Adjust the alignment accordingly and select the next $X_{i}$ to update
>
>We usually process the variables in random order, using a different random permutation at each cycle, therefore we might get different results even starting from the same initial point. Therefore we iterate several times with different starting profiles.

>[!tip] MAX-CUT
>A famous NP hard problem for combinatorial optimisation is the MAX-CUT problem. Given a graph $(V,E)$, find the partition of the vertices/nodes into two subsets $V_{1},V_{2}$ with $V_{1}\cup V_{2}=V$ such that the number of edges connecting one node of the first partition $V_{1}$ with one node of the second partition $V_{2}$ is maximised.
>
>![[Pasted image 20260612133008.png|center]]
>
>Thus, given $k$ vertices, each solution can be described by $k$ binary variables $X_{1},X_{2},\cdots, X_{k}$ indicating to which of the two subsets each node belongs. Therefore $f(X_{1},\cdots, X_{k})$ is the objective function to be optimised.
>
>To solve this we initialise each binary variable $X_{i}$ randomly, and variables are updated one by one, from $X_{1}$ to $X_{k}$.
>
>We the iteratively calculate $f(x_{1},\cdots, \overline{x}_{i},\cdots, x_{k})$, and if it's greater then the previous solution it becomes the new solution. This is repeated over all the variables $X_{i}$ until a complete cycle over all the variables is completed without any improvement.

>[!tip] Stochastic optimisation
>In stochastic search, given $X_{1},\cdots, X_{k}$ binary variables, at step $i$, we compare the value $X_{i}$ with the logical not of its current value, and we the compute: $$\Delta f= f(x_{1},\cdots, \overline{x}_{i},\cdots, x_{k})-f(x_{1},\cdots, x_{i},\cdots, x_{k})$$
>If $\Delta f>0$ we move to the next solution, if not we might accept the change and update the current solution with probability: $$P[\text{update }X_{i}]= e^{\frac{\Delta t}{T}}$$
>Where $T$ is a parameter that allows for "fine tuning" the degree of randomness. The value of $T$ can be dynamically changed during the execution of the algorithm, for example, we can start with a high value and decrease it at each iteration (thus making the algorithm more and more deterministic).
>
>Applying it to our problem, we have $k$ variables $X_{1},\cdots, X_{k}$ defining the substrings' starting points in the original input, and at step $t$, we choose to replace the value of the variable $X_{i}$. For all possible starting positions $p_{j}$ that the variable can assume we consume the $S(p_{j})$ score for each possible $p_{j}$, and then we compute the sum of all the scores for all possible substrings for $X_{i}$: $$S(S_{i})=\sum\limits_{j=1}^{n-w+1}S(p_{j})$$
>We accept a new value for variable $X_{i}$ with a random choice of probability: $$\text{Pr}[X_{i}=p_{j}]= \frac{S(p_{j})}{S(S_{i})}$$
>This strategy is known as Gibbs sampling.