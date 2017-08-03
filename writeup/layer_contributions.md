*Note:* compile this to html with
```r
templater::render_template("layer_contributions.md")
```

# Contributions of different layers

The covariance between two populations, $i$, and $j$,
can be written as follows.
(Note that this is the *usual* covariance, not the genetic covariance, which subtracts allele frequencies first;
but *after* we have randomly polarized alleles.)
Let $A$ and $B$ be a randomly chosen alleles from the two populations, respectively,
at a randomly chosen site.
(If the two populations are the same, then choose without replacement.)
Then, if we let $X=2(A-1/2)$ and $Y=2(B-1/2)$,
since $X$ and $Y$ take the values $\pm 1$,
and thanks to the random polarization, $\E[A] = \E[A] = 1/2$,
$$\begin{aligned}
    4 \Sigma_{ij} 
        &= \E[XY] \\
        &= \P\{ X=Y \} - \P\{ X \neq Y \} \\
        &= 2 \P\{ X=Y \} - 1 \\
        &= 1 - 2 \P\{ X \neq Y \} .
\end{aligned}$$
To translate, $\P\{ X=Y \}$ is the probability that the alleles from our two focal samples agree with each other,
while $\P\{ X \neq Y \}$ is the probability that they disagree.
For the first thing to happen, $i$ and $j$ must coalesce before coalescing with the reference samples,
so $X$ and $Y$ to agree they must be more closely related than typical in the population.
The last part, $\P\{ XY \neq 0\}$, is the probability that both samples disagree from the references.
(This implies that $4 \Sigma_{ij} = 1 - 2 D_{ij}$, 
where $D_{ij}$ is the probability that two randomly chosen alleles differ, 
which is the genetic divergence.) 

Now, here is a generative model that gives us the form of the covariance we have postulated.
To decide whether or not $X$ and $Y$ will agree,
first each sample randomly chooses a layer: call these $I$ and $J$.
The probability that $X$ chooses layer $k$ is $\P{I=k} = w_i^{(k)}$, and likewise for $Y$.
If they do not choose the same layer, the probability that they agree is $p_\gamma$.
If they do choose the same layer, then they agree with a probability $p_\gamma + q^{(k)}_{ij}$ that depends on their distance apart.
Here
$$\begin{aligned}
    p_\gamma &= (2 (\gamma + \delta_{ij} \eta_i) + 1/2) \\
    q^{(k)}_{ij} &= \left(2 \alpha_0^{(k)} \exp\left( - \left(\alpha_D^{(k)} D_{ij}\right)^{\alpha_2^{(k)}} \right) + 2 \mu^{(k)}  + 1/2 \right) .
\end{aligned}$$

One way to summarize the contribution of each layer is to partition the probability of agreement
into contributions due to agreement "in" each layer.
So, the contribution from layer $k$ to agreement between $i$ and $j$ is $w_i^{(k)} w_j^{(k)} q^{(k)}_{ij} / (p_\gamma + \sum_m w_i^{(m)} w_j^{(m)} q^{(m)}_{ij})$,
which is the probability that they agree thanks to layer $k$ given that they agree.
Similarly, the "background" contribution is $p_\gamma / (p_\gamma + \sum_m w_i^{(m)} w_j^{(m)} q^{(m)}_{ij})$.
*Note:* our signal comes from *variation* in covariance, so maybe we want to omit the $p_\gamma$ terms.

This suggests defining the overall contribution of layer $k$ to agreement
to be the average of that quantity of $i$ and $j$.
Should we include the diagonal?  
Yes, but note that in that case, the $\eta_i$ term is grouped with $\gamma$, and so therefore perhaps omitted.



## Numerical check

The above implies that $\Sigma_{ij} = (1 - 2 D_{ij})/4$, and so $1-D_{ij} = (2 \Sigma_{ij} + 1/2)$,
where $D_{ij}$ is divergence, the proportion of sites at which the two samples differ.
Let's check.
```r
nindivs <- 10
nloci <- 10000
freqs <- rbeta(nloci, 1, 5)
G <- matrix(rbinom(nindivs*nloci, 2, freqs), nrow=nloci)
# randomly polarize
flips <- (runif(nloci)>0.5)
G[flips,] <- 1-G[flips,]

cov_G <- cov(G)
div_G <- cov_G
div_G[] <- NA
for (i in 1:nindivs) {
    for (j in 1:nindivs) {
        div_G[i,j] <- mean(G[,i] * (1-G[,j]) + G[,j] * (1-G[,i]))   
    }
}

layout(t(1:2))
plot(as.vector(cov_G), as.vector(1-2*div_G)/4, main='all entries')
abline(0,1)
plot(as.vector(cov_G)[upper.tri(cov_G,diag=FALSE)], 
     as.vector(1-2*div_G)[upper.tri(cov_G,diag=FALSE)]/4,
     main='offdiags')
abline(0,1)

```
