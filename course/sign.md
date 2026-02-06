<head>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.27/dist/katex.min.css" integrity="sha384-Pu5+C18nP5dwykLJOhd2U4Xen7rjScHN/qusop27hdd2drI+lL5KvX7YntvT8yew" crossorigin="anonymous">
    <!-- The loading of KaTeX is deferred to speed up page rendering -->
    <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.27/dist/katex.min.js" integrity="sha384-2B8pfmZZ6JlVoScJm/5hQfNS2TI/6hPqDZInzzPc8oHpN5SgeNOf4LzREO6p5YtZ" crossorigin="anonymous"></script>
    <!-- To automatically render math in text elements, include the auto-render extension: -->
    <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.27/dist/contrib/auto-render.min.js" integrity="sha384-hCXGrW6PitJEwbkoStFjeJxv+fSOOQKOPbJxSfM6G5sWZjAyWhXiTIIAmQqnlLlh" crossorigin="anonymous"
        onload="renderMathInElement(document.body);"></script>
    <link rel="stylesheet" type="text/css" href="https://tikzjax.com/v1/fonts.css">
    <script src="https://tikzjax.com/v1/tikzjax.js"></script>
</head>

# The Annotated ML DSA Algorithms
```{=latex}
\newcommand{\TComment}[1]{\qquad #1}
\newcommand{\LComment}[1]{\text{---} #1}
```

## Conversion Between Data Types
```{=latex}
\begin{algorithm}
\caption*{\textbf{Algorithm 14} \text{CoeffFromThreeBytes}($b_0: \mathbb{B}^1, b_1: \mathbb{B}^1, b_2: \mathbb{B}^1) \rightarrow z: Z_q \cup \bot$}
\vspace*{2pt}
\normalsize{Generates an element of $\{0, 1, 2, \ldots, q-1\} \cup \bot$.}

\begin{algorithmic}[1]
\State $b'_2 \gets b_2$
\If {$b'_2 > 127$}
  \State $b'_2 \gets b'_2 - 128$ \Comment{set the top bit of $b'_2$ to zero}
\EndIf
\State {$z \gets 2^{16} \cdot b'_2 + 2^8 \cdot b_1 + b_0$} {\Comment{$0 \leq z \leq z^{23}-1}$}
    \Statex {\TComment{$\triangleright \ \mathrm{max}(z) = 2^{16} \cdot (2^7 - 1) + 2^8 \cdot (2^8 - 1) + (2^8 - 1)$}}
    \Statex {\TComment{$\triangleright \hspace*{1.4cm}  = 2^{23} - 1$}}
\If {$z < q$} {\Return z} {\Comment{rejection sampling}}
\Else \ \Return \bot
\EndIf
\end{algorithmic}
\end{algorithm}
```


```{=latex}
\begin{algorithm}
\caption*{\textbf{Algorithm 15} \text{CoeffFromHalfByte}($b: [0, 15]) \rightarrow z: [-\eta, \eta] \cup \bot$}
\vspace*{2pt}
\normalsize{$\mathrm{Let} \ \eta \in \{2, 4\}. \text{Generates an element of} \ \{-\eta, -\eta+1, \ldots, \eta\} \cup \{\bot\}$.}

\begin{algorithmic}[1]
\If {$\eta = 2 \ \mathbf{and} \ b < 15$} {\Return {$2 - (b \ \mathrm{mod} \ 5)$}}
    \Comment{rejection sampling from $\{−2, \ldots, 2\}$}
    \Statex {\TComment{$\triangleright$ \ case 1: ML-DSA-44 and ML-DSA-87}}
\Else
    \Statex {\TComment{$\triangleright$ \ case 2: ML-DSA-65}}
    \If {$\eta = 4 \ \mathbf{and} \ b < 9$} {\Return $4 - b$} \Comment{rejection sampling from $\{−4, \ldots, 4\}$}
    \Else \ {\Return $\bot$}
    \EndIf
\EndIf
\end{algorithmic}
\end{algorithm}
```


## Pseudorandom Sampling


```{=latex}
\begin{algorithm}
\caption*{\textbf{Algorithm 29} \text{SampleInBall}($\rho: \mathbb{B}^{\lambda/4}$) $\rightarrow c: R_{[-1, 0, 1]}$}
\vspace*{2pt}
\normalsize{Samples a polynomial $c \in R$ with coefficients from $\{-1, 0, 1\}$ and \textit{Hamming weight} $\tau \le 64$.}

\begin{algorithmic}[1]
\State $c \gets 0$
\State $\mathrm{ctx} \gets \mathrm{H.Init()}$ \Comment{\text{H $\doteq$ SHAKE256}}
\State $\mathrm{ctx} \gets \mathrm{H.Absorb}(\mathrm{ctx}, \rho)$
    \Statex \TComment{\text{$\triangleright$ $len(\rho)$ is 32 bytes in ML-DSA-44, 48 in ML-DSA-65, and 64 in ML-DSA-87.}}
\State $(\mathrm{ctx}, s:\mathbb{B}^8) \gets \mathrm{H.Squeeze}(\mathrm{ctx}, 8)$
\State $h: \{0, 1\}^{64} \gets \mathrm{BytesToBits}(s)$
    \Comment {$h$ is a bit string of length 64}
\For{$i$ from $256 - \tau$ to 255}
    \Statex {\TComment{$\triangleright \ \tau = 39, 49, 60$ in ML-DSA-44, ML-DSA-65, and ML-DSA-87, respectively.}}
    \State $(\mathrm{ctx}, \, j: \mathbb{B}^1) \gets \mathrm{H.Squeeze(ctx, 1)}$
    \While {$j > i$} \Comment{rejection sampling in $\{0, \ldots, i\}$}
        \State $(\mathrm{ctx}, \, j: \mathbb{B}^1) \gets \mathrm{H.Squeeze(ctx, 1)}$
    \EndWhile \Comment{$j$ is a pseudorandom byte that is $\le i$}
    \State {$c_i \gets c_j$}
        \Statex {\TComment{$\triangleright \ c_j$ is a smaller-degree coefficient, pseudorandomly selected.}}
        \Statex {\TComment {$\triangleright \ c_i$ is a larger degree coefficient. $c_i$ receives the value of $c_j$.}}
    \State {$c_j \gets (-1)^{h[i+\tau-256]}$}
        \Statex {\TComment{$\triangleright$ access pattern: $h[0], h[1], \ldots, h[\tau] \ where \ 39 \le \tau \le 60, and \ h:\{0, 1\}^{64}$.}}
        \Statex {\TComment{$\triangleright$ This pseudorandom shuffling is performed $\tau$ times in total.}}
\EndFor
\State \Return c

\end{algorithmic}
\end{algorithm}
```


```{=latex}
\begin{algorithm}
\caption*{\textbf{Algorithm 30} \text{RejNTTPoly}($\rho: \mathbb{B}^{34}$) $\rightarrow \hat{a}: T_q$}
\vspace*{2pt}
\normalsize{Samples a polynomial $\hat{a} \in T_q$.}
\begin{algorithmic}[1]
\State $j \gets 0$
\State $\mathrm{ctx} \gets \mathrm{G.Init()}$ \Comment{\text{G $\doteq$ SHAKE128}}
\State $\mathrm{ctx} \gets \mathrm{G.Absorb}(\mathrm{ctx}, \rho)$
\While{$j < 256$}
    \State $(\mathrm{ctx}, s:\mathbb{B}^3) \gets \mathrm{G.Squeeze}(\mathrm{ctx}, 3)$
    \State $\hat{a}_j \gets \mathrm{CoeffFromThreeBytes}(s[0], s[1], s[2])$
    \If{$\hat{a}[j] \neq \bot$}
    \State $j \gets j+1$
    \EndIf
\EndWhile

\State \Return $\hat{a}$

\end{algorithmic}
\end{algorithm}
```

```{=latex}
\begin{algorithm}
\caption*{\textbf{Algorithm 31} \text{RejBoundedPoly}($\rho: \mathbb{B}^{66}$) $\rightarrow a: R_{[-\eta, \eta]}$}
\vspace*{2pt}
\normalsize{Samples an element $a \in R$ with coefficients in $[-\eta, \eta]$ computed via rejection sampling from $\rho$.}
\begin{algorithmic}[1]
\State $j \gets 0$
\State $\mathrm{ctx} \gets \mathrm{H.Init()}$ \Comment{\text{H $\doteq$ SHAKE256}}
\State $\mathrm{ctx} \gets \mathrm{H.Absorb}(\mathrm{ctx}, \rho)$

\While{$j < 256$}
    \State {$z: \mathbb{B}^1 \gets \mathrm{H.Squeeze(ctx, 1)}$}
    \State {$z_0 \gets \mathrm{CoeffFromHalfByte}(z \mathrm{mod 16})$}
    \State {$z_1 \gets \mathrm{CoeffFromHalfByte}(z/16)$}
    \If {$z_0 \neq \bot$}
        \State {$a_j \gets z_0$}
        \State {$j \gets j+1$}
    \EndIf
    \If {$z_1 \neq \bot \ \mathbf{and} \ j < 256 $}
        \State {$a_j \gets z_1$}
        \State {$j \gets j+1$}
    \EndIf
\EndWhile

\State \Return $a$

\end{algorithmic}
\end{algorithm}
```



## The Key Generation Algorithm
```{=latex}
\begin{algorithm}
\caption*{\textbf{Algorithm 6} \text{ML-DSA.KeyGen\_internal}($\xi$)}
\begin{algorithmic}[1]
\State {$(\rho, \rho', K) \gets \text{H}(\xi \, \Vert \, \text{Integer2Bytes}(k, 1) \, \Vert \, \text{Integer2Bytes}(\ell, 1), 128)$}
\State {\Comment{expand seed}}
\State {$\hat{\textbf{A}} \gets \text{ExpandA}(\rho)$}
\State {$(\textbf{s}_1, \textbf{s}_2) \gets \text{ExpandS}(\rho')$}
\State {$t \gets \text{NTT}^{-1}(\hat{\textbf{A}} \circ \text{NTT}(\textbf{s}_1)) + \textbf{s}_2$}
\State {$(\textbf{t}_1, \textbf{t}_0) \gets \text{Power2Round}(\textbf{t})$} \Comment{compress t}
    \State \Comment{Power2Round is applied componentwise}
\State {$pk \gets \text{pkEncode}(\rho, \textbf{t}_1)$}
\State {$tr \gets \text{H}(pk, 64)$}
\State {$sk \gets \text{skEncode}(\rho, K, tr, \textbf{s}_1, \textbf{s}_2, \textbf{t}_0)$}
\State {$\textbf{return} \ (pk, sk)$}
\end{algorithmic}
\end{algorithm}
```

## The Annotated Key Generation Algorithm
```{=latex}
% \newcommand{\TComment}[1]{\qquad #1}
\begin{algorithm}
\caption*{\textbf{Annotated Algorithm 6} \text{ML-DSA.KeyGen\_internal}($\xi$)}
\begin{algorithmic}[1]
% line 1
\State {$(\rho: \mathbb{B}^{32}, \rho':\mathbb{B}^{64}, K: \mathbb{B}^32) \gets \text{H}(\xi \, \Vert \, \text{Integer2Bytes}(k, 1) \, \Vert \, \text{Integer2Bytes}(\ell, 1), 128)$}
% line 2
\State {\Comment{expand seed}}
% line 3
\State {$\hat{\textbf{A}}:T_q^{k \times \ell} \gets \text{ExpandA}(\rho)$}
% line 4
\State {$(\textbf{s}_1:R_m^{\ell}, \textbf{s}_2:R_m^{k}) \gets \text{ExpandS}(\rho')$ }
    \Statex \TComment{$ m \in [ -2, 2 ] \ \textit{if} \ \text{ML-DSA-44} \ \textit{or} \ \text{ML-DSA-87} $} \Comment{$\eta = 2$}
    \Statex \TComment{$ m \in [ -4, 4 ] \ \textit{otherwise} $} \Comment{$\eta = 4$ }

% line 5
\State {$t:R_q^k \gets \text{NTT}^{-1}(\hat{\textbf{A}} \circ \text{NTT}(\textbf{s}_1)) + \textbf{s}_2$ }

% line 6
\State {$(\textbf{t}_1: R_{q_1}^k, \textbf{t}_0: R_{q_0}^k) \gets \text{Power2Round}(\textbf{t})$}
    \Comment{compress t }
    \Statex \TComment{$ t_{q_1} \in [0, 1023] $} \Comment{10-bit value }
    \Statex \TComment{$ t_{q_0} \in [ -4095, 4096 ] $} \Comment{$\text{mod}^{\pm} 2^d, \, d = 13, \, ([ -2^{12}+1, 2^{12} ])$}
% line 7
\State \Comment{Power2Round is applied componentwise }

% line 8
\State {$pk \gets \text{pkEncode}(\rho, \textbf{t}_1)$}
% line 9
\State {$tr: \mathbb{B}^{64} \gets \text{H}(pk, 64)$}
% line 10
\State {$sk \gets \text{skEncode}(\rho, K, tr, \textbf{s}_1, \textbf{s}_2, \textbf{t}_0)$}
% line 11
\State {$\textbf{return} \ (pk, sk)$}

\end{algorithmic}
\end{algorithm}
```


## The Annotated Signature Algorithm

The ML-DSA Signing (Internal) algorithm named \textbf{ML-DSA.Sign\_internal} in FIPS 204 standard is reproduced below for reference. The line numbers and the pseudocode matches exactly with the original version.

```{=latex}
\begin{algorithm}
\caption*{\textbf{Algorithm 7} \text{ML-DSA.Sign\_internal}(\textit{sk, M', rnd})}
\begin{algorithmic}[1]
\State $(\rho, K, tr, s_1, s_2, t_0) = \text{skDecode}(sk)$
\State $\hat{\textbf{s}}_1 \gets \text{NTT}(s_1)$
\State $\hat{\textbf{s}}_2 \gets \text{NTT}(s_2)$
\State $\hat{\textbf{t}}_0 \gets \text{NTT}(t_0)$
\State $\hat{\textbf{A}} \gets \text{ExpandA}(\rho)$
\State $\mu \gets \text{H}(\text{BytesToBits}(tr \,\Vert\, M', 64))$
\State $\rho'' \gets \text{H}(K \,\Vert\, rnd \,\Vert\, \mu, 64)$
\State $\kappa \gets 0$
\State $(\textbf{z}, \textbf{h}) \gets \bot$

\While{(\textbf{z}, \textbf{h}) = $\bot$}
    \State $\text{y} \in R_q^{\ell} \gets \text{ExpandMask}(\rho'', \kappa)$
    \State $\text{w} \gets \text{NTT}^{-1}(\hat{\textbf{A}} \circ \text{NTT}(\text{y}))$
    \State $\text{w}_1 \gets \text{HighBits}(\text{w})$
        \State \Comment{HighBits is applied componentwise}
    \State $\tilde{c} \gets \text{H}(\mu \,\Vert\, \text{w1Encode}(\text{w}_1), \lambda/4)$
    \State $c \in R_q \gets \text{SampleInBall}(\tilde{c})$
    \State $\hat{c} \gets \text{NTT}(c)$
    \State $\langle\!\langle c\textbf{s}_1 \rangle\!\rangle \gets \text{NTT}^{-1}(\hat{c} \circ \hat{\textbf{s}}_1)$
    \State $\langle\!\langle c\textbf{s}_2 \rangle\!\rangle \gets \text{NTT}^{-1}(\hat{c} \circ \hat{\textbf{s}}_2)$
    \State $\textbf{z} \gets \textbf{y} + \langle\!\langle c\textbf{s}_1 \rangle\!\rangle$
    \State $\textbf{r}_0 \gets \text{LowBits}(\textbf{w} - \langle\!\langle c\textbf{s}_2 \rangle\!\rangle)$
        \State \Comment{LowBits is applied componentwise}
    \If{$ \Vert \textbf{z} \Vert_\infty \ge \gamma_1 - \beta \ \, \textbf{or} \ \, \Vert \textbf{r}_0 \Vert_\infty \ge \gamma_2 - \beta  $}
        $(\textbf{z}, \textbf{h}) \gets \bot$
    \Else
        \State $\langle\!\langle c\textbf{t}_0 \rangle\!\rangle \gets \text{NTT}^{-1}(\hat{c} \circ \hat{\text{t}}_0)$
        \State $\textbf{h} \gets \text{MakeHint}(-\langle\!\langle c\textbf{t}_0\rangle\!\rangle, \textbf{w} - \langle\!\langle c\textbf{s}_2\rangle\!\rangle + \langle\!\langle c\textbf{t}_0\rangle\!\rangle )$
            \State \Comment{MakeHint is applied componentwise}
        \If{$\Vert \langle\!\langle c\textbf{t}_0\rangle\!\rangle \Vert_\infty \ge \gamma_2 \ \, \textbf{or} \ \, \text{the number of 1's in} \ \textbf{h} \ \text{is greater than} \ \omega$} $(\textbf{z}, \textbf{h}) \gets \bot $
        \EndIf
    \EndIf
    \State $\kappa \gets \kappa + \ell$
 \EndWhile

 \State $\sigma \gets \text{sigEncode}(\tilde{c}, \textbf{z} \ \text{mod}^{\pm} \ q, \textbf{h})$
 \State \Return $\sigma$

\end{algorithmic}
\end{algorithm}
```

```{=latex}
\begin{algorithm}
\caption*{\textbf{Annotated Algorithm 7} \text{ML-DSA.Sign\_internal}(\textit{sk, M', rnd})}
\begin{algorithmic}[1]
\State $(\rho:\mathbb{B}^{32}, \ K:\mathbb{B}^{32}, \ tr:\mathbb{B}^{64}, \ s_1:R_m^{\ell}, \ s_2:R_m^{k}, \ t_0:R_t^k) = \text{skDecode}(sk)$
    \Statex \TComment{$ m \in [ -2, 2 ] \ \textit{if} \ \text{ML-DSA-44} \ \textit{or} \ \text{ML-DSA-87} $} \Comment{$\eta = 2$}
    \Statex \TComment{$ m \in [ -4, 4 ] \ \textit{otherwise} $} \Comment{$\eta = 4$}
    \Statex \TComment{$ t \in [ -2^{12}+1, 2^{12}-1 ] $} \Comment{$d = 13$}
\State $\hat{\textbf{s}}_1: T_q^{\ell} \gets \text{NTT}(s_1)$
\State $\hat{\textbf{s}}_2: T_q^k \gets \text{NTT}(s_2)$
\State $\hat{\textbf{t}}_0: T_q^k \gets \text{NTT}(t_0)$
\State $\hat{\textbf{A}}: T_q^{k \times \ell} \gets \text{ExpandA}(\rho)$
\State $\mu:\mathbb{B}^{64} \gets \text{H}(\text{BytesToBits}(tr \,\Vert\, M', 64))$
\State $\rho'':\mathbb{B}^{64} \gets \text{H}(K \,\Vert\, rnd \,\Vert\, \mu, 64)$
\State $\kappa: \text{u16} \gets 0$ \Comment{unsigned 16-bit integer}

\State $(\textbf{z}: R_{mod^{\pm} q}^l, \ \textbf{h}: R_2^k) \gets \bot$

\While{(\textbf{z}, \textbf{h}) = $\bot$} \Comment{Rejection sampling loop}
    \State {$\textbf{y}: R_q^{\ell} \gets \text{ExpandMask}(\rho'', \kappa)$} \Comment{Sample a fresh mask $y$ mixing $\rho''$ and new $\kappa$ }

    \State $\textbf{w}: R_q^k \gets \text{NTT}^{-1}(\hat{\textbf{A}} \circ \text{NTT}(\textbf{y}))$

    \State $\textbf{w}_1: R_{m_1}^k \gets \text{HighBits}(\textbf{w})$
        \Statex \TComment{$ m_1 \in [ 0, 21 ] \ \textit{if} \ \text{ML-DSA-44} $}
        \Statex \TComment{$ m_1 \in [ 0, 15 ] \ \textit{if} \ \text{ML-DSA-65} \ \textit{or} \ \text{ML-DSA-87}$}
        \State \Comment{\text{HighBits is applied componentwise}}

    \State $\tilde{c}:\mathbb{B}^{\lambda/4} \gets \text{H}(\mu \,\Vert\, \text{w1Encode}(\textbf{w}_1), \lambda/4)$
        \Comment{$ \lambda \ \text{is the collission strength of} \ \tilde{c}$}
        \Statex \TComment{$\triangleright \ \lambda/4$ = 32 (ML-DSA-44) or 48 (ML-DSA-65) or 64 (ML-DSA-87)}

    \State $c: R_q \gets \text{SampleInBall}(\tilde{c})$
        \Comment{polynomial in $R$ with coefficients from \{-1, 0, 1\} }

    \State $\hat{c}: T_q \gets \text{NTT}(c)$

    \State $\langle\!\langle c\textbf{s}_1 \rangle\!\rangle: R_q^{\ell} \gets \text{NTT}^{-1}(\hat{c} \circ \hat{\textbf{s}}_1)$ \Comment{$\hat{c}$ is multiplied with each of $\hat{\textbf{s}}_1^0, \ldots, \hat{\textbf{s}}_1^{l-1}$}

    \State $\langle\!\langle c\textbf{s}_2 \rangle\!\rangle: R_q^k \gets \text{NTT}^{-1}(\hat{c} \circ \hat{\textbf{s}}_2)$ \Comment{$\hat{c}$ is multiplied with each of $\hat{\textbf{s}}_2^0, \ldots, \hat{\textbf{s}}_2^{k-1}$}

    \State $\textbf{z}: R_{mod^{\pm} q}^l \gets \textbf{y} + \langle\!\langle c\textbf{s}_1 \rangle\!\rangle$
        \Comment{$\textbf{z}: R_{mod^{\pm} q}^l \gets \textbf{y}: R_q^{\ell} + \langle\!\langle c\textbf{s}_1 \rangle\!\rangle: R_q^{\ell} $}

    \State $\textbf{r}_0: R_q^k \gets \text{LowBits}(\textbf{w} - \langle\!\langle c\textbf{s}_2 \rangle\!\rangle)$
    \State \Comment{\text{LowBits is applied componentwise}}

    \If{$ \Vert \textbf{z} \Vert_\infty \ge \gamma_1 - \beta \ \, \textbf{or} \ \, \Vert \textbf{r}_0 \Vert_\infty \ge \gamma_2 - \beta $}
        $(\textbf{z}, \textbf{h}) \gets \bot$
    \Else
        \State $\langle\!\langle c\textbf{t}_0 \rangle\!\rangle: T_q^k \gets \text{NTT}^{-1}(\hat{c} \circ \hat{\textbf{t}}_0)$
        \Comment{$\hat{c}$ is multiplied with each of $\hat{\textbf{t}}_0^0, \ldots, \hat{\textbf{t}}_0^{k-1} $}

        \State $\textbf{h}: R_2^k \gets \text{MakeHint}(-\langle\!\langle c\textbf{t}_0\rangle\!\rangle, \textbf{w} - \langle\!\langle c\textbf{s}_2\rangle\!\rangle + \langle\!\langle c\textbf{t}_0\rangle\!\rangle )$
        \State \Comment{MakeHint is applied componentwise  }

        \If{$\Vert \langle\!\langle c\textbf{t}_0\rangle\!\rangle \Vert_\infty \ge \gamma_2 \ \, \textbf{or} \ \, \text{the number of 1's in} \ \textbf{h} \ \text{is greater than} \ \omega$} $(\textbf{z}, \textbf{h}) \gets \bot $
        \EndIf
    \EndIf
    \State $\kappa \gets \kappa + \ell$
 \EndWhile
 \State $\sigma \gets \text{sigEncode}(\tilde{c}, \textbf{z} \ \text{mod}^{\pm} q, \textbf{h})$

 \State \Return $\sigma$

\end{algorithmic}
\end{algorithm}
```
