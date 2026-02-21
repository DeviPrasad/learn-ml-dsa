# The Annotated ML DSA Algorithms

## Conversion Between Data Types
```{=latex}
\begin{algorithm}[H]
\caption*{\textbf{Algorithm 14}
    $\CoeffFromThreeBytes(\bZero: \OneByte, \bOne: \OneByte, \bTwo: \OneByte)
                                                            \rightarrow z: Z_q \cup \bot $ }
\vspace*{2pt}
\normalsize{Generates an element of $\{0, 1, 2, \ldots, q-1\} \cup \bot$.}
\begin{algorithmic}[1]
%line 1
\State {$ \bTwoDash \is \bTwo $}

%line 2
\If {$ \bTwoDash > 127 $}

%line 3
    \State {$ \bTwoDash \is \bTwoDash - 128 $}
    \Comment{set the top bit of $\bTwoDash$ to zero}

%line 4
\EndIf

%line 5
\State {$z \is 2^{16} \cdot \bTwoDash + 2^8 \cdot \bOne + \bZero $}
    {\Comment{$ 0 \leq z \leq z^{23}-1} $}
    \Statex {\TComment\LComment{$\mathrm{max}(z) = 2^{16} \cdot (2^7 - 1) + 2^8 \cdot (2^8 - 1) + (2^8 - 1)$}}
    \Statex {\TComment\LComment{\hspace*{1.3cm} = $2^{23} - 1$ }}

%line 6
\If {$z < q$} {\Return z} {\Comment{rejection sampling}}

%line 7
\Else \ \Return \bot

%line 8
\EndIf

\end{algorithmic}
\end{algorithm}
```


```{=latex}
\begin{algorithm}[H]
\caption*{\textbf{Algorithm 15} $\CoeffFromHalfByte(b: [0, 15]) \rightarrow z: [-\eta, \eta] \cup \bot $}
\vspace*{2pt}
{Generates an element of $\{-\eta, -\eta+1, \ldots, \eta\} \cup \{\bot\}. $}\\
{$\eta \in \{2, 4\}. $}

\begin{algorithmic}[1]
%line 1
\If {$ \eta = 2 \ \mathbf{and} \ b < 15 $} {\Return {$ 2 - (b \ \Mod{5}) $} }
    \Comment{rejection sampling from $\{−2, \ldots, 2\}$}
    \Statex {\TComment\LComment{ case 1: ML-DSA-44 and ML-DSA-87 }}
%line 2
\Else
    \Statex {\TComment\LComment{ case 2: ML-DSA-65 }}
%line 3
    \If {$\eta = 4 \ \mathbf{and} \ b < 9$} {\Return $4 - b$}
        \Comment{rejection sampling from $\{−4, \ldots, 4\}$}
%line 4
    \Else \ {\Return $\bot$}
%line 5
    \EndIf
%line 6
\EndIf

\end{algorithmic}
\end{algorithm}
```


## Pseudorandom Sampling

```{=latex}
\newcommand{\Ctx}{\mathrm{ctx}}
\newcommand{\HInit}{\mathrm{H.Init}}
\newcommand{\HAbsorb}[2]{\mathrm{H.Absorb}(#1, #2)}
\newcommand{\HSqueeze}[2]{\mathrm{H.Squeeze}(#1, #2)}
\newcommand{\SixtyFourBits}{\ensuremath{\{0, 1\}^{64}}}
\newcommand{\Ci}{\ensuremath{c_i}}
\newcommand{\Cj}{\ensuremath{c_j}}

\begin{algorithm}[H]
\caption*{\textbf{Algorithm 29} $\SampleInBall(\rho: \mathbb{B}^{\lambda/4}) \rightarrow c: R_{[-1, 0, 1]} $}
\vspace*{2pt}
{Samples a polynomial $c \in R$ with coefficients from $\{-1, 0, 1\}$ and \textit{Hamming weight} $\tau \le 64$.}

\begin{algorithmic}[1]

%line 1
\State $c \is 0$

%line 2
\State {$ \mathrm{ctx} \is \HInit() $}
    \Comment{$\NH \doteq$ SHAKE256}

%line 3
\State {$ \Ctx \is \HAbsorb{\Ctx}{\rho} $}
    \Statex \TComment\LComment{$len(\rho)$ is 32 bytes in ML-DSA-44, 48 in ML-DSA-65, and 64 in ML-DSA-87.}

%line 4
\State {$ (\Ctx, s:\BytesEight) \is \HSqueeze{\Ctx}{8} $}

%line 5
\State {$h: \SixtyFourBits \is \BytesToBits(s) $}
    \Comment {$h$ is a bit string of length 64}

%line 6
\For{$i$ from $256 - \tau$ to 255}
    \Statex {\TComment\LComment{$\tau = 39, 49, 60$ in ML-DSA-44, ML-DSA-65, and ML-DSA-87, respectively.}}

%line 7
    \State {$ (\Ctx, \, j: \OneByte) \is \HSqueeze{\Ctx}{1} $}

%line 8
    \While {$j > i$} \Comment{rejection sampling in $\{0, \ldots, i\}$}

%line 9
        \State {$ (\Ctx, \, j: \OneByte) \is \HSqueeze{\Ctx}{1} $}

%line 10
    \EndWhile \Comment{$j$ is a pseudorandom byte that is $\le i$}

%line 11
    \State {$ \Ci \is \Cj$}
        \Statex {\TComment\LComment{$\Cj$ is a smaller-degree coefficient, pseudorandomly selected.}}
        \Statex {\TComment\LComment{$\Ci$ is a larger degree coefficient. $\Ci$ receives the value of $Cj$.}}

%line 12
    \State {$ \Cj \is (-1)^{h[i+\tau-256]} $}
        \Statex {\TComment\LComment{access pattern: $h[0], h[1], \ldots, h[\tau] \ where \ 39 \le \tau \le 60$, \textit{and} $h:\{0, 1\}^{64}$.} }
        \Statex {\TComment\LComment{ This pseudorandom shuffling is performed $\tau$ times in total.}}

%line 13
\EndFor

%line 14
\State \Return c

\end{algorithmic}
\end{algorithm}
```


```{=latex}
\newcommand{\GInit}{\mathrm{G.Init}}
\newcommand{\GAbsorb}[2]{\mathrm{G.Absorb}(#1, #2)}
\newcommand{\GSqueeze}[2]{\mathrm{G.Squeeze}(#1, #2)}

\begin{algorithm}[H]
\caption*{\textbf{Algorithm 30} \text{RejNTTPoly}($\rho: \mathbb{B}^{34}$) $\rightarrow \hat{a}: T_q$}
\vspace*{2pt}
Samples a polynomial $\HA \in T_q$.

\begin{algorithmic}[1]
%line 1
\State {$ j \is 0 $}
%line 2
\State {$ \Ctx \is \GInit() $} \Comment{G $\doteq$ SHAKE128}
%line 3
\State {$ \Ctx \is \GAbsorb{\Ctx}{\rho} $}
%line 4
\While {$ j < 256 $}
%line 5
    \State {$ (\Ctx, s:\ThreeBytes) \is \GSqueeze{\Ctx}{3} $}
%line 6
    \State {$ \HA[j] \is \CoeffFromThreeBytes(s[0], s[1], s[2]) $}
%line 7
    \If {$ \HA[j] \neq \bot $}
%line 8
    \State {$ j \is j+1 $}
%line 9
    \EndIf
%line 10
\EndWhile

%line 11
\State \Return $\HA$

\end{algorithmic}
\end{algorithm}
```


```{=latex}
\begin{algorithm}[H]
\caption*{\textbf{Algorithm 31} $\RejBoundedPoly(\rho: \SixtySixBytes) \rightarrow a: R_{[-\eta, \eta]} $}
\vspace*{2pt}
{Samples an element $a \in R$ with coefficients in $[-\eta, \eta]$ computed via rejection sampling from $\rho$.}

\begin{algorithmic}[1]
%line 1
\State {$ j \is 0 $}

%line 2
\State {$ \Ctx \is \HInit() $} \Comment{H $\doteq$ SHAKE256}

%line 3
\State {$ \Ctx \is \HAbsorb{\Ctx}{\rho} $}

%line 4
\While {$ j < 256 $}
    \State {$ z: \OneByte \is \HSqueeze{\Ctx}{1} $}
    \State {$ \zZero \is \CoeffFromHalfByte(z \Mod{16}) $}
    \State {$ \zOne \is \CoeffFromHalfByte(z\,/\,16) $}
    \If {$ \zZero \neq \bot $}
        \State {$ a_j \is \zZero $}
        \State {$ j \is j+1 $}
    \EndIf
    \If {$ \zOne \neq \bot \ \mathbf{and} \ j < 256 $}
        \State {$ a_j \is \zOne $}
        \State {$ j \is j+1 $}
    \EndIf
\EndWhile

\State \Return $a$

\end{algorithmic}
\end{algorithm}
```



## The Key Generation Algorithm
```{=latex}
\begin{algorithm}[H]
\caption*{\textbf{Algorithm 6} \text{ML-DSA.KeyGen\_internal}($\xi$)}
\begin{algorithmic}[1]

% line 1
\State {$(\rho, \rho', K) \is \NH(\xi \, \Vert \, \IntegerToBytes(k, 1) \, \Vert \, \IntegerToBytes(\ell, 1), 128)$}

% line 2
\State {\Comment{expand seed}}

% line 3
\State {$ \HatA \is \ExpandA(\rho) $}

% line 4
\State {$(\SOne, \STwo) \is \ExpandS(\RhoDash) $}

% line 5
\State {$\Bold{t} \is \NTTInv{\HatA \circ \NTT{\SOne}} + \STwo $}

% line 6
\State {$ (\TOne, \TZero) \is \PowerTwoRound(\Bold{t}) $} \Comment{compress $\Bold{t}$ }

% line 7
\State \Comment{Power2Round is applied componentwise}

% line 8
\State {$ pk \is \pkEncode(\rho, \TOne) $}

% line 9
\State {$ tr \is \NH(pk, 64) $}

% line 10
\State {$ sk \is \skEncode( \rho, K, tr, \SOne, \STwo, \TZero) $}

% line 11
\State {\Return $(pk, sk)$ }

\end{algorithmic}
\end{algorithm}
```

## The Annotated Key Generation Algorithm
```{=latex}
\begin{algorithm}[H]
\caption*{\textbf{Annotated Algorithm 6} \text{ML-DSA.KeyGen\_internal}($ \xi: \BytesThirtyTwo $)}
\begin{algorithmic}[1]

% line 1
\State {$(\rho: \BytesThirtyTwo, \rho': \BytesSixtyFour, K: \BytesThirtyTwo) \is \NH(\xi \, \Vert \, \IntegerToBytes(k, 1) \, \Vert \, \IntegerToBytes(\ell, 1), 128)$}

% line 2
\State {\Comment{expand seed}}

% line 3
\State {$\HatA: \RingTqKxL \is \ExpandA(\rho)$}

% line 4
\State {$(\SOne: \RingRmL, \STwo: \RingRmK) \is \ExpandS(\RhoDash)$ }
    \Statex \TComment{$m \in [ -2, 2 ]$ \textit{if} ML-DSA-44 or ML-DSA-87 } \Comment{$\eta = 2$}
    \Statex \TComment{$m \in [ -4, 4 ]$ \textit{otherwise}} \Comment{$\eta = 4$ }

% line 5
\State {$\Bold{t}: \RingRqK \is \NTTInv{\HatA \circ \NTT{\SOne}} + \STwo $}

% line 6
\State {$ (\TOne: \RingRqOneK, \TZero: \RingRqZeroK) \is \PowerTwoRound(\Bold{t}) $}
    \Comment{compress $\Bold{t}$ }
    \Statex \TComment{$ \TOne \in [0, 1023] $} \Comment{10-bit value }
    \Statex \TComment{$ \TZero \in [ -4095, 4096 ] $} \Comment{$\ModPlusMinus{2^d}, \, d = 13, \, ([ -2^{12}+1, 2^{12} ])$}

% line 7
\State \Comment{Power2Round is applied componentwise }

% line 8
\State {$ pk \is \pkEncode(\rho, \TOne)$}

% line 9
\State {$ tr: \BytesSixtyFour \is \NH(pk, 64) $}

% line 10
\State {$ sk \is \skEncode(\rho, K, tr, \SOne, \STwo, \TZero)$}

% line 11
\State {$\textbf{return} \ (pk, sk)$}

\end{algorithmic}
\end{algorithm}
```


## The Annotated Signature Algorithm

The ML-DSA Signing (Internal) algorithm named \textbf{ML-DSA.Sign\_internal} in FIPS 204 standard is reproduced below for reference. The line numbers and the pseudocode matches exactly with the original version.

```{=latex}
\begin{algorithm}[H]
\caption*{\textbf{Algorithm 7} \text{ML-DSA.Sign\_internal}($ sk, M', rnd $)}
\begin{algorithmic}[1]
\State {$ (\rho, K, tr, \SOne, \STwo, \TZero) \is \skDecode(sk) $}
\State {$ \HatSOne  \is \NTT{\SOne} $}
\State {$ \HatSTwo  \is \NTT{\STwo} $}
\State {$ \HatTZero \is \NTT{\TZero} $}
\State {$ \HatA     \is \ExpandA(\rho) $}
\State $\mu \is \NH(\text{BytesToBits}(tr \,\Vert\, M', 64))$
\State $\rho'' \is \NH(K \,\Vert\, rnd \,\Vert\, \mu, 64)$
\State $\kappa \is 0$
\State $(\z, \h) \is \bot$

\While{(\z, \h) = $\bot$}
    \State {$ \y \in \RingRqL \is \ExpandMask(\rho'', \kappa) $}

    \State {$ \w \is \NTTInv{\HatA \circ \NTT{\y}} $}

    \State {$ \wOne \is \HighBits(\w) $}
        \State \Comment{HighBits is applied componentwise}

    \State {$ \tildeC \is \NH(\mu \,\Vert\, \wOneEncode(\wOne), \lambda/4) $}

    \State {$ \nc \in \Rq \is \SampleInBall(\tildeC) $}

    \State {$ \HatC \is \NTT{c} $}

    \State {$ \CSOne \is \NTTInv{\HatC \circ \HatSOne} $}

    \State {$ \CSTwo \is \NTTInv{\HatC \circ \HatSTwo} $}

    \State {$ \z \is \y +\, \CSOne $}

    \State {$ \rZero \is \LowBits(\w -\, \CSTwo) $}
        \State \Comment{LowBits is applied componentwise}

    \If { $\InfNormOfZ \ge \Diff{\gamma_1}{\beta}$ \, \textbf{or} \, $\InfNormOfRZero \ge \Diff{\gamma_2}{\beta}$}
        $(\z, \h) \is \bot$
    \Else
        \State {$ \CTZero \is \NTTInv{\HatC \circ \HatTZero} $}
        \State {$ \h \is \MakeHint(-\CTZero, \w - \CSTwo + \CTZero) $}
            \State \Comment{MakeHint is applied componentwise}
        \If{ $\InfNormOfCTZero \ge \gamma_2$ \, \textbf{or} \, the number of 1's in $\h$ \ is greater than $\omega$ }
            {$ (\z, \h) \is \bot $}
        \EndIf
    \EndIf
    \State {$ \kappa \is \kappa + \ell $}
 \EndWhile

 \State {$ \sigma \is \sigEncode(\tildeC, \z \ModPlusMinus{q}, \h) $}
 \State \Return $\sigma$

\end{algorithmic}
\end{algorithm}
```

```{=latex}
\begin{algorithm}[H]
\caption*{\textbf{Annotated Algorithm 7} \text{ML-DSA.Sign\_internal}(\textit{sk, M', rnd})}
\begin{algorithmic}[1]

%line 1
\State {$ ( \rho:   \BytesThirtyTwo, \,
            K:      \BytesThirtyTwo, \,
            tr:     \BytesSixtyFour, \,
            \SOne:  \RmL, \,
            \STwo:  \RmK, \,
            \TZero: \RtK) = \skDecode(sk) $}
    \Statex{\TComment\LComment{$m \in [ -2, 2 ]$ \textit{if} ML-DSA-44 \textit{or} ML-DSA-87}}
        \Comment{$\eta = 2$}
    \Statex{\TComment\LComment{$ m \in [ -4, 4 ]$ \textit{otherwise}}}
        \Comment{$\eta = 4$}
    \Statex \TComment\LComment{$ t \in [ -2^{12}+1, 2^{12}-1 ] $}
        \Comment{$d = 13$}

%line 2
\State {$ \HatSOne: \TqL \is \NTT{\SOne} $}

%line 3
\State {$ \HatSTwo: \TqK \is \NTT{\STwo} $}

%line 4
\State {$ \HatTZero: \TqK \is \NTT{\TZero} $}

%line 5
\State {$ \HatA: \RingTqKxL \is \ExpandA(\rho) $}

%line 6
\State {$ \mu: \BytesSixtyFour \is \NH(\BytesToBits(tr \,\Vert\, M', 64)) $}

%line 7
\State $\rho'': \BytesSixtyFour \is \NH(K \,\Vert\, rnd \,\Vert\, \mu, 64)$

%line 8
\State {$ \kappa: \text{u16} \is 0 $}
    \Comment{unsigned 16-bit integer}

%line 9
\State {$ (\z: R_{mod^{\pm} q}^l, \, \h: R_2^k) \is \bot $}

%line 10
\While{ (\z, \h) = $\bot$ } \Comment{Rejection sampling loop}

%line 11
    \State {$ \y: \RingRqL \is \ExpandMask(\rho'', \kappa) $}
        \Comment{Sample a fresh mask $\y$ mixing $\rho''$ and new $\kappa$ }

%line 12
    \State {$ \w: \RqK \is \NTTInv{\HatA \circ \NTT{\y}} $}

%line 13
    \State {$ \wOne: R_{m_1}^k \is \HighBits(\w) $}
        \Statex \TComment\LComment{$ m_1 \in [ 0, 21 ]$ \textit{if} ML-DSA-44}
        \Statex \TComment\LComment{$ m_1 \in [ 0, 15 ]$ \textit{if} ML-DSA-65 \,or\, ML-DSA-87 }

%line 14
        \State \Comment{\text{HighBits is applied componentwise}}

%line 15
    \State {$ \tildeC: \mathbb{B}^{\lambda/4} \is \NH(\mu \Vert \wOneEncode(\wOne), \lambda/4) $}
        \Comment{ $\lambda$ is the collision strength of $\tildeC$ }
        \Statex \TComment\LComment{$\lambda/4$ = 32 (ML-DSA-44) or 48 (ML-DSA-65) or 64 (ML-DSA-87) }

%line 16
    \State {$ \nc: \Rq \is \SampleInBall(\tildeC) $}
        \Comment{polynomial in $R$ with coefficients from \{-1, 0, 1\} }

%line 17
    \State {$\HatC: T_q \is \NTT{c} $}

%line 18
    \State {$ \CSOne: \RingRqL \is \NTTInv{\HatC \circ \HatSOne} $}
        \Comment{$\HatC$ is multiplied with each of $\HatSOne^0, \ldots, \HatSOne^{\ell-1}$ }

%line 19
    \State {$ \CSTwo: \RingRqK \is \NTTInv{\HatC \circ \HatSTwo} $}
        \Comment{$\HatC$ is multiplied with each of $\HatSTwo^0, \ldots, \HatSTwo^{k-1} $}

%line 20
    \State {$ \z: R_{\ModPlusMinus{q}}^{\ell} \is \y + \,\CSOne $}
        \Comment{$\z: R_{\ModPlusMinus{q}}^{\ell} \is \y: \RingRqL + \CSOne: \RingRqL $}

%line 21
    \State {$\rZero : R_q^k \is \LowBits(\w -\, \CSTwo) $}

%line 22
    \State \Comment{\text{LowBits is applied componentwise}}

%line 23
    \If{$\InfNormOfZ \ge \Diff{\gamma_1}{\beta}$ \, \textbf{or} \, $\InfNormOfRZero \ge \Diff{\gamma_2}{\beta}$}
        $(\z, \h) \is \bot$

%line 24
    \Else

%line 25
        \State {$ \CTZero: T_q^k \is \NTTInv{\HatC \circ \HatTZero} $}
        \Comment{$\HatC$ is multiplied with each of $\HatBold{t}_0^0, \ldots, \HatBold{t}_0^{k-1}$ }

%line 26
        \State {$ \h: R_2^k \is \MakeHint(-\CTZero, \w - \CSTwo + \CTZero) $}
%line 27
        \State \Comment{MakeHint is applied componentwise  }
%line 28
        \If{$\Vert \langle\!\langle c\textbf{t}_0\rangle\!\rangle \Vert_\infty \ge \gamma_2 \ \, \textbf{or} \ \, \text{the number of 1's in} \ \textbf{h} \ \text{is greater than} \ \omega$} $(\textbf{z}, \textbf{h}) \is \bot $
%line 29
        \EndIf
%line 30
    \EndIf
%line 31
    \State {$ \kappa \is \kappa + \ell $}
%line 32
 \EndWhile
%line 33
 \State {$ \sigma \is \sigEncode(\tildeC, \z \ModPlusMinus{q}, \h) $}

%line 34
 \State \Return $\sigma$

\end{algorithmic}
\end{algorithm}
```
