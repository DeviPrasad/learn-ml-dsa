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

# Efficient Implementation of Polynomial Arithmetic
Number Theoretic Transforms (NTTs) are integral to the ML DSA specification. The NTT representation enables efficient polynomial addition, subtraction, and multiplication. In this section, we develop all the key ideas behind NTTs in a step-by-step fashion. To illustrate these concepts, we define *Tiny DSA*, a simplified system that follows the ML-DSA structure. This development draws heavily from Prof. Alfred Menezes' excellent notes and videos.

## Polynomial Reduction
Let $a(x) \in R_q$ be a polynomial of degree $\le$ 31. We can reorganize the terms in $a(x)$ into two distinct polynomials. The first polynomial, $a_L(x)$, contains terms with degree $0 \le d \le 15$:
$$
    {\color{NavyBlue}{\mathbf{a_L(x)}}} = a_0 + a_1 x + a_2 x^2 + \ldots + a_{15} x^{15} \in R_q
$$

The second polynomial, $a_r(x)$ has terms with degree $16 \le d \le 31$:
$$
\begin{array}{lll}
    a_r(x) & = & a_{16}x^{16} + a_{17} x^{17} + a_{18} x^{18} + \ldots + a_{31} x^{31} \in R_q \\
           & = & (a_{16} x + a_{17} x^{1} + a_{18} x^{2} + \ldots + a_{31} x^{15}) x^{16} \in R_q \\
           & = & {\color{NavyBlue}{\mathbf{a_R(x)}}} x^{16} \in R_q \ \ \ where \ {\color{NavyBlue}{\mathbf{a_R(x)}}} = (a_{16} x + a_{17} x^{1} + \ldots + a_{31} x^{15}) x^{16} \\
\end{array}
$$

We use this method to repeatedly decompose a polynomial $a(x) \in R_q$ of degree $2n$ into two sub-polynomials each of degree $d$. We stop reducing when we reach factors of degree 1.

$$
    a(x) = a_L(x) + a_R(x)x^{16} \ \ \ \ a_L(x) \in R_q, \ a_R(x) \in R_q, \ a_L(x) \ \text{and} \ a_R(x) \ \text{have  degree} \le \text{15}.
$$

## Tiny DSA - Factoring the Reduction Polynomial $x^{32} + 1$
In this section we derive the irreducible factors of the reduction polynomial $x^{32} + 1$.

The first interesting fact to notice is the following equivalence under the polynomial ring $R_q$
$$
\begin{array}{ll}
    x^{32} + 1 & = x^{32} + \ {\boxed{\color{magenta}\zeta^{0}}} \\
               & = x^{32} \ {\boxed{\color{magenta} - \zeta^{32 + 0}}} \\
               & = x^{32} - \zeta^{32}
\end{array}
$$

In general, we replace a term of the form
$$
    x^{m} + \zeta^n \in R_q  \qquad m \ge 1, n \ge 0
$$
by
$$
    x^{m} - \zeta^{32 + n}
$$
and break this down further till we reach an irreducible term of degree 1, $x \pm z^k \ \, 1 \le k \le 63$. In Tiny DSA, we obtain 32 irreducible terms

$$
\begin{array}{ll}
    x \pm z^k & \ \, k = 2i+1, \ 0 \le i \le 15. \\
              & \text{In other words,} \ k \in \{ 1, 3, 5, 7, 9, \ldots, 23, 25, 27, 29, 31 \}.
\end{array}
$$

In ML DSA, there are 256 irreducible factors of degree 1:
$$
\begin{array}{ll}
    x \pm z^k & \ k = 2i+1, \ 0 \le i \le 127. \\
              & \text{In other words,} \ k \in \{ 1, 3, 5, 7, 9, \ldots, 23, 25, 27, 29, 31 \}.
\end{array}
$$

### Reducing Polynomial $a(x)$ by the factors of $x^{32} - \zeta^{32}$

This is the first (level 1) reduction of the given polynomial $a(x)$.

```{=latex}
\noindent
\begin{minipage}[t]{0.49\textwidth}
  \begin{tikzpicture}[
        scale=0.9,
        transform shape,
        level distance=2cm,
        sibling distance=4cm,
        every label/.style={text=blue, font=\small}
    ]
    % root of the expression tree
    \node (root) {$x^{32} - \zeta^{32}$}
        % left
        child {
            % node[label={[]below:{1}}] (left0) {$x^{16} - \zeta^{16}$}
            node (left0) {$x^{16} - \zeta^{16}$}
            edge from parent node[left, xshift=-2mm, red] {0}
        }
        % right
        child {
            % node (right) {$x^{16} - \zeta^{48}$}
            node (right1) {$x^{16} + \zeta^{16}$}
            edge from parent node[right, xshift=2mm, red] {1}
        };
  \end{tikzpicture}
  \captionof{figure}{Step 1 of Polynomial Factoring}
  \label{fig:ntt-poly-decompose-l1}
\end{minipage}
\hfill
\begin{minipage}[t]{0.49\textwidth}
  \centering
  \begin{tikzpicture}[
        scale=0.9,
        transform shape,
        level distance=2cm,
        sibling distance=5cm,
        every label/.style={text=blue, font=\small}
    ]
    \node (root) {$a(x)$}
        child {
            node [label={above:{$a_0(x)$}}] (left0) {$a(x) \ mod \ (x^{16} - \zeta^{16})$}
            edge from parent node[left, xshift=-2mm, red] {0}
        }
        child {
            node [label={above:{$a_1(x)$}}] (right1) {$a(x) \ mod \ (x^{16} + \zeta^{16})$}
            edge from parent node[right, xshift=2mm, red] {1}
        };
  \end{tikzpicture}
  \captionof{figure}{Step 1 of NTT Evaluation}
  \label{fig:ntt-eval-l1}
\end{minipage}
```


```{=html}
<div class="centered-picture">
<script type="text/tikz">
\begin{tikzpicture}[
        scale=1.5,
        transform shape,
        level distance=2cm,
        sibling distance=4cm,
        every label/.style={text=blue, font=\large},
        every node/.style={font=\small}
    ]
    \node (root) {$a(x)$}
        child {
            node [label={[xshift=-3mm, yshift=-1.5mm] above:{$a_0(x)$}}] (left0) {$a(x) \ mod \ (x^{16} - \zeta^{16})$}
            edge from parent node[left, xshift=-2mm, red] {0}
        }
        child {
            node [label={[xshift=3mm, yshift=-1.5mm] above:{$a_1(x)$}}] (right1) {$a(x) \ mod \ (x^{16} + \zeta^{16})$}
            edge from parent node[right, xshift=2mm, red] {1}
        };
\end{tikzpicture}
</script>
</div>

<div class="centered-picture" style="margin-top: 1.5cm;">
<script type="text/tikz">
\begin{tikzpicture}[
        scale=1.5,
        transform shape,
        level distance=2cm,
        sibling distance=4cm,
        every label/.style={text=blue, font=\large},
        every node/.style={font=\large}
    ]
    \node (root) {$x^{32} - \zeta^{32}$}
        child {
            % node[label={[yshift=0.5mm]below:{1}}] (left) {$x^{16} - \zeta^{16}$}
            node (left0) {$x^{16} - \zeta^{16}$}
            edge from parent node[left, xshift=-2mm, red] {0}
        }
        child {
            % node (right1) {$x^{16} - \zeta^{48}$}
            % node [label={right:{= $x^{16} - \zeta^{48}$}}] (right1) {$x^{16} + \zeta^{16}$}
            node (right1) {$x^{16} + \zeta^{16}$}
            edge from parent node[right, xshift=2mm, red] {1}
        };
\end{tikzpicture}
</script>
</div>
<div style="clear: both;"></div>
<div style="float: none;"></div>
```


$$
\begin{array}{lll}
    a_0(x) & = & a(x) \ \text{mod} \ (x^{16} - \zeta^{16}) \\
           & = & (a_L(x) + a_R(x) \, x^{16}) \ \text{mod} \ (x^{16} - \zeta^{16}) \\
           & = & a_l(x) + a_R(x) \, \zeta^{16} \ \ \ \ \ \textit{(polynomial remainder theorem)}\\
           \\
    a_1(x) & = & a(x) \ \text{mod} \ (x^{16} + \zeta^{16}) \\
           & = & (a_L(x) + a_R(x) \, x^{16}) \ \text{mod} \ (x^{16} + \zeta^{16}) \\
           & = & (a_L(x) + a_R(x) \, x^{16}) \ \text{mod} \ (x^{16} - (-\zeta^{16})) \\
           & = & a_l(x) + a_R(x) \cdot (-\zeta^{16}) \ \ \ \ \ \textit{(polynomial remainder theorem)} \\
           & = & a_l(x) - a_R(x) \, \zeta^{16}\\
\end{array}
$$



### Level 2 - Factors of $x^{16} - \zeta^{16}$ and $x^{16} - \zeta^{48}$
We will now factor the smaller-degree polynomials obtained in the previous step:

```{=latex}
\begin{center}
\begin{tikzpicture}[
        scale=0.9,
        transform shape,
        level distance=2cm,
        level 1/.style={sibling distance=8cm},
        level 2/.style={sibling distance=4cm},
        every label/.style={text=blue, font=\small}
    ]
    % root of the expression tree
    \node (root) {$x^{32} - \zeta^{32}$}
        % left - 0 - x^16-1
        child {
            % node[label={[]below:{1}}] (left) {$x^{16} - \zeta^{16}$}
            node (left0) {$x^{16} - \zeta^{16}$}
            % left - 00
            child {
                % node[label={[]below:{2}}] {$x^{8} - \zeta^{8}$}
                node (left00) {$x^{8} - \zeta^{8}$}
                edge from parent node[left, xshift=-2mm, red] {00}
            }
            % right - 01
            child {
                node (right01) {$x^{8} + \zeta^{8}$}
                edge from parent node[right, xshift=2mm, red] {01}
            }
            edge from parent node[left, xshift=-3mm, red] {0}
        }
        % right -1
        child {
            % node (right1) {$x^{16} - \zeta^{48}$}
            node (right1) {$x^{16} + \zeta^{16}$}
            % left - 10
            child {
                % node[label={[]below:{3}}] (left10) {$x^{8} - \zeta^{24}$}
                node (left10) {$x^{8} - \zeta^{24}$}
                edge from parent node[left, xshift=-2mm, red] {10}
            }
            % right - 11
            child {
                % node (right11) {$x^{8} - \zeta^{56}$}
                node (right11) {$x^{8} + \zeta^{24}$}
                edge from parent node[right, xshift=2mm, red] {11}
            }
            edge from parent node[right, xshift=3mm, red] {1}
        };
\end{tikzpicture}
\end{center}
\captionof{figure}{Step 2 of Polynomial Factoring}
\label{fig:ntt-poly-decompose-l2}
```

```{=latex}
\begin{center}
\begin{tikzpicture}[
        scale=0.9,
        transform shape,
        level distance=2cm,
        level 1/.style={sibling distance=10cm},
        level 2/.style={sibling distance=5cm},
        every label/.style={text=blue, font=\small}
    ]
    \node (root) {$a(x)$}
        child {
            node (left0) {$a_0(x)$}
            child {
                node [label={[xshift=-3mm, yshift=-1.5mm] above:{$a_{00}(x)$}}] (left00) {$a_0(x) \ mod \ (x^{8} - \zeta^{8})$}
                edge from parent node[left, xshift=-2mm, red] {00}
            }
            child {
                node [label={[xshift=3mm, yshift=-1.5mm] above:{$a_{01}(x)$}}] (left01) {$a_0(x) \ mod \ (x^{8} + \zeta^{8})$}
                edge from parent node[right, xshift=2mm, red] {01}
            }
            edge from parent node[left, xshift=-2mm, red] {0}
        }
        child {
            node (right1) {$a_1(x)$}
            child {
                node [label={[xshift=-3mm, yshift=-1.5mm] above:{$a_{10}(x)$}}] (right10) {$a_1(x) \ mod \ (x^{8} - \zeta^{24})$}
                edge from parent node[left, xshift=-2mm, red] {10}
            }
            child {
                node [label={[xshift=3mm, yshift=-1.5mm] above:{$a_{11}(x)$}}] (right11) {$a_1(x) \ mod \ (x^{8} + \zeta^{24})$}
                edge from parent node[right, xshift=2mm, red] {11}
            }
            edge from parent node[right, xshift=2mm, red] {1}
        };
\end{tikzpicture}
\end{center}
\captionof{figure}{Step 2 of NTT Evaluation}
\label{fig:ntt-eval-l2}
```


```{=html}
<div class="centered-picture">
<script type="text/tikz">
\begin{tikzpicture}[
        scale=1.5,
        transform shape,
        level distance=2cm,
        level 1/.style={sibling distance=7.5cm},
        level 2/.style={sibling distance=4cm},
        every label/.style={text=blue, font=\large},
        every node/.style={font=\small}
    ]
    \node (root) {$a(x)$}
        child {
            node (left0) {$a_0(x)$}
            child {
                node [label={[xshift=-3mm, yshift=-1.5mm] above:{$a_{00}(x)$}}] (left00) {$a_0(x) \ mod \ (x^{8} - \zeta^{8})$}
                edge from parent node[left, xshift=-2mm, red] {00}
            }
            child {
                node [label={[xshift=3mm, yshift=-1.5mm] above:{$a_{01}(x)$}}] (left01) {$a_0(x) \ mod \ (x^{8} + \zeta^{8})$}
                edge from parent node[right, xshift=2mm, red] {01}
            }
            edge from parent node[left, xshift=-2mm, red] {0}
        }
        child {
            node (right1) {$a_1(x)$}
            child {
                node [label={[xshift=-3mm, yshift=-1.5mm] above:{$a_{10}(x)$}}] (right10) {$a_1(x) \ mod \ (x^{8} - \zeta^{24})$}
                edge from parent node[left, xshift=-2mm, red] {10}
            }
            child {
                node [label={[xshift=3mm, yshift=-1.5mm] above:{$a_{11}(x)$}}] (right11) {$a_1(x) \ mod \ (x^{8} + \zeta^{24})$}
                edge from parent node[right, xshift=2mm, red] {11}
            }
            edge from parent node[right, xshift=2mm, red] {1}
        };
\end{tikzpicture}
</script>
</div>

<div style="clear: both;"></div>
<div style="float: none; margin-top: 2cm;"></div>

<div class="centered-picture">
<script type="text/tikz">
\begin{tikzpicture}[
        scale=1.5,
        transform shape,
        level distance=2cm,
        level 1/.style={sibling distance=6cm},
        level 2/.style={sibling distance=3.6cm},
        every label/.style={text=blue, font=\large},
        every node/.style={font=\large}
    ]
    % root of the expression tree
    \node (root) {$x^{32} - \zeta^{32}$}
        % left - 0 - x^16-1
        child {
            % node[label={[]below:{1}}] (left) {$x^{16} - \zeta^{16}$}
            node (left0) {$x^{16} - \zeta^{16}$}
            % left - 00
            child {
                node (left00) {$x^{8} - \zeta^{8}$}
                edge from parent node[left, xshift=-2mm, red] {00}
            }
            % right - 01
            child {
                node (right01) {$x^{8} + \zeta^{8}$}
                edge from parent node[right, xshift=2mm, red] {01}
            }
            edge from parent node[left, xshift=-3mm, red] {0}
        }
        % right -1
        child {
            node (right1) {$x^{16} + \zeta^{16}$}
            % left - 10
            child {
                node (left10) {$x^{8} - \zeta^{24}$}
                edge from parent node[left, xshift=-2mm, red] {10}
            }
            % right - 11
            child {
                node (right11) {$x^{8} + \zeta^{24}$}
                edge from parent node[right, xshift=2mm, red] {11}
            }
            edge from parent node[right, xshift=3mm, red] {1}
        };
\end{tikzpicture}
</script>
</div>
```


### Level 3 - Factors of $x^{8} \pm \zeta^{8}$ and $x^{8} \pm \zeta^{24}$

We use the equalities $x^{8} + \zeta^{8} = x^{8} - \zeta^{40}$ and $x^{8} + \zeta^{24} = x^{8} - \zeta^{56}$ in evaluating the NTT of polynomial $a(x)$.

```{=latex}
\begin{center}
\begin{tikzpicture}[
        scale=0.9,
        transform shape,
        level distance=2cm,
        level 1/.style={sibling distance=7.5cm},
        level 2/.style={sibling distance=3.75cm},
        level 3/.style={sibling distance=1.8cm},
        every label/.style={text=blue, font=\scriptsize},
    ]
    % root of the expression tree
    \node (root) {$x^{32} - \zeta^{32}$}
        % left - 0 - x^16-1
        child {
            node (left) {$x^{16} - \zeta^{16}$}
            % left - 00
            child {
                node {$x^{8} - \zeta^{8}$}
                child {
                    % left - x^4 - \zeta^4 - 000
                    node {$x^{4} - \zeta^{4}$}
                    edge from parent node[left, xshift=-1mm, red] {000}
                }
                child {
                    % left - x^4 - \zeta^34 - 001
                    node (right001) {$x^{4} - \zeta^{36}$}
                    edge from parent node[right, xshift=1mm, red] {001}
                }
                edge from parent node[left, xshift=-2mm, red] {00}
            }
            % right - 01
            child {
                node (right01) {$x^{8} - \zeta^{40}$}
                child {
                    % left - x^4 - \zeta^4 - 010
                    node {$x^{4} - \zeta^{20}$}
                    edge from parent node[left, xshift=-1mm, red] {010}
                }
                child {
                    % left - x^4 - \zeta^34 - 011
                    node (right011) {$x^{4} - \zeta^{52}$}
                    edge from parent node[right, xshift=1mm, red] {011}
                }
                edge from parent node[right, xshift=2mm, red] {01}
            }
            edge from parent node[left, xshift=-3mm, red] {0}
        }
        % right -1
        child {
            node (right1) {$x^{16} - \zeta^{48}$}
            % left - 10
            child {
                node (left10) {$x^{8} - \zeta^{24}$}
                child {
                    % left - x^4 - \zeta^12 - 100
                    node {$x^{4} - \zeta^{12}$}
                    edge from parent node[left, xshift=-1mm, red] {100}
                }
                child {
                    % left - x^4 - \zeta^34 - 101
                    node (right101) {$x^{4} - \zeta^{44}$}
                    edge from parent node[right, xshift=1mm, red] {101}
                }
                edge from parent node[left, xshift=-2mm, red] {10}
            }
            % right - 11
            child {
                % Renamed from (right) to (right11)
                node (right11) {$x^{8} - \zeta^{56}$}
                child {
                    % left - x^4 - \zeta^28 - 110
                    node {$x^{4} - \zeta^{28}$}
                    edge from parent node[left, xshift=-1mm, red] {110}
                }
                child {
                    % left - x^4 - \zeta^60 - 111
                    node (right111) {$x^{4} - \zeta^{60}$}
                    edge from parent node[right, xshift=1mm, red] {111}
                }
                edge from parent node[right, xshift=2mm, red] {11}
            }
            edge from parent node[right, xshift=3mm, red] {1}
        };
\end{tikzpicture}
\end{center}
\captionof{figure}{Step 3 of Polynomial Factoring}
\label{fig:ntt-poly-decompose-l3}
```

```{=latex}
\begin{center}
\begin{tikzpicture}[
        scale=0.9,
        transform shape,
        level distance=2cm,
        level 1/.style={sibling distance=10cm},
        level 2/.style={sibling distance=5cm},
        level 3/.style={sibling distance=2.5cm},
        every label/.style={text=blue, font=\small}
    ]
    \node (root) {$a(x)$}
        child {
            node (left0) {$a_0(x)$}
            child {
                node (left00) {$a_{00}(x)$}
                child {
                    node (left000) {$a_{000}(x)$}
                    edge from parent node[left, xshift=-2mm, red] {000}
                }
                child {
                    node (right001) {$a_{001}(x)$}
                    edge from parent node[right, xshift=2mm, red] {001}
                }
                edge from parent node[left, xshift=-2mm, red] {00}
            }
            child {
                node (left01) {$a_{01}(x)$}
                child {
                    node (left000) {$a_{010}(x)$}
                    edge from parent node[left, xshift=-2mm, red] {010}
                }
                child {
                    node (right001) {$a_{011}(x)$}
                    edge from parent node[right, xshift=2mm, red] {011}
                }
                edge from parent node[right, xshift=2mm, red] {01}
            }
            edge from parent node[left, xshift=-2mm, red] {0}
        }
        child {
            node (right1) {$a_1(x)$}
            child {
                node (right10) {$a_1(x)$}
                child {
                    node (left100) {$a_{100}(x)$}
                    edge from parent node[left, xshift=-2mm, red] {100}
                }
                child {
                    node (right101) {$a_{101}(x)$}
                    edge from parent node[right, xshift=2mm, red] {101}
                }
                edge from parent node[left, xshift=-2mm, red] {10}
            }
            child {
                node (right11) {$a_{11}(x)$}
                child {
                    node (left110) {$a_{110}(x)$}
                    edge from parent node[left, xshift=-2mm, red] {110}
                }
                child {
                    node (right111) {$a_{111}(x)$}
                    edge from parent node[right, xshift=2mm, red] {111}
                }
                edge from parent node[right, xshift=2mm, red] {11}
            }
            edge from parent node[right, xshift=2mm, red] {1}
        };
\end{tikzpicture}
\end{center}
\captionof{figure}{Step 3 of NTT Evaluation}
\label{fig:ntt-eval-l3}
```

```{=html}
<div class="centered-picture">
<script type="text/tikz">
\begin{tikzpicture}[
        scale=1.5,
        transform shape,
        level distance=2cm,
        level 1/.style={sibling distance=7.5cm},
        level 2/.style={sibling distance=3.75cm},
        level 3/.style={sibling distance=1.8cm},
        every label/.style={text=blue, font=\large},
        every node/.style={font=\large}
]
    % root of the expression tree
    \node (root) {$x^{32} - \zeta^{32}$}
        % left - 0 - x^16-1
        child {
            node (left) {$x^{16} - \zeta^{16}$}
            % left - 00
            child {
                node {$x^{8} - \zeta^{8}$}
                child {
                    % left - x^4 - \zeta^4 - 000
                    node {$x^{4} - \zeta^{4}$}
                    edge from parent node[left, xshift=-1mm, red] {000}
                }
                child {
                    % left - x^4 - \zeta^34 - 001
                    node (right001) {$x^{4} - \zeta^{36}$}
                    edge from parent node[right, xshift=1mm, red] {001}
                }
                edge from parent node[left, xshift=-2mm, red] {00}
            }
            % right - 01
            child {
                node (right01) {$x^{8} - \zeta^{40}$}
                child {
                    % left - x^4 - \zeta^4 - 010
                    node {$x^{4} - \zeta^{20}$}
                    edge from parent node[left, xshift=-1mm, red] {010}
                }
                child {
                    % left - x^4 - \zeta^34 - 011
                    node (right011) {$x^{4} - \zeta^{52}$}
                    edge from parent node[right, xshift=1mm, red] {011}
                }
                edge from parent node[right, xshift=2mm, red] {01}
            }
            edge from parent node[left, xshift=-3mm, red] {0}
        }
        % right -1
        child {
            node (right1) {$x^{16} - \zeta^{48}$}
            % left - 10
            child {
                node (left10) {$x^{8} - \zeta^{24}$}
                child {
                    % left - x^4 - \zeta^12 - 100
                    node {$x^{4} - \zeta^{12}$}
                    edge from parent node[left, xshift=-1mm, red] {100}
                }
                child {
                    % left - x^4 - \zeta^34 - 101
                    node (right101) {$x^{4} - \zeta^{44}$}
                    edge from parent node[right, xshift=1mm, red] {101}
                }
                edge from parent node[left, xshift=-2mm, red] {10}
            }
            % right - 11
            child {
                % Renamed from (right) to (right11)
                node (right11) {$x^{8} - \zeta^{56}$}
                child {
                    % left - x^4 - \zeta^28 - 110
                    node {$x^{4} - \zeta^{28}$}
                    edge from parent node[left, xshift=-1mm, red] {110}
                }
                child {
                    % left - x^4 - \zeta^60 - 111
                    node (right111) {$x^{4} - \zeta^{60}$}
                    edge from parent node[right, xshift=1mm, red] {111}
                }
                edge from parent node[right, xshift=2mm, red] {11}
            }
            edge from parent node[right, xshift=3mm, red] {1}
        };
\end{tikzpicture}
</script>
</div>
```


### Level 4

```{=latex}
\newpage
\begin{landscape}
\thispagestyle{empty}
\vspace*{3mm}\text{Number Theoretic Transforms in ML DSA} \hfill \text{Jan-2026}\vspace{-4mm}\\
\noindent\rule{\linewidth}{0.4pt}
\vfill
\centerline{\begin{tikzpicture}[
        scale=0.81,
        transform shape,
        every label/.style={text=blue, font=\small},
        level distance=2.3cm,
        level 1/.style={sibling distance=17.6cm},
        level 2/.style={sibling distance=9cm},
        level 3/.style={sibling distance=4.5cm},
        level 4/.style={sibling distance=2.3cm},
        %level 5/.style={sibling distance=1.5cm},
        level 5/.style={
            sibling distance=1cm,
        },
        level 5 left/.style={
            level distance=2cm,
            edge from parent path={(\tikzparentnode.south) -- (\tikzchildnode.north)}
        },
        level 5 right/.style={
            level distance=4cm,
            edge from parent path={(\tikzparentnode.south) -- (\tikzchildnode.north)}
        },
    ]
    % root of the expression tree
    \node (root) {$x^{32} - \zeta^{32}$}
        % left - 0 - x^16-1
        child {
            node[label={below:{1}}] (left) {$x^{16} - \zeta^{16}$}
            % left - 00
            child {
                node[label={below:{2}}] {$x^{8} - \zeta^{8}$}
                child {
                    % left - x^4 - \zeta^4 - 000
                    node[label={below:{4}}] {$x^{4} - \zeta^{4}$}
                    child {
                        % left - x^2 - \zeta^2 - 0000
                        node[label={above:{8}}] {$x^{2} - \zeta^{2}$}
                        child [level 5 left] {
                            % left: x - \zeta - 00000
                            node[label={[xshift=-2mm, yshift=-2mm]above:{16}}] {$x - \zeta$}
                            edge from parent node[left, xshift=-0.3mm, red] {00000}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{33} - 00001
                            node (right00001) {$x- \zeta^{33}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {00001}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {0000}
                    }
                    child {
                        % left - x^2 - \zeta^34 - 0001
                        node (right0001) {$x^{2} - \zeta^{34}$}
                        child [level 5 left] {
                            % left: x - \zeta^{17} - 00010
                            node[label={[xshift=-2mm, yshift=-2mm]above:{17}}] {$x - \zeta^{17}$}
                            edge from parent node[left, xshift=-0.3mm, red] {00010}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{49} - 00011
                            node (right00011) {$x- \zeta^{33}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {00011}
                        }
                        edge from parent node[right, xshift=0.3mm, yshift=0.3mm, red] {0001}
                    }
                    edge from parent node[left, xshift=-1mm, red] {000}
                }
                child {
                    % left - x^4 - \zeta^36 - 001
                    node (right001) {$x^{4} - \zeta^{36}$}
                    child {
                        % left - x^2 - \zeta^18 - 0010
                        node[label={above:{9}}] {$x^{2} - \zeta^{18}$}
                        child [level 5 left] {
                            % left: x - \zeta^{9} - 00100
                            node[label={[xshift=-2mm, yshift=-2mm]above:{18}}] {$x - \zeta^9$}
                            edge from parent node[left, xshift=-0.3mm, red] {00100}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{41} - 00101
                            node (right00101) {$x- \zeta^{41}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {00101}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {0010}
                    }
                    child {
                        % left - x^2 - \zeta^50 - 0011
                        node (right0011) {$x^{2} - \zeta^{50}$}
                        child [level 5 left] {
                            % left: x - \zeta^{25} - 00110
                            node[label={[xshift=-2mm, yshift=-2mm]above:{19}}] {$x - \zeta^{25}$}
                            edge from parent node[left, xshift=-0.3mm, red] {00110}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{57} - 00111
                            node (right00111) {$x- \zeta^{57}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {00111}
                        }
                        edge from parent node[right, xshift=0.3mm, yshift=0.3mm, red] {0011}
                    }
                    edge from parent node[right, xshift=1mm, red] {001}
                }
                edge from parent node[left, xshift=-2mm, red] {00}
            }
            % right - 01
            child {
                node (right01) {$x^{8} - \zeta^{40}$}
                child {
                    % left - x^4 - \zeta^4 - 010
                    node[label={[]below:{5}}] {$x^{4} - \zeta^{20}$}
                    child {
                        % left - x^2 - \zeta^10 - 0100
                        node[label={above:{10}}] {$x^{2} - \zeta^{10}$}
                        child [level 5 left] {
                            % left: x - \zeta^{5} - 01000
                            node[label={[xshift=-2mm, yshift=-2mm]above:{20}}] {$x - \zeta^{5}$}
                            edge from parent node[left, xshift=-0.3mm, red] {01000}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{37} - 01001
                            node (right01001) {$x- \zeta^{37}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {01001}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {0100}
                    }
                    child {
                        % left - x^2 - \zeta^42 - 0101
                        node (right0101) {$x^{2} - \zeta^{42}$}
                        child [level 5 left] {
                            % left: x - \zeta^{21} - 01010
                            node[label={[xshift=-2mm, yshift=-2mm]above:{21}}] {$x - \zeta^{21}$}
                            edge from parent node[left, xshift=-0.3mm, red] {01010}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{53} - 01011
                            node (right01011) {$x- \zeta^{53}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {01011}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {0101}
                    }
                    edge from parent node[left, xshift=-1mm, red] {010}
                }
                child {
                    % left - x^4 - \zeta^52 - 011
                    node (right011) {$x^{4} - \zeta^{52}$}
                    child {
                        % left - x^2 - \zeta^26 - 0110
                        node[label={above:{11}}] {$x^{2} - \zeta^{26}$}
                        child [level 5 left] {
                            % left: x - \zeta^{13} - 01100
                            node[label={[xshift=-2mm, yshift=-2mm]above:{22}}] {$x - \zeta^{13}$}
                            edge from parent node[left, xshift=-0.3mm, red] {01100}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{45} - 01101
                            node (right01101) {$x-\zeta^{45}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {01101}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {0110}
                    }
                    child {
                        % left - x^2 - \zeta^58 - 0111
                        node (right0111) {$x^{2} - \zeta^{58}$}
                        child [level 5 left] {
                            % left: x - \zeta^{29} - 01110
                            node[label={[xshift=-2mm, yshift=-2mm]above:{23}}] {$x-\zeta^{29}$}
                            edge from parent node[left, xshift=-0.3mm, red] {01110}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{61} - 01111
                            node (right01111) {$x-\zeta^{61}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {01111}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {0111}
                    }
                    edge from parent node[right, xshift=1mm, red] {011}
                }
                edge from parent node[right, xshift=2mm, red] {01}
            }
            edge from parent node[left, yshift=3mm, red] {0}
        }
        % right -1
        child {
            node (right1) {$x^{16} - \zeta^{48}$}
            % left - 10
            child {
                node[label={[]below:{3}}] (left10) {$x^{8} - \zeta^{24}$}
                child {
                    % left - x^4 - \zeta^12 - 100
                    node[label={[]below:{6}}] {$x^{4} - \zeta^{12}$}
                    child {
                        % left - x^2 - \zeta^6 - 1000
                        node[label={above:{12}}] {$x^{2} - \zeta^{6}$}
                        child [level 5 left] {
                            % left: x - \zeta^{3} - 10000
                            node[label={[xshift=-2mm, yshift=-2mm]above:{24}}] {$x-\zeta^{3}$}
                            edge from parent node[left, xshift=-0.3mm, red] {10000}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{35} - 10001
                            node (right10001) {$x-\zeta^{35}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {10001}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {1000}
                    }
                    child {
                        % left - x^2 - \zeta^38 - 1001
                        node (right1001) {$x^{2} - \zeta^{38}$}
                        child [level 5 left] {
                            % left: x - \zeta^{19} - 10010
                            node[label={[xshift=-2mm, yshift=-2mm]above:{25}}] {$x-\zeta^{19}$}
                            edge from parent node[left, xshift=-0.3mm, red] {10010}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{51} - 10001
                            node (right10011) {$x-\zeta^{51}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {10011}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {1001}
                    }
                    edge from parent node[left, xshift=-1mm, red] {100}
                }
                child {
                    % left - x^4 - \zeta^34 - 101
                    node (right101) {$x^{4} - \zeta^{44}$}
                    child {
                        % left - x^2 - \zeta^22 - 1000
                        node[label={above:{13}}] {$x^{2} - \zeta^{22}$}
                        child [level 5 left] {
                            % left: x - \zeta^{11} - 10000
                            node[label={[xshift=-2mm, yshift=-2mm]above:{26}}] {$x-\zeta^{11}$}
                            edge from parent node[left, xshift=-0.3mm, red] {10000}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{43} - 10001
                            node (right10001) {$x-\zeta^{43}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {10001}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {1010}
                    }
                    child {
                        % left - x^2 - \zeta^54 - 1001
                        node (right1011) {$x^{2} - \zeta^{54}$}
                        child [level 5 left] {
                            % left: x - \zeta^{27} - 10010
                            node[label={[xshift=-2mm, yshift=-2mm]above:{27}}] {$x-\zeta^{27}$}
                            edge from parent node[left, xshift=-0.3mm, red] {10010}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{59} - 10011
                            node (right10011) {$x-\zeta^{59}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {10011}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {1011}
                    }
                    edge from parent node[right, xshift=1mm, red] {101}
                }
                edge from parent node[left, xshift=-2mm, red] {10}
            }
            % right - 11
            child {
                % Renamed from (right) to (right11)
                node (right11) {$x^{8} - \zeta^{56}$}
                child {
                    % left - x^4 - \zeta^28 - 110
                    node[label={[]below:{7}}] {$x^{4} - \zeta^{28}$}
                    child {
                        % left - x^2 - \zeta^14 - 1100
                        node[label={above:{14}}] {$x^{2} - \zeta^{14}$}
                        child [level 5 left] {
                            % left: x - \zeta^{7} - 11000
                            node[label={[xshift=-2mm, yshift=-2mm]above:{28}}] {$x-\zeta^{7}$}
                            edge from parent node[left, xshift=-0.3mm, red] {11000}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{39} - 11001
                            node (right11001) {$x-\zeta^{39}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {11001}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {1100}
                    }
                    child {
                        % left - x^2 - \zeta^46 - 1101
                        node (right1101) {$x^{2} - \zeta^{46}$}
                        child [level 5 left] {
                            % left: x - \zeta^{23} - 11010
                            node[label={[xshift=-2mm, yshift=-2mm]above:{29}}] {$x-\zeta^{23}$}
                            edge from parent node[left, xshift=-0.3mm, red] {11010}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{55} - 11011
                            node (right11011) {$x-\zeta^{55}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {11011}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {1101}
                    }
                    edge from parent node[left, xshift=-1mm, red] {110}
                }
                child {
                    % left - x^4 - \zeta^60 - 111
                    node (right111) {$x^{4} - \zeta^{60}$}
                    child {
                        % left - x^2 - \zeta^30 - 1110
                        node[label={above:{15}}] {$x^{2} - \zeta^{30}$}
                        child [level 5 left] {
                            % left: x - \zeta^{15} - 11100
                            node[label={[xshift=-2mm, yshift=-2mm]above:{30}}] {$x-\zeta^{15}$}
                            edge from parent node[left, xshift=-0.3mm, red] {11100}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{47} - 11101
                            node (right11101) {$x-\zeta^{47}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {11101}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {11101}
                    }
                    child {
                        % left - x^2 - \zeta^62 - 1111
                        node (right1111) {$x^{2} - \zeta^{62}$}
                        child [level 5 left] {
                            % left: x - \zeta^{31} - 11110
                            node[label={[xshift=-2mm, yshift=-2mm]above:{31}}] {$x-\zeta^{31}$}
                            edge from parent node[left, xshift=-0.3mm, red] {11110}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{63} - 11111
                            node (right11111) {$x-\zeta^{63}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {11111}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {1111}
                    }
                    edge from parent node[right, xshift=1mm, red] {111}
                }
                edge from parent node[right, xshift=2mm, red] {11}
            }
            edge from parent node[right, yshift=3mm, red] {1}
        };
\end{tikzpicture}}
\vfill
\noindent\rule{\linewidth}{0.4pt}\\
\raggedleft\thepage\vspace{1cm}
\restoregeometry
\end{landscape}
```


```{=html}
<script type="text/tikz">
    \begin{tikzpicture}[
        scale=1.25,
        transform shape,
        every label/.style={text=blue, font=\small},
        level distance=2.3cm,
        level 1/.style={sibling distance=17.6cm},
        level 2/.style={sibling distance=9cm},
        level 3/.style={sibling distance=4.5cm},
        level 4/.style={sibling distance=2.3cm},
        %level 5/.style={sibling distance=1.5cm},
        level 5/.style={
            sibling distance=1cm,
        },
        level 5 left/.style={
            level distance=2cm,
            edge from parent path={(\tikzparentnode.south) -- (\tikzchildnode.north)}
        },
        level 5 right/.style={
            level distance=4cm,
            edge from parent path={(\tikzparentnode.south) -- (\tikzchildnode.north)}
        },
    ]
    % root of the expression tree
    \node (root) {$x^{32} - \zeta^{32}$}
        % left - 0 - x^16-1
        child {
            node[label={below:{1}}] (left) {$x^{16} - \zeta^{16}$}
            % left - 00
            child {
                node[label={below:{2}}] {$x^{8} - \zeta^{8}$}
                child {
                    % left - x^4 - \zeta^4 - 000
                    node[label={below:{4}}] {$x^{4} - \zeta^{4}$}
                    child {
                        % left - x^2 - \zeta^2 - 0000
                        node[label={above:{8}}] {$x^{2} - \zeta^{2}$}
                        child [level 5 left] {
                            % left: x - \zeta - 00000
                            node[label={[xshift=-2mm, yshift=-2mm]above:{16}}] {$x - \zeta$}
                            edge from parent node[left, xshift=-0.3mm, red] {00000}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{33} - 00001
                            node (right00001) {$x- \zeta^{33}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {00001}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {0000}
                    }
                    child {
                        % left - x^2 - \zeta^34 - 0001
                        node (right0001) {$x^{2} - \zeta^{34}$}
                        child [level 5 left] {
                            % left: x - \zeta^{17} - 00010
                            node[label={[xshift=-2mm, yshift=-2mm]above:{17}}] {$x - \zeta^{17}$}
                            edge from parent node[left, xshift=-0.3mm, red] {00010}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{49} - 00011
                            node (right00011) {$x- \zeta^{33}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {00011}
                        }
                        edge from parent node[right, xshift=0.3mm, yshift=0.3mm, red] {0001}
                    }
                    edge from parent node[left, xshift=-1mm, red] {000}
                }
                child {
                    % left - x^4 - \zeta^36 - 001
                    node (right001) {$x^{4} - \zeta^{36}$}
                    child {
                        % left - x^2 - \zeta^18 - 0010
                        node[label={above:{9}}] {$x^{2} - \zeta^{18}$}
                        child [level 5 left] {
                            % left: x - \zeta^{9} - 00100
                            node[label={[xshift=-2mm, yshift=-2mm]above:{18}}] {$x - \zeta^9$}
                            edge from parent node[left, xshift=-0.3mm, red] {00100}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{41} - 00101
                            node (right00101) {$x- \zeta^{41}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {00101}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {0010}
                    }
                    child {
                        % left - x^2 - \zeta^50 - 0011
                        node (right0011) {$x^{2} - \zeta^{50}$}
                        child [level 5 left] {
                            % left: x - \zeta^{25} - 00110
                            node[label={[xshift=-2mm, yshift=-2mm]above:{19}}] {$x - \zeta^{25}$}
                            edge from parent node[left, xshift=-0.3mm, red] {00110}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{57} - 00111
                            node (right00111) {$x- \zeta^{57}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {00111}
                        }
                        edge from parent node[right, xshift=0.3mm, yshift=0.3mm, red] {0011}
                    }
                    edge from parent node[right, xshift=1mm, red] {001}
                }
                edge from parent node[left, xshift=-2mm, red] {00}
            }
            % right - 01
            child {
                node (right01) {$x^{8} - \zeta^{40}$}
                child {
                    % left - x^4 - \zeta^4 - 010
                    node[label={[]below:{5}}] {$x^{4} - \zeta^{20}$}
                    child {
                        % left - x^2 - \zeta^10 - 0100
                        node[label={above:{10}}] {$x^{2} - \zeta^{10}$}
                        child [level 5 left] {
                            % left: x - \zeta^{5} - 01000
                            node[label={[xshift=-2mm, yshift=-2mm]above:{20}}] {$x - \zeta^{5}$}
                            edge from parent node[left, xshift=-0.3mm, red] {01000}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{37} - 01001
                            node (right01001) {$x- \zeta^{37}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {01001}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {0100}
                    }
                    child {
                        % left - x^2 - \zeta^42 - 0101
                        node (right0101) {$x^{2} - \zeta^{42}$}
                        child [level 5 left] {
                            % left: x - \zeta^{21} - 01010
                            node[label={[xshift=-2mm, yshift=-2mm]above:{21}}] {$x - \zeta^{21}$}
                            edge from parent node[left, xshift=-0.3mm, red] {01010}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{53} - 01011
                            node (right01011) {$x- \zeta^{53}$}
                            edge from parent node[right, xshift=2mm, yshift=-1cm, red] {01011}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {0101}
                    }
                    edge from parent node[left, xshift=-1mm, red] {010}
                }
                child {
                    % left - x^4 - \zeta^52 - 011
                    node (right011) {$x^{4} - \zeta^{52}$}
                    child {
                        % left - x^2 - \zeta^26 - 0110
                        node[label={above:{11}}] {$x^{2} - \zeta^{26}$}
                        child [level 5 left] {
                            % left: x - \zeta^{13} - 01100
                            node[label={[xshift=-2mm, yshift=-2mm]above:{22}}] {$x - \zeta^{13}$}
                            edge from parent node[left, xshift=-0.3mm, red] {01100}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{45} - 01101
                            node (right01101) {$x-\zeta^{45}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {01101}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {0110}
                    }
                    child {
                        % left - x^2 - \zeta^58 - 0111
                        node (right0111) {$x^{2} - \zeta^{58}$}
                        child [level 5 left] {
                            % left: x - \zeta^{29} - 01110
                            node[label={[xshift=-2mm, yshift=-2mm]above:{23}}] {$x-\zeta^{29}$}
                            edge from parent node[left, xshift=-0.3mm, red] {01110}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{61} - 01111
                            node (right01111) {$x-\zeta^{61}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {01111}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {0111}
                    }
                    edge from parent node[right, xshift=1mm, red] {011}
                }
                edge from parent node[right, xshift=2mm, red] {01}
            }
            edge from parent node[left, yshift=3mm, red] {0}
        }
        % right -1
        child {
            node (right1) {$x^{16} - \zeta^{48}$}
            % left - 10
            child {
                node[label={[]below:{3}}] (left10) {$x^{8} - \zeta^{24}$}
                child {
                    % left - x^4 - \zeta^12 - 100
                    node[label={[]below:{6}}] {$x^{4} - \zeta^{12}$}
                    child {
                        % left - x^2 - \zeta^6 - 1000
                        node[label={above:{12}}] {$x^{2} - \zeta^{6}$}
                        child [level 5 left] {
                            % left: x - \zeta^{3} - 10000
                            node[label={[xshift=-2mm, yshift=-2mm]above:{24}}] {$x-\zeta^{3}$}
                            edge from parent node[left, xshift=-0.3mm, red] {10000}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{35} - 10001
                            node (right10001) {$x-\zeta^{35}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {10001}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {1000}
                    }
                    child {
                        % left - x^2 - \zeta^38 - 1001
                        node (right1001) {$x^{2} - \zeta^{38}$}
                        child [level 5 left] {
                            % left: x - \zeta^{19} - 10010
                            node[label={[xshift=-2mm, yshift=-2mm]above:{25}}] {$x-\zeta^{19}$}
                            edge from parent node[left, xshift=-0.3mm, red] {10010}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{51} - 10001
                            node (right10011) {$x-\zeta^{51}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {10011}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {1001}
                    }
                    edge from parent node[left, xshift=-1mm, red] {100}
                }
                child {
                    % left - x^4 - \zeta^34 - 101
                    node (right101) {$x^{4} - \zeta^{44}$}
                    child {
                        % left - x^2 - \zeta^22 - 1000
                        node[label={above:{13}}] {$x^{2} - \zeta^{22}$}
                        child [level 5 left] {
                            % left: x - \zeta^{11} - 10000
                            node[label={[xshift=-2mm, yshift=-2mm]above:{26}}] {$x-\zeta^{11}$}
                            edge from parent node[left, xshift=-0.3mm, red] {10000}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{43} - 10001
                            node (right10001) {$x-\zeta^{43}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {10001}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {1010}
                    }
                    child {
                        % left - x^2 - \zeta^54 - 1001
                        node (right1011) {$x^{2} - \zeta^{54}$}
                        child [level 5 left] {
                            % left: x - \zeta^{27} - 10010
                            node[label={[xshift=-2mm, yshift=-2mm]above:{27}}] {$x-\zeta^{27}$}
                            edge from parent node[left, xshift=-0.3mm, red] {10010}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{59} - 10011
                            node (right10011) {$x-\zeta^{59}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {10011}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {1011}
                    }
                    edge from parent node[right, xshift=1mm, red] {101}
                }
                edge from parent node[left, xshift=-2mm, red] {10}
            }
            % right - 11
            child {
                % Renamed from (right) to (right11)
                node (right11) {$x^{8} - \zeta^{56}$}
                child {
                    % left - x^4 - \zeta^28 - 110
                    node[label={[]below:{7}}] {$x^{4} - \zeta^{28}$}
                    child {
                        % left - x^2 - \zeta^14 - 1100
                        node[label={above:{14}}] {$x^{2} - \zeta^{14}$}
                        child [level 5 left] {
                            % left: x - \zeta^{7} - 11000
                            node[label={[xshift=-2mm, yshift=-2mm]above:{28}}] {$x-\zeta^{7}$}
                            edge from parent node[left, xshift=-0.3mm, red] {11000}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{39} - 11001
                            node (right11001) {$x-\zeta^{39}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {11001}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {1100}
                    }
                    child {
                        % left - x^2 - \zeta^46 - 1101
                        node (right1101) {$x^{2} - \zeta^{46}$}
                        child [level 5 left] {
                            % left: x - \zeta^{23} - 11010
                            node[label={[xshift=-2mm, yshift=-2mm]above:{29}}] {$x-\zeta^{23}$}
                            edge from parent node[left, xshift=-0.3mm, red] {11010}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{55} - 11011
                            node (right11011) {$x-\zeta^{55}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {11011}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {1101}
                    }
                    edge from parent node[left, xshift=-1mm, red] {110}
                }
                child {
                    % left - x^4 - \zeta^60 - 111
                    node (right111) {$x^{4} - \zeta^{60}$}
                    child {
                        % left - x^2 - \zeta^30 - 1110
                        node[label={above:{15}}] {$x^{2} - \zeta^{30}$}
                        child [level 5 left] {
                            % left: x - \zeta^{15} - 11100
                            node[label={[xshift=-2mm, yshift=-2mm]above:{30}}] {$x-\zeta^{15}$}
                            edge from parent node[left, xshift=-0.3mm, red] {11100}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{47} - 11101
                            node (right11101) {$x-\zeta^{47}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {11101}
                        }
                        edge from parent node[left, xshift=-0.3mm, red] {11101}
                    }
                    child {
                        % left - x^2 - \zeta^62 - 1111
                        node (right1111) {$x^{2} - \zeta^{62}$}
                        child [level 5 left] {
                            % left: x - \zeta^{31} - 11110
                            node[label={[xshift=-2mm, yshift=-2mm]above:{31}}] {$x-\zeta^{31}$}
                            edge from parent node[left, xshift=-0.3mm, red] {11110}
                        }
                        child [level 5 right] {
                            % right: x - \zeta^{63} - 11111
                            node (right11111) {$x-\zeta^{63}$}
                            edge from parent node[right, xshift=1mm, yshift=-1cm, red] {11111}
                        }
                        edge from parent node[right, xshift=0.3mm, red] {1111}
                    }
                    edge from parent node[right, xshift=1mm, red] {111}
                }
                edge from parent node[right, xshift=2mm, red] {11}
            }
            edge from parent node[right, yshift=3mm, red] {1}
        };
    \end{tikzpicture}
</script>
```

# Tiny DSA and Number Theoretic Transforms
Tiny DSA uses the following constants, symbols, and mathematical expressions:

$$
\begin{array}{rll}
&q & \text{The prime number} \ q = 2^{10} - 2^8 + 1 = 769.\\\\
&\mathbb{N} & \text{The set of natural numbers}.\\\\
&\mathbb{Z} & \text{The ring of integers}.\\
&\mathbb{Z}_q & \text{The ring of integers modulo} \ q \ \text{whose set of elements is} \ \{0, 1, \ldots, q-1\}.\\
&\mathbb{Z}^n_q & \text{The set of $n$-tuples over $Z_q$}.\\\\
&T_q & \text{The ring} \ \Pi_{j=0}^{32}\,\mathbb{Z}_q. \\\\
R & Z[X]/(X^{32}+1) &  \text{The ring of single-variable polynomials over} \ \mathbb{Z} \ \text{modulo} \ X^{32}+1. \\
    && \text{The coefficients of polynomials in} \ R \ \text{belong to the ring} \ \mathbb{Z}.\\
    && \text{The highest-degree term is at most} \ x^{31}.\\\\
R_q & Z_q[X]/(X^{32}+1) &  \text{The ring of single-variable polynomials over} \ \mathbb{Z}_q \ \text{modulo} \ X^{32}+1. \\
    && \text{The coefficients of polynomials in} \ R_q \ \text{belong to the ring} \ \mathbb{Z}_q.\\
    && \text{The highest-degree term is at most} \ x^{31}. \\\\
    & (X^{32}+1) & \text{The} \ reduction \ polynomial. \\\\
 \omega, \ z \, & \zeta & \text{The} \ \mathbf{primitive \ 64^{th} \ root} \ \text{of unity in} \ \mathbb{Z}_q. \\
    && \zeta^{64} \equiv 1 \ \text{mod} \ q.\\
    && \zeta^k \not\equiv 1 \ \text{mod} \ q \ for \ all \ k < 64.\\
    && \zeta = 12 \ \text{in Tiny DSA}.\\\\
& \lbrack a, b \rbrack & \text{For two integers} \ a \le b, \ [a, b] \ \text{denotes the set of integers} \ \{a, a+1, \ldots, b \}.\\\\
\end{array}
$$

## Examples of Polynomial Rings
Polynomial ring $Z[x]$ is the set of all polynomials with coefficients in $Z$ = {-$\infty$, $\ldots$, -1, 0, 1, $\ldots$, $\infty$}.

Polynomial ring $Z_q[x]$ is the set of all polynomials with coefficients in $Z_q$ = {0, 1, $\ldots$, $q-1$}.


*Example 1* - Polynomial $a(x) \in \mathbb{Z}[x]$.
$$
    a(x) = -12 + 20367x - 2081x^{56} + 10x^{1023}
$$

$\colorbox{gray}{}$


*Example 2* - Polynomial $b(x) \in \mathbb{Z}_q[x]$.

The coefficients of $a(x)$ from *Example 1* are mapped to $\mathbb{Z}_q$ in polynomial $b(x)$ below.

Keep in mind, Tiny DSA defines $q = 769$.

$$
  \begin{array}{lll}
    a(x) & = & -12 + 20367x - 2081x^{56} + 10x^{1023}\\
    b(x) & = & (-12 \ \text{mod} \ 769) \\
        &   & \, + \, (20367 \ \text{mod} \ 769) x \\
        &   & \, + \, (-(2081 \ \text{mod} \ 769) \ \text{mod} \ 769)x^{56} \\
        &   & \, + \, 10x^{1023} \\
        & = & \, (769-12) + 373x + (769-543)x^{56} + 10x^{1023}\\
        & = & \, 757 + 373x + 226x^{56} + 10x^{1023}\\
\end{array}
$$

$\colorbox{gray}{}$


## Examples of Polynomials in Tiny DSA
If polynomial $a(x) \in R[x]$, its coefficients are in the ring $\mathbb{Z}$ and its highest-degree term is at most $x^{31}$.

If polynomial $a(x) \in R_q[x]$, its coefficients are in the ring $\mathbb{Z}_q$ and its highest-degree term is at most $x^{31}$.

*Example 1* - Polynomial $a(x) \in R_q[x]$.

In the polynomial $a(x) = 7 + 23x^{31}$, first term is the lowest-degree term, and the second, the highest-degree term.

$\colorbox{gray}{}$

*Example 2* - Polynomial $a(x) \in R_q$.

Consider the polynomial $a(x) = 63 + 159x + 48x^2 + 746x^{28} + x^{30}$. Its coefficients are in $\mathbb{Z_q}$, and the highest-degree term is $x^{30}$. Therefore, $a(x) \in R_q$.

$\colorbox{gray}{}$

*Example 3* - Polynomial $a(x) \in R \ \text{and} \ a(x) \not\in R_q$.

Consider the polynomial $a(x) = 2163 - 169x + 1048x^2 - 1746x^{28} + 2x^{30} \in R$.
Note $a(x) \not\in R_q$ because ***not all coefficients*** are in $Z_q$. The coefficients of all terms except the last one are in $\mathbb{Z}$.

$\colorbox{gray}{}$

*Example 4* - Transform $a(x) \in R$ to $b(x) \in R_q$.

For a number $z \in \ \mathbb{Z}$, $z \ \text{mod} \ q$ maps it to $\mathbb{Z}_q$. Therefore, given a polynomial $a(x) \in R$, transforming its coefficients *mod 769* transforms the polynomial to $R_q$.

$$
\begin{array}{llll}
a(x) & = & 2163 - 169x + 1048x^2 - 1746x^{28} + 2x^{30} & \in R\\
b(x) & = & (2163 \ mod \ 769) \\
    &   & \, + \, (-169 \ mod \ 769) x \\
    &   & \, + \, (1048 \ mod \ 769)x^2 \\
    &   & \, + \, ((-(1746 \ mod \ 769) \ mod \ 769))x^{28} \\
    &   & \, + \, (2 \ mod \ 769) x^{30} \\
   & = & 625 + 600x + 279x^2 +  561x^{28} + 2x^{30} & \in R_q\\
\end{array}
$$

$\colorbox{gray}{}$


<a id="primitives-roots-of-unity"></a>

## Primitive Roots of Unity {#primitives-roots-of-unity}
Let $n = 2^k$, and let $q$ be a prime such that $q-1$ is divisible by $2n$.

$$
\begin{array}{rlll}
n & = \ 32 & = 2^5\\
2n & = \ 64 & = 2^6 \\
q & = \ 2^{10} - 2^{8} + 1 & = 769\\
q-1 & = \ 2^8(2^2 - 1) & = 768\\
(q-1)/2n & = \ 2^8(2^2 - 1)/2^6 & = 12\\
\end{array}
$$

Let $\zeta \in \mathbb{Z}_q$ be an element of order $2n$. The order of of $\zeta$ is the smallest positive integer $t$ such that $\zeta^t = 1$.

In Tiny DSA, $\zeta = 12$, $\zeta^{2n} \equiv 1  \ \text{mod} \ q$, and $\zeta^n \equiv -1 \ \text{mod} \ q$.

```{.number-lines .python}
# python3
n = 32
q = 769
z = 12
assert pow(z, 2*n, q) == 1
assert pow(z, n, q) == q-1
# following two expressions evaluate to the same value.
assert pow(-1, 1, q) == q-1
assert -1 % q == q-1
# this follows directly from the above expressions.
assert -1 % q == pow(-1, 1, q)
# for all k < 64, z^k != 1 mod q.
# in other words, there is no k in [1, 63] such that z^k = 1 mod q
assert not any([pow(z, k, q)==1 for k in range(1, 2*n)])
```

## NTT Definition

Let $a(x) \in R_q$.

Define NTT(a) = $\hat{a}$ = (a($\zeta$), a($\zeta^3$), a($\zeta^5$), $\ldots$, a($\zeta^{2n-1}$)) $\in Z^n_q$. NTT(a) is a ***polynomial evaluation*** of $a(x)$ at $\zeta$, $\zeta^3$, $\zeta^5$, $\ldots$, $\zeta^{2n-1}$.



<figure style="text-align: center;">
![The Irreducible Factors of $X^{32}+1$](./images/factors-of-reduct-poly.png)
</figure>

# References

V4: The Number-Theoretic Transform (NTT). $\copyright$ Alfred Menezes. August 2024.
<https://cryptography101.ca/wp-content/uploads/2024/12/V4-slides-Kyber-and-Dilithium.pdf>
