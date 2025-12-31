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

# The Polynomial Remainder Theorem
We will develop the constructive proof briefly mentioned in the wikipedia article on [polynomial remainder theorem](https://en.wikipedia.org/wiki/Polynomial_remainder_theorem).

Let
\begin{equation}
\begin{array}{lll}
    f(x) & = & a_n \, x^n + a_{n-1} \, x^{n-1} + \dots + a_1 \, x + a_0
\end{array}
\end{equation}

be a polynomial with coefficients in a ring, and let $r$ be any element of the ring. Then there exists a polynomial $Q(x)$ of degree at most $(n-1)$ such that

\begin{equation}
\begin{array}{lll}
    f(x) & = & (x - r) \ Q(x) \ + \ f(r).
\end{array}
\end{equation}

The remainder of the division of $f(x)$ by the linear polynomial $(x - r)$ is exactly $f(r)$.

## Proof

### The polynomial and its evaluation at $r$
Let

\begin{equation}
\begin{array}{lll}
    f(x) & = & a_n x^n + a_{n-1} x^{n-1} + \dots + a_1 x + a_0.\\
\end{array}
\end{equation}

Evaluate $f(r)$:

\begin{equation}
\begin{array}{lll}
    f(r) & = & a_n r^n + a_{n-1} r^{n-1} + \dots + a_1 r + a_0.\\
\end{array}
\end{equation}

### The difference $f(x) - f(r)$

\begin{equation}
\begin{array}{lll}
f(x) - f(r) & = & (\color{blue}{a_n x^n + a_{n-1} x^{n-1} + \dots + a_1 x + a_0})
                    - (\color{Mulberry}{a_n r^n + a_{n-1} r^{n-1} + \dots + a_1 r + a_0}). \ \ \ \ \ \ \ \ \\
            & = & ({\color{blue}{a_n x^n}} - {\color{Mulberry}a_n r^n}) + \,
                  ({\color{blue}a_{n-1} x^{n-1}} - {\color{Mulberry}a_{n-1} r^{n-1}}) + \, \ldots
                  + \, ({\color{blue}a_1 x} - {\color{Mulberry}a_1 r})
                  + \, ({\color{blue}a_o} - {\color{Mulberry}a_0}) \\
            & = & \underbrace{a_n(x^n - r^n)}_{\text{TERM} \ n}
                    + \, \underbrace{a_{n-1}x^{n-1} - r^{n-1}}_{\text{TERM} \ n-1}
                    + \ \dots
                    + \, \underbrace{a_1(x - r)}_{\text{TERM 1}}
\end{array}
\end{equation}


### Recall the algebraic identity
\begin{equation}
\begin{array}{lll}
    x^k - r^k & = & (x - r)\bigl(x^{k-1} + x^{k-2}r + \dots + xr^{k-2} + r^{k-1} \bigr) \ \ \ for \ an \ integer \ k \ge 1. \\
    a_k (x^k - r^k) & = & a_k\,(x - r)\,\bigl(x^{k-1} + x^{k-2}r + \dots + xr^{k-2} + r^{k-1} \bigr) \ \ \ for \ k = 1, \dots, n.\\
\end{array}
\end{equation}

### Apply the identity to each term in $f(x) - f(r)$

Consider TERM 5 - $a_n(x^n - r^n)$, in (5) above and substitute
\begin{equation}
   a_k(x^k - r^k) = (x - r)\,a_k\bigl( x^{k-1} + x^{k-2}r + \dots + r^{k-1} \bigr) \ \ \ 1 \le k \le n.
\end{equation}

Rewriting the terms in (5) gives
\begin{equation}
\begin{array}{lll}
    a_n(x^n - r^n) & = & (x - r)\,a_n(x^{n-1} + x^{n-2}r + \dots + r^{\,n-1} )
    \\
    a_n(x^{n-1} - r^{n-1}) & = & (x - r)\,a_{n-1}( x^{n-2} + x^{n-3}r + \dots + r^{\,n-2} )
    \\
    \vdots
    \\
    a_1(x - r) & = & (x - r)\,a_1
\end{array}
\end{equation}


### Factor out $(x - r)$.

\begin{equation}
\begin{array}{ll}
        f(x) - f(r) &= (x - r)\,a_n( x^{n-1} + x^{n-2}r + \dots + r^{\,n-1} ) \\
                    &\quad + (x - r)\,a_{n-1}( x^{n-2} + x^{n-3}r + \dots + r^{\,n-2} ) \\
                    &\quad + \dots + (x - r)\,a_1.\\
\end{array}
\end{equation}


Factor $(x - r)$ from the entire sum (the right-hand side expression):

\begin{equation}
    f(x) - f(r) = (x - r)\Bigl(a_n( x^{n-1} + x^{n-2}r + \dots + r^{\,n-1}) + \dots + a_1 \Bigr).
\end{equation}

### Define the polynomial $Q(x)$ and Complete the Proof
Let
\begin{equation}
    Q(x) = a_n( x^{n-1} + x^{n-2}r + \dots + r^{\,n-1} )
      + a_{n-1}( x^{n-2} + x^{n-3}r + \dots + r^{\,n-2} )
      + \dots + a_1.
\end{equation}
where
\begin{equation}
    \deg(Q) \le n - 1
\end{equation}

We use polynomial $Q$, and write

\begin{equation}
    f(x) - f(r) = (x - r) \, Q(x).
\end{equation}

from which it follows that

\begin{equation}
   f(x) = (x - r) \, Q(x) + f(r).
\end{equation}

$\colorbox{gray}{}$