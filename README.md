\documentclass[11pt]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath, amssymb, amsthm}
\usepackage{hyperref}
\usepackage{microtype}
\usepackage{booktabs}

\title{\textbf{Sofyen's Theorem on Universal Non-Ordinarity \\ for the Genus-3 Family $y^2 = x^7 + \sigma x + 1$}}
\author{\textbf{Mohamed Sofyen} \\ 
Chief Technical Officer, Office National de la Télédiffusion \\ 
Tunis, Tunisia \quad \text{---} \quad \texttt{mohamedsofiene@gmail.com}}
\date{September 2026}

\begin{document}

\maketitle

\begin{abstract}
This repository contains the official theoretical framework, mathematical proofs, and multi-engine computational certificates establishing the \textbf{Sofyen Theorem} for the one-parameter family of genus-3 hyperelliptic curves $C_\sigma: y^2 = x^7 + \sigma x + 1$ over a prime field $\mathbb{F}_p$ ($p$ odd) [0.1.1]. 
The computational framework guarantees that the Cartier--Manin matrix determinant $\det M_p(\sigma)$ vanishes identically in $\mathbb{F}_p[\sigma]$ if and only if $p$ belongs to the finite exceptional set $\mathcal{E} = \{3, 7, 11, 23\}$ [0.1.1, 0.1.2].
\end{abstract}

\hrule
\vspace{0.4cm}

\section{Project Architecture \& Methodology}
To bypass the exponential slowdown of direct algebraic expansions (which typically stall around $p \sim 150$--$400$ [0.1.4]), this work operates on a modular, four-part scientific architecture [0.1.4]:
\begin{itemize}
    \item \textbf{Reachability Analysis (Principle 1):} An $O(1)$ structural check that determines if a matrix entry possesses any non-zero $\sigma$-term without evaluating the full polynomial [0.1.4].
    \item \textbf{Extremal-Exponent Forcing (Principle 2):} Isolates the unique dominant permutations in the Leibniz expansion of the determinant to certify non-vanishing directly [0.1.4].
    \item \textbf{The $p$-adic Gamma Shortcut (Principle 3):} Reduces complex multinomial coefficient ratios within a residue class $a = p \pmod 7$ to a single rational invariant $R_a \in \mathbb{Q}$, valid for all large primes simultaneously [0.1.4].
    \item \textbf{Log-Retained Replication (Principle 4):} Multi-engine replication across independent environments with complete retained execution outputs [0.1.4].
\end{itemize}

\section{Repository Structure}
The repository is strictly divided between documentation and its cross-verified computational engines [0.1.1]:
\begin{itemize}
    \item \texttt{/paper/} -- Academic preprint, algebraic derivations, and tables of $R_a$ constants.
    \item \texttt{/code/python/} -- Pure Python/SymPy standalone certificate performing the exhaustive $p < 10^8$ structural scan and exact evaluations [0.1.3].
    \item \texttt{/code/sage\_local/} -- First SageMath engine reproducing the $R_a$ rational table and crossing results with raw multinomial extractions up to $p < 20,000$ [0.1.3].
    \item \texttt{/code/sage\_global/} -- Independent, generic-$d$ alternative SageMath script generalizing the framework architecture to higher-order trinomials $x^d + \sigma x + 1$ [0.1.4, 0.1.12].
\end{itemize}

\section{Academic Citation}
Until a formal digital object identifier (DOI) is registered, scholarly works, dissertations, or software distributions relying on these proofs or computational invariants should use the following bibliographic reference [0.1.3, 0.1.5]:

\begin{quote}
\small
\textbf{Sofyen, Mohamed (2026).} \emph{The Four Exceptional Primes: The Sofyen Theorem on Universal Non-Ordinarity for the Genus-3 Family $y^2 = x^7 + \sigma x + 1$}. Complete Edition --- Theorem, Method, and Source Code [0.1.1].
\end{quote}

\section{Licensing \& Attribution}
\begin{itemize}
    \item \textbf{Mathematical Content \& Text:} Governed by the \texttt{Sofyen Theorem Academic License v1.0} (see the root \texttt{LICENSE} file). This license explicitly protects the authorship and historical attribution of the proof and formulas while preserving the freedom of independent mathematical rediscovery.
    \item \textbf{Software \& Scripts:} All computational implementations found within the \texttt{/code/} directories are dual-licensed under the highly permissive \textbf{MIT License} to encourage adaptation, checking, and redistribution within larger computer algebra suites [0.1.1].
\end{itemize}

\end{document}
