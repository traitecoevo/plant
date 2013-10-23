## Converting density from one size unit to another

Population density is explicitly modelled in relation to a given size unit (Eq  $\ref{eq:PDE}$). But what if we want to express density in relation to another size unit? A relation between the two can be derived by noting that the total number of individuals within a given size range must be equal. So let's say density is expressed in units of size $m$, but we want density in units of size $h$. First we require a one-to-one function which h for a given $m$: $h = \hat{h}(m)$. Then the following must hold
\begin{equation} \label{eq:n_conversion} \int_{m_1}^{m_2} n(x,m,a) \; \textrm{d}m =  \int_{\hat{h}(m_1)}^{\hat{h}(m_2)} n^\prime(x,m,a) \; \textrm{d}h  \end{equation}

For very small size intervals, this equation is equivalent to
\begin{equation} \left(m_2- m_1 \right) \; n(x,m_1,a) = \left( \hat{h}(m_2) - \hat{h}(m_1)\right) \; n^\prime(x, \hat{h}(m_1),a). \end{equation}

Rearranging gives
\begin{equation}  n^\prime(x, \hat{h}(m_1),a) = n(x, m_1,a) \; \frac{m_2- m_1}{\hat{h}(m_2) - \hat{h}(m_1)}  \end{equation}

Noting that the second term on the RHS is simply the definition of $\frac{\delta m}{\delta h}$ evaluated at $m_1$, we have
\begin{equation} \label{eq:n_conversion2} n^\prime(x, h, a) = n(x, m,a) \; \frac{\delta m}{\delta h}. \end{equation}
