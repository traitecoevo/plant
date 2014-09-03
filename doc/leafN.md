% Leaf N in competitive canopies
% Daniel Falster

\newcommand{\ud}{\mathrm{d}}
\newcommand{\tm}{\textrm}

# Background

# Model

We want to calculate total assimilation for a single plant $A_c$, which is obtained by integrating the instantaneous assimilation rates $A$ over the height $h$ of the plant:

\begin{equation}\label{eq:Ac} A_\tm{c} = a_\tm{l} \, \int_0^h A\left(I(z), A_\tm{max}(z)\right) \, q(z) \, \ud z. \end{equation}

Here $I(z), A_\tm{max}(z), q(z)$ are the light intensity, maximum leaf photosynthetic capacity per area, and density at height $z$.
Light assumed to follow Beer-Lambert law
\begin{equation} I(z) = I_0 \exp(-c_\tm{ext} \, L(z)),\end{equation}
where $L(z)$ is cumulative leaf area above height $z$ and $-c_\tm{ext}$ is the light extinction coefficient. It will be helpful below to use the substitution
\begin{equation} E(z) = \exp(-c_\tm{ext} \, L(z)),\end{equation}
such that
\begin{equation}\label{eq:I} I(z) = I_0 E(z).\end{equation}
Moreover, the average canopy openness within the crown is
\begin{equation}\label{eq:Ebar}\bar{E}(h) = \int_0^h E(z) q(z) \, \ud z \end{equation}

Following @sands_modelling_1995, it can be shown (using calculus of variations) that if the plant has a given amount nitrogen to allocate towards photosynthetic proteins, the optimal distribution is for leaf N to be distributed in proportion to light, i.e.
\begin{equation}\label{eq:Nz} N(z,h) = N(h) \frac{E(z)}{E(h)},\end{equation}
where $N(h)$ is the nitrogen per leaf area in the top leaf of a plant height $h$. The average nitrogen per unit leaf area within the crown is then

\begin{equation}\bar{N} = \int_0^h N(z) q(z) \, \ud z.\end{equation}
Combing with eqs \ref{eq:Ebar} and \ref{eq:Nz}, we can solve for $N(h)$ as
\begin{equation}N(h) = \bar{N}\frac{E(h)}{\bar{E}(h)}.\end{equation}

This gives

\begin{equation}\label{eq:Nz} N(z,h) = \bar{N}\frac{E(z)}{\bar{E}(h)}.\end{equation}

So in the case where there is no light gradient then there is also no gradient in leaf N, i.e. $E(z) = \bar{E}(h) \rightarrow N(z) = \bar{N}$.

Let's now look at assimilation. Following Sands 1995, we use a non-rectangular hyperbola to model leaf assimilation:

\begin{equation}\label{eq:PLRC1} \theta \, A^2 - (\alpha I + A_\tm{max}) A + \alpha I A_\tm{max}=0,\end{equation}

where $A_\tm{max}$ is the photosynthetic capacity, $\alpha$ is the quantum yield, and $\theta$ describes the curvature of the response. This has solution

\begin{equation}\label{eq:PLRC} A = \frac{\alpha I + A_\tm{max}+ \sqrt{(\alpha I + A_\tm{max})^2 - 4 \theta \alpha I A_\tm{max}}}{2\theta}.\end{equation}

In addition we assume a linear relation between maximum photosynthetic rate and leaf N:
\begin{equation} A_\tm{max} = c_\tm{N} \, N.\end{equation}

Substituting from equation \ref{eq:Nz}, we obtain an equation for the change in photosynthetic capacity through the canopy:

\begin{equation} \label{eq:amax} A_\tm{max}(z) = c_\tm{N} \, \bar{N}\frac{E(z)}{\bar{E}(h)}. \end{equation}

Substituting eqs. \ref{eq:I} and \ref{eq:amax} into \ref{eq:PLRC} then gives the following for instantaneous assimilation:

\begin{equation}\label{eq:PLRC2} A(I(z), N(z)) = E(z) A_0(\bar{N},\bar{E}_\tm{h}),\end{equation}

where
\begin{equation}\label{eq:A0} A_0(\bar{N}, \bar{E}) = \frac{\alpha I_0 + c_\tm{N} \frac{\bar{N}}{\bar{E}} + \sqrt{(\alpha I_0 + c_\tm{N} \frac{\bar{N}}{\bar{E}})^2 - 4 \theta \alpha I_0 c_\tm{N} \frac{\bar{N}}{\bar{E}}}}{2\theta}\end{equation}
is the photosynthesis of a leaf with leaf nitrogen of $\frac{\bar{N}}{\bar{E}}$ in full light. Whole plant assimilation then becomes

\begin{equation}\begin{split}\label{eq:Ac2} A_\tm{c} &= a_\tm{l} \, A_0(\bar{N},\bar{E}_\tm{h}) \int_0^h E(z) \, q(z) \, \ud z \\&= a_\tm{l} \, A_0(\bar{N},\bar{E}_\tm{h}) \, \bar{E}_\tm{h}. \end{split}\end{equation}

which is beautifully simple, given what is being achieved.

So, the recipe for estimating photosynthesis (eq. \ref{eq:Ac}), for a plant with known $\bar{N}$ and $h$ and given the function $E$, is

1. Calculate $\bar{E}_\tm{h}$ via eq. \ref{eq:Ebar}
2. Calculate $A_0(\bar{N}, \bar{E})$ via eq. \ref{eq:A0}
3. Calculate $A\tm{c}$ via eq \ref{eq:Ac2}.


## Next steps

1. Show that solution in \ref{eq:Ac2} works, by comparing to numerical integration of equations \ref{eq:Ac}, \ref{eq:I}, \ref{eq:Ebar}, \ref{eq:amax}.
2. Sow that commonly-used whole canopy expression as special case of above equations
3. Calculate surface of $A_0$ with respect to $\bar{E}, \bar{N}$, assess how fast it is to calculate
4. Look at canopy gradients.

# References
