% Details about derivations
\newcommand{\ud}{\mathrm{d}}

# Moving from leaf mass to height as key variable

The focus on leaf mass as the key variable is fairly arbitrary.
Instead, consider a change of variables to use "height" as the key
size variable.  This is useful because we need to know about plant
height elsewhere, and it's a common measuring stick that is fairly
independent of the key traits.

To convert from leaf mass ($m_l$) to height ($h$), use equations 3 and
4 in the paper:

$$
h = \alpha_1 \left(\frac{m_l}{\phi}\right)^{\beta_1}
$$

where $\phi$ is the leaf mass per unit area (\textsc{lma}) and
$\alpha_1$ and $\beta_1$ are scaling parameters relating leaf area
($m_l/\phi$) to height.

What we need for the main calculations is the derivative of height $h$
with respect to time.  Because $m_l$ is a function of time, we have

$$
\frac{\ud h}{\ud t} = \frac{\ud}{\ud t}[m_l(t)]
\frac{\alpha_1 \beta_1 \left(\frac{m_l(t)}{\phi}\right)^{\beta_1 - 1}}{\phi}
$$

where $\frac{\ud}{\ud t}[m_l(t)]$ comes from equation 19:

$$
\frac{\ud}{\ud t}[m_l(t)] = (1 - r)P\frac{\ud m_l}{\ud m_t}
$$

where $P$ is total production, $r$ is allocation to reproduction and
$\ud m_l/\ud m_t$ is the fraction of whole plant growth that is leaf.



