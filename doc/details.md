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
\frac{\alpha_1 \beta_1}{\phi} \left(\frac{m_l(t)}{\phi}\right)^{\beta_1 - 1}
$$

where $\frac{\ud}{\ud t}[m_l(t)]$ comes from equation 19:

$$
\frac{\ud}{\ud t}[m_l(t)] = (1 - r)P\frac{\ud m_l}{\ud m_t}
$$

where $P$ is total production, $r$ is allocation to reproduction and
$\ud m_l/\ud m_t$ is the fraction of whole plant growth that is leaf.

So overall, we have

$$
\frac{\ud h}{\ud t} =
\frac{\ud m_t}{\ud t} \frac{\ud m_l}{\ud m_t} \frac{\ud h}{\ud m_l}
$$

---

Here, I'm just explicitly writing down the derivations of the pieces
needed for this equation, and for the equivalent equation used in the
R version of the growth model.

From equation (18) we have

$$
\frac{\ud m_l}{\ud m_t} =
\frac{1}{1 +
	\frac{\ud m_s}{\ud m_l} +
	\frac{\ud m_b}{\ud m_l} +
	\frac{\ud m_h}{\ud m_l} +
	\frac{\ud m_r}{\ud m_l}}
$$

(because $m_t = m_l + m_s + m_b + mh + m_r$, and $\ud m_l/\ud m_t = 1
/(\ud m_t / \ud m_l)$).

Using appendix S2 of Falster et al. (2010).

$$
\frac{\ud m_s}{\ud m_l} = \frac{\ud}{\ud m_l}\left[
\frac{\rho \eta_c}{\theta} \omega(m_l) h(m_l)\right] =
\frac{\ud}{\ud m_l}\left[
\frac{\rho \eta_c}{\theta} \frac{m_l}{\phi}
\alpha_1\left(\frac{m_l}{\phi}\right)^{\beta_1}\right] =
\frac{\ud}{\ud m_l}\left[
\frac{\rho \eta_c}{\theta}
\alpha_1\left(\frac{m_l}{\phi}\right)^{\beta_1 + 1}\right] =
% r n / t (m/p) a (m/p) ^ b
% D[r n / t (m/p) a (m/p) ^ b, m] // Simplify
\frac{\rho \eta_c}{\theta\phi}\alpha_1(\beta_1+1)
\left(\frac{m_l}{\phi}\right)^{\beta_1}
$$
and
$$
\frac{\ud m_b}{\ud m_l} = \frac{\ud}{\ud m_l}[b m_s] =
b \frac{\ud m_s}{\ud m_l} =
b \frac{\rho \eta_c}{\theta\phi}\alpha_1(\beta_1+1)
\left(\frac{m_l}{\phi}\right)^{\beta_1}
$$
and
$$
\frac{\ud m_h}{\ud m_l} =
\frac{\ud}{\ud m_l}\left[
\rho\eta_c\alpha_2 \omega(m_l)^{\beta_2}
\right] =
\frac{\ud}{\ud m_l}\left[
\rho\eta_c\alpha_2 \left(\frac{m_l}{\phi}\right)^{\beta_2}
\right] =
\frac{\rho\eta_c\alpha_2}{\phi}\beta_2
\left(\frac{m_l}{\phi}\right)^{\beta_2-1} =
% D[r n a (m / p)^b, m]
% D[r n a (m / p)^b, m] // Simplify
\frac{\rho\eta_c\alpha_2}{m_l}\beta_2
\left(\frac{m_l}{\phi}\right)^{\beta_2}
$$
and
$$
\frac{\ud m_r}{\ud m_l} =
\frac{\ud}{\ud m_l}\left[\alpha_3 \omega(m_l)\right] =
\frac{\ud}{\ud m_l}\left[\alpha_3 \frac{m_l}{\phi}\right] =
\frac{\alpha_3}{\phi}
$$

So, the denomininator of $\ud m_l / \ud m_t$ looks like:

$$
1 +
% d m_s / d m_l
\frac{\rho \eta_c}{\theta\phi}\alpha_1(\beta_1+1)
\left(\frac{m_l}{\phi}\right)^{\beta_1} +
% d m_b / d m_l
b \frac{\rho \eta_c}{\theta\phi}\alpha_1(\beta_1+1)
\left(\frac{m_l}{\phi}\right)^{\beta_1} +
% d m_h / d m_l
\frac{\rho\eta_c\alpha_2}{m_l}\beta_2
\left(\frac{m_l}{\phi}\right)^{\beta_2} +
% d m_r / d m_l
\frac{\alpha_3}{\phi}
$$

Repeating this with the version in the R verison, working in terms of
leaf area.
$$
\frac{\ud \omega}{\ud m_t} =
\frac{1}{
	\frac{\ud m_l}{\ud \omega} +
	\frac{\ud m_s}{\ud \omega} +
	\frac{\ud m_b}{\ud \omega} +
	\frac{\ud m_h}{\ud \omega} +
	\frac{\ud m_r}{\ud \omega}}
$$

$$
\frac{\ud m_l}{\ud \omega} =
\frac{\ud}{\ud \omega}\left[\omega\phi\right] =
\phi
$$
and
$$
\frac{\ud m_s}{\ud \omega} =
\frac{\ud}{\ud \omega}\left[
\frac{\rho \eta_c}{\theta} \omega h(\omega)\right] =
\frac{\ud}{\ud \omega}\left[
\frac{\rho \eta_c}{\theta} \alpha_1 \omega^{\beta_1+1}\right] =
\frac{\rho \eta_c}{\theta}\alpha_1(\beta_1+1)\omega^{\beta_1}
$$
and
$$
\frac{\ud m_b}{\ud \omega} =
\frac{\ud}{\ud \omega}\left[b m_s\right] =
b \frac{\ud m_s}{\ud \omega} =
b \frac{\rho \eta_c}{\theta}\alpha_1(\beta_1+1)\omega^{\beta_1}
$$
and
$$
\frac{\ud m_h}{\ud \omega} =
\frac{\ud}{\ud \omega}\left[\rho\eta_c\alpha_2 \omega^{\beta_2}
\right] =
\rho\eta_c\alpha_2\beta_2\omega^{\beta_2 - 1}
$$
and
$$
\frac{\ud m_r}{\ud \omega} =
\frac{\ud}{\ud \omega}\left[\alpha_3 \omega\right] = \alpha_3
$$
