\newcommand{\ud}{\mathrm{d}}

# Equations

For the sake of brevity, dependencies of functions are shown only in the Symbol column. The variable $x$ refers to a vector of four traits that are varied in the model: $x = (\phi, \rho, h_m, s)$. Subscripts for size variables are: l = leaves, s = sapwood, b = bark and phloem, h = heartwood, r = fine roots. All mass measurements are in terms of dry mass.

## Traits

### Leaf mass per unit area (LMA)

Symbol: $\varphi$

Units: kg m$^{-2}$

### Stem tissue density
Symbol: $\rho$

Units: kg m$^{-3}$

### Height at maturation

Symbol: $h_m$

Units: m

### Seed size

Symbol: $s$

Units: kg

## Individual size

### Height

Symbol: $h(x, a)$

Units: m

Determination: Plant height is obtained by solving the IVP with initial conditions: $h(x, a) = h_0$, and with rates of change: $\frac{\textrm{d}h}{\textrm{d}t} = g(x, h(x, a))$. Thus, $h(x, a) = h_0 + \int_{0}^{a} g(x, h(x, a^\prime)) \, \rm{d}a^\prime$. The initial plant height at germination, $h_{0}$ is obtained by finding a value that satisfies the equation $s=m_\textrm{t}(x, h)$.

### Leaf area

Symbol: $a_\textrm{l}(x, h)$

Units: m$^2$

Determination: $a_\textrm{l}=\alpha_1 \, h^{\beta_1}$, $\frac{\textrm{d}h}{\textrm{d}a_\textrm{l}}= -\beta_1\big(\frac{a_\textrm{l}}{\alpha_1}\big)^{-(\beta_1+1)}$


### Area of sapwood

Symbol: $a_{ss}(x, a_\textrm{l})$

Units: kg

Determination: $a_\textrm{ss}=\theta^{-1} \, a_\textrm{l}$, $\frac{\textrm{d}a_\textrm{ss}}{\textrm{d}t} =\theta^{-1} \, \frac{\textrm{d}a_\textrm{l}}{\textrm{d}t}$

### Mass of sapwood

Symbol: $m_{ss}(x, a_\textrm{l})$

Units: kg

Determination: $m_\textrm{sb}=\theta^{-1} \, \rho \, \eta_c \, a_\textrm{l} \, h$, $\frac{\textrm{d}m_\textrm{s}}{\textrm{d}a_\textrm{l}}=\theta^{-1}\, \rho\, \eta_c\, \big( h + a_\textrm{l}\, \frac{\textrm{d}h}{\textrm{d}a_\textrm{l}} \big)$

### Area of bark

Symbol: $a_{sb}(x, a_\textrm{l})$

Units: kg

Determination: $a_\textrm{sb}= b \, \theta^{-1} \, a_\textrm{l}$, $\frac{\textrm{d}a_\textrm{sb}}{\textrm{d}t} =b \, \theta^{-1} \, \frac{\textrm{d}a_\textrm{l}}{\textrm{d}t}$

### Mass of bark

Symbol: $m_{sb}(x, a_\textrm{l})$

Units: kg

Determination: $m_\textrm{sb}=\theta^{-1} \, \rho \, \eta_c \, a_\textrm{l} \, h$, $\frac{\textrm{d}m_\textrm{sb}}{\textrm{d}a_\textrm{l}}=\theta^{-1}\, \rho\, \eta_c\, \big( h + a_\textrm{l}\, \frac{\textrm{d}h}{\textrm{d}a_\textrm{l}} \big)$

### Area of heartwood

Symbol: $a_\textrm{sh}(x, a)$

Units: kg

Determination: Area of heartwood is obtained by solving the IVP with initial conditions: $a_\textrm{sh} =0$, and  rates of change
$\frac{\textrm{d}a_\textrm{sh}}{\textrm{d}t}=k_\textrm{s} \, a_\textrm{ss}$. Thus, $a_\textrm{sh}=\int_0^t \frac{\textrm{d}a_\textrm{sh}}{\textrm{d}t}(t^\prime) \, dt^\prime$.

### Mass of heartwood

Symbol: $m_\textrm{sh}(x, a)$

Units: kg

Determination: Mass  of heartwood is obtained by solving the IVP with initial conditions: $m_\textrm{sh} =0$, and  rates of change $\frac{\textrm{d}m_\textrm{sh}}{\textrm{d}t}=k_\textrm{s} \, m_\textrm{ss}$. Thus, $m_\textrm{sh}=\int_0^t \frac{\textrm{d}m_\textrm{sh}}{\textrm{d}t}(t^\prime) \, dt^\prime$

### Mass of fine roots

Symbol: $m_\textrm{r}(x, a_\textrm{l})$

Units: kg

Determination: $m_\textrm{r} = \alpha_3 a_\textrm{l}(x, h)$, $\frac{\textrm{d}m_\textrm{r}}{\textrm{d}a_\textrm{l}}= \alpha_3$

### Total mass

Symbol: $m_\textrm{t}(x, h)$

Units: kg

Determination: $m_\textrm{t} = m_\textrm{l}+m_\textrm{ss}+m_\textrm{sb}+m_\textrm{sh}+m_\textrm{r}$

## Competitive environment

### Probability density of leaf area at height z for an individual of height h

Symbol: $q(z, h)$

Units: m$^{-1}$

Determination: $q = 2\eta(1-z^\eta h^{-\eta}) z^{\eta-1} h^{-\eta}$ if $z \leq h$, otherwise 0.

### Fraction of leaf area above height z for an individual of height h

Symbol: $Q(z, h)$

Units: Dimensionless (0 to 1)

Determination: $Q = \int_z^h q(z^\prime,h) \; \textrm{d}z^\prime$ if $z \leq h$, otherwise 0.

Equation: 10

### Canopy openness at height z in a patch of age a

Symbol: $E(z, a)$

Units: Dimensionless (0 to 1)

Determination: $E = \exp\left(-c_\textrm{ext} \int_0^\infty Q(z, h(m_l)) \;a_\textrm{l}(x, h)  \;n(x,h,a) \; \textrm{d}m_l \right)$

Equation: 11

## Mass production

All rates are per plant.

### Gross annual CO2 assimilation

Symbol: $A(x,h_\textrm{l}, E(\cdot,a))$

Units: mol yr$^{-1}$

Determination: $A =a_\textrm{l}(x, h) \int_0^{h(m_l)} A_{\textrm{lf}}(A_0 \upsilon,E(z,a)) \; q(z, h(m_l)) \; \textrm{d} z$

Equation: 12

$A_{\textrm{lf}}(A_0\upsilon,E(z,a))$ is the gross annual CO2 assimilation per unit leaf area at canopy openness $E(z, a)$ for a leaf with maximum capacity $A0\upsilon$, determined by integrating instantaneous rates of assimilation (described by a rectangular hyperbola) over the diurnal solar cycles throughout the year.

### Total maintenance respiration

Symbol: $R(x,h})$

Units: mol yr$^{-1}$

Determination: $R = a_\textrm{l}(x, h) \upsilon c_\textrm{R,l} + \frac{m_s(x,h) + 2 m_b(x,h)}{\rho}  c_\textrm{R,s} +m_r(x, h)  c_\textrm{R,l}$

Equation: 13

### Total turnover

Symbol: $T(x,h)$

Units: kg yr$^{-1}$

Determination: $T=m_l\alpha_4 \phi^{\beta_4} + m_b(x, m_l) k_b +m_r(x,m_l)k_r$

Equation: 14

### Net production

Symbol: $P(x, h,E(\cdot, a))$

Units: kg yr$^{-1}$

Determination: $P = c_{\textrm{bio}}Y \left[ A\left(x,h,E(\cdot, a)\right) - R(x,h)\right]-T(x,h)$

Equation: 15

### Fraction of production allocated to reproduction

Symbol: $r(x,h)$

Units: Dimensionless (0 to 1)

Determination: $r = c_{\textrm{r},1}\left(1+\exp\left(c_{\textrm{r},2} \left(1-\frac{h(x,m_l)}{h_m}\right)\right)\right)^{-1}$

Equation: 16

### Rate of offspring production

Symbol: $f(x, h, E(\cdot ,a))$

Units: yr$^{-1}$

Determination: $f = \frac{r(x,h) \;P(x,h,E(\cdot, a))}{c_{\textrm{acc}}s}$ if $P(x,m_{\textrm{l},0},E(\cdot, a))>0$, otherwise 0

Equation: 17


### Fraction of whole-plant growth that is leaf

Symbol: $\frac{\textrm{d}m_l}{\textrm{d}m_t} (x, m_l)$

Units: Dimensionless (0 to 1)

Determination: $\frac{\textrm{d}m_l}{\textrm{d}m_t} =  \frac{ 1}{1+  \frac{\textrm{d}m_s}{\textrm{d}m_l} (x, m_l) +  \frac{\textrm{d}m_b}{\textrm{d}m_l} (x, m_l) +  \frac{\textrm{d}m_h}{\textrm{d}m_l} (x, m_l)+  \frac{\textrm{d}m_r}{\textrm{d}m_l} (x, m_l)}$

Equation: 18

Notes: The derivatives on the right-hand side of eqn 18 can be calculated directly from eqn 4–7. For solutions see Appendix S2.

### Growth rate in leaf mass

Symbol: $g(x,h,E(\cdot ,a))$

Units: kg yr$^{-1}$

Determination:  $g =(1-r(x, m_l)) \; P(x, m_l, E(\cdot, a)) \; \frac{\textrm{d}m_l}{\textrm{d}m_t} (x, m_l)$ if $P(x,m_{\textrm{l},0},E(\cdot, a)) > 0$, otherwise 0.

Equation: 19

## Mortality

### Survival of seedlings during germination

Symbol: $\Pi_1(x, m_{\textrm{l},0}, E(\cdot, a))$

Units: Dimensionless (0 to 1)

Determination: $\Pi_1 = \left(1+ \left(\frac{c^2_{\textrm{s},0}}{\frac{P(x, m_{\textrm{l},0} E(\cdot, a))}{a_\textrm{l}(x, h_0)}} \right)^{2} \right)^{-1}$ if $P(x,m_{\textrm{l},0},E(\cdot, a)) > 0$, otherwise 0

Equation: 20

### Instantaneous mortality rate

Symbol: $d(x, h, E(\cdot, a))$


Units: yr$^{-1}$


Determination: $d = c_{\textrm{d}0}\exp\left(-c_{\textrm{d}1}\rho\right) + c_{\textrm{d}2} \exp\left(-c_{\textrm{d}3} \frac{P(x, m_l, E(\cdot, a))}{a_\textrm{l}(x, h)} \right)$

Equation: 21

### Development of size distribution within patches

### Density per ground area of individuals with traits x and size h in a patch of age a

Symbol: $n(x,h, a)$

Units: kg$^{-1} m$^{-2}$

Determination:

$\frac{\partial}{\partial a} n(x,h, a) =  - d(x, m_l, E(\cdot, a))\;n(x,h,a) - \frac{\partial}{\partial m_l}  \left[ g(x, m_l, E(\cdot, a))\;n(x,h,a) \right]$

$ n(x,h, 0) = 0$

$n(x, h_0, a) =  \frac{pi_1(x, m_{\textrm{l},0}, E(\cdot, a))}{ g(x, m_{\textrm{l},0}, E(\cdot, a))} \int_0^\infty p(\tau) \int_0^\infty \Pi_0  f(x, m_l,  E(\cdot, \tau))\;n(x, h, \tau)\;\textrm{d}h \textrm{d}\tau$

if $P(x,m_{\textrm{l},0},E(\cdot, a)) > 0$, otherwise 0

Equation: 22

## Metapopulation dynamics

### Probability density of patch age a in the metapopulation

Symbol: $p(a)$

Units: yr$^{-1}$

Determination: $p= \frac{1}{\hat a} \exp\left(- \frac{\pi}{4} \left(\frac{a}{\hat a }\right)^2 \right)$

Equation: 23

Notes: $\hat a$ is the mean interval between disturbances. The probability of patch disturbance is assumed to increase linearly with patch age, and can be expressed as a function of mean disturbance interval, $\gamma(a) = \pi \frac{a}{2\hat a^2}$. For more details see Appendix S1.

### Emergent properties of vegetation

Averages over all patches in the metapopulation, calculated as inline image, where K(a) is the considered vegetation property at patch age a.

### Average height of leaf area

Symbol: $H(a)$

Units: m

Determination: $H = \frac{1}{L(a)} \int_0^\infty \int_0^\infty a_\textrm{l}(x, h) q(z, h(x,m_l)) n(x, h, a)  \;\textrm{d}h \textrm{d}z$

Equation: 24

### Leaf-area index

Symbol: $L(a)$

Units: Dimensionless

Determination: $L = \int_0^\infty a_\textrm{l}(x, h) n(x, h, a)  \;\textrm{d}h$

Equation: 25

### Net primary production

Symbol: $N(a)$

Units: kg m$^{-2}$ yr$^{-1}$

Determination: $N =c_{\textrm{bio}}Y  \int_0^\infty \left[ A\left(x,h,E(\cdot, a)\right) - R(x,h)\right] n(x, h, a)  \;\textrm{d}h$

Equation: 26

### Biomass density

Symbol: $B(a)$

Units: kg m$^{-2}$

Determination: $B = \int_0^\infty m_\textrm{t}(h) n(x, h, a)  \;\textrm{d}h$

# Appendix 1: Moving from leaf mass to height as key variable

The focus on leaf mass as the key variable in @falster_influence_2011 is fairly arbitrary.
Instead, consider a change of variables to use "height" as the key
size variable.  This is useful because we need to know about plant
height elsewhere, and it's a common measuring stick that is fairly
independent of the key traits.

To convert from leaf mass ($m_l$) to height ($h$), use equations 3 and
4 in @falster_influence_2011:

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

$$
1 +
% d m_s / d m_l + d m_b / d m_l
(1 + b) \frac{\rho \eta_c}{\theta\phi}\alpha_1(\beta_1+1)
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

Other useful derivatives:

Because
$$
h = \alpha_1 \left(\frac{m_l}{\phi}\right)^{\beta_1}
\qquad
m_l = \phi\left(\frac{h}{\alpha_1}\right)^{1/\beta_1}
$$

we have

$$
\frac{\ud h}{\ud m_l} = \frac{\alpha_1 \beta_1}{\phi}
\left(\frac{m_l}{\phi}\right)^{\beta_1 - 1}
\qquad
\frac{\ud m_l}{\ud h} =
\frac{\phi}{\alpha_1 \beta_1}\left(\frac{h}{\alpha_1}\right)^
{\frac{1}{\beta_1} - 1}
$$

# References
