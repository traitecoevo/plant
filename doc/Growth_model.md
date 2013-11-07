

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

Determination: $s=m_\textrm{t}(x, m_{\textrm{l},0})$

Equation: 1

Notes: Leaf mass at germination, $m_{\textrm{l},0}$ is obtained by finding a value that satisfies equation 1 and varies as a function of $\phi$, $\rho$ and $s$.

## Individual size

### Mass of leaves

Symbol: $m_\textrm{l}$

Units: kg

### Leaf area

Symbol: $\omega(x, m_\textrm{l})$

Units: m$^2$

Determination: $\omega=m_\textrm{l}\varphi^{-1}$

Equation: 2

### Height

Symbol: $h(x, m_\textrm{l})$

Units: m

Determination: $h=\alpha_1 \omega(x, m_\textrm{l})^{\beta_1}$

Equation: 3

### Mass of sapwood

Symbol: $m_s(x, m_\textrm{l})$

Units: kg

Determination: $m_s = \rho \eta_c \theta^{-1} \omega(x, m_\textrm{l}) h(x, m_\textrm{l})$

Equation: 4

### Mass of bark

Symbol: $m_\textrm{b}(x, m_\textrm{l})$

Units: kg

Determination: $m_\textrm{b} = b m_s(x, m_\textrm{l})$

Equation: 5

### Mass of heartwood

Symbol: $m_\textrm{h}(x, m_\textrm{l})$

Units: kg

Determination: $m_h = \rho \eta_c \alpha_2 \omega(x, m_\textrm{l})^{\beta_2}$

Equation: 6

### Mass of fine roots

Symbol: $m_\textrm{r}(x, m_\textrm{l})$

Units: kg

Determination: $m_\textrm{r} = \alpha_3 \omega(x,m_\textrm{l})$

Equation: 7

### Total mass

Symbol: $m_\textrm{t}(x, m_\textrm{l})$

Units: kg

Determination: $m_\textrm{t} = m_\textrm{l}+m_s+m_\textrm{b}+m_\textrm{h}+m_\textrm{r}$

Equation: 8

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

Determination: $E = \exp\left(-c_\textrm{ext} \int_0^\infty Q(z, h(m_l)) \;\omega (x, m_l)  \;n(x, m_l,a) \; \textrm{d}m_l \right) $

Equation: 11

## Mass production

All rates are per plant.

### Gross annual CO2 assimilation

Symbol: $A(x,m_\textrm{l}, E(\cdot,a))$

Units: mol yr$^{-1}$

Determination: $A = \omega(x,m_l) \int_0^{h(m_l)} A_{\textrm{lf}}(A_0 \upsilon,E(z,a)) \; q(z, h(m_l)) \; \textrm{d} z$

Equation: 12

$A_{\textrm{lf}}(A_0\upsilon,E(z,a))$ is the gross annual CO2 assimilation per unit leaf area at canopy openness $E(z, a)$ for a leaf with maximum capacity $A0\upsilon$, determined by integrating instantaneous rates of assimilation (described by a rectangular hyperbola) over the diurnal solar cycles throughout the year.

### Total maintenance respiration

Symbol: $R(x,m_\textrm{l})$

Units: mol yr$^{-1}$

Determination: $R = \omega(x,m_l) \upsilon c_\textrm{R,l} + \frac{m_s(x,m_l) + 2 m_b(x,m_l)}{\rho}  c_\textrm{R,s} +m_r(x,m_l)  c_\textrm{R,l}$

Equation: 13

### Total turnover

Symbol: $T(x,m_\textrm{l})$

Units: kg yr$^{-1}$

Determination: $T=m_l\alpha_4 \phi^{\beta_4} + m_b(x, m_l) k_b +m_r(x,m_l)k_r$

Equation: 14

### Net production

Symbol: $P(x, m_\textrm{l},E(\cdot, a))$

Units: kg yr$^{-1}$

Determination: $P = c_{\textrm{bio}}Y \left[ A\left(x,m_\textrm{l},E(\cdot, a)\right) - R(x,m_\textrm{l})\right]-T(x,m_\textrm{l})$

Equation: 15

### Fraction of production allocated to reproduction

Symbol: $r(x,m_\textrm{l})$

Units: Dimensionless (0 to 1)

Determination: $r = c_{\textrm{r},1}\left(1+\exp\left(c_{\textrm{r},2} \left(1-\frac{h(x,m_l)}{h_m}\right)\right)\right)^{-1}$

Equation: 16

### Rate of offspring production

Symbol: $f(x, m_\textrm{l}, E(\cdot ,a))$

Units: yr$^{-1}$

Determination: $f = \frac{r(x,m_\textrm{l}) \;P(x,m_\textrm{l},E(\cdot, a))}{c_{\textrm{acc}}s}$ if $P(x,m_{\textrm{l},0},E(\cdot, a))>0$, otherwise 0

Equation: 17


### Fraction of whole-plant growth that is leaf

Symbol: $\frac{\textrm{d}m_l}{\textrm{d}m_t} (x, m_l)$

Units: Dimensionless (0 to 1)

Determination: $\frac{\textrm{d}m_l}{\textrm{d}m_t} =  \frac{ 1}{1+  \frac{\textrm{d}m_s}{\textrm{d}m_l} (x, m_l) +  \frac{\textrm{d}m_b}{\textrm{d}m_l} (x, m_l) +  \frac{\textrm{d}m_h}{\textrm{d}m_l} (x, m_l)+  \frac{\textrm{d}m_r}{\textrm{d}m_l} (x, m_l)}$

Equation: 18

Notes: The derivatives on the right-hand side of eqn 18 can be calculated directly from eqn 4–7. For solutions see Appendix S2.

### Growth rate in leaf mass

Symbol: $g(x,m_\textrm{l},E(\cdot ,a))$

Units: kg yr$^{-1}$

Determination:  $g =(1-r(x, m_l)) \; P(x, m_l, E(\cdot, a)) \; \frac{\textrm{d}m_l}{\textrm{d}m_t} (x, m_l)$ if $P(x,m_{\textrm{l},0},E(\cdot, a)) > 0$, otherwise 0.

Equation: 19

## Mortality

### Survival of seedlings during germination

Symbol: $\Pi_1(x, m_{\textrm{l},0}, E(\cdot, a))$

Units: Dimensionless (0 to 1)

Determination: $\Pi_1 = \left(1+ \left(\frac{c^2_{\textrm{s},0}}{\frac{P(x, m_{\textrm{l},0} E(\cdot, a))}{\omega(x, m_{l,0})}} \right)^{2} \right)^{-1}$ if $P(x,m_{\textrm{l},0},E(\cdot, a)) > 0$, otherwise 0

Equation: 20

### Instantaneous mortality rate

Symbol: $d(x, m_\textrm{l}, E(\cdot, a))$


Units: yr$^{-1}$


Determination: $d = c_{\textrm{d}0}\exp\left(-c_{\textrm{d}1}\rho\right) + c_{\textrm{d}2} \exp\left(-c_{\textrm{d}3} \frac{P(x, m_l, E(\cdot, a))}{\omega(x, m_l)} \right) $

Equation: 21

### Development of size distribution within patches

### Density per ground area of individuals with traits x and size m_\textrm{l} in a patch of age a

Symbol: $n(x, m_\textrm{l}, a)$

Units: kg$^{-1} m$^{-2}$

Determination:

$\frac{\partial}{\partial a} n(x, m_\textrm{l}, a) =  - d(x, m_l, E(\cdot, a))\;n(x, m_l,a) - \frac{\partial}{\partial m_l}  \left[ g(x, m_l, E(\cdot, a))\;n(x, m_l,a) \right] $

$ n(x, m_\textrm{l}, 0) = 0$

$n(x, m_{\textrm{l},0}, a) =  \frac{pi_1(x, m_{\textrm{l},0}, E(\cdot, a))}{ g(x, m_{\textrm{l},0}, E(\cdot, a))} \int_0^\infty p(\tau) \int_0^\infty \Pi_0  f(x, m_l,  E(\cdot, \tau))\;n(x, m_l, \tau)\;\textrm{d}m_l \textrm{d}\tau$

if $P(x,m_{\textrm{l},0},E(\cdot, a)) > 0$, otherwise 0

Equation: 22

## Metapopulation dynamics

### Probability density of patch age a in the metapopulation

Symbol: $p(a)$

Units: yr$^{-1}$

Determination: $p= \frac{1}{\hat a} \exp\left(- \frac{\pi}{4} \left(\frac{a}{\hat a }\right)^2 \right) $

Equation: 23

Notes: $\hat a$ is the mean interval between disturbances. The probability of patch disturbance is assumed to increase linearly with patch age, and can be expressed as a function of mean disturbance interval, $\gamma(a) = \pi \frac{a}{2\hat a^2}$. For more details see Appendix S1.

### Emergent properties of vegetation

Averages over all patches in the metapopulation, calculated as inline image, where K(a) is the considered vegetation property at patch age a.

### Average height of leaf area

Symbol: $H(a)$

Units: m

Determination: $H = \frac{1}{L(a)} \int_0^\infty \int_0^\infty \omega(x,m_l) q(z, h(x,m_l)) n(x,m_l,a) \;\textrm{d}m_l \textrm{d}z$

Equation: 24

### Leaf-area index

Symbol: $L(a)$

Units: Dimensionless

Determination: $L = \int_0^\infty \omega(x,m_l) n(x,m_l,a) \;\textrm{d}m_l$

Equation: 25

### Net primary production

Symbol: $N(a)$

Units: kg m$^{-2}$ yr$^{-1}$

Determination: $N =c_{\textrm{bio}}Y  \int_0^\infty \left[ A\left(x,m_\textrm{l},E(\cdot, a)\right) - R(x,m_\textrm{l})\right] n(x,m_l,a) \;\textrm{d}m_l$

Equation: 26

### Biomass density

Symbol: $B(a)$

Units: kg m$^{-2}$

Determination: $B = \int_0^\infty m_\textrm{t}(m_\textrm{l}) n(x,m_l,a) \;\textrm{d}m_l$
