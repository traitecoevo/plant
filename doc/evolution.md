% Trait evolution and community assembly in the TREE model
% Daniel Falster; Rich FitzJohn

# Background

Our goal with the evolutionary element of TREE is to predict the trait mixture favoured by a given ecological model. By *ecological model* we mean a particular set of functions, environmental conditions and parameter values that together determine the population dynamics of the system. Given such a model, we want to know whether it can maintain a mixture of different types and if so, what mixture is likely to be favoured under competitive evolutionary assembly. Comparing the evolutionary mixtures obtained from different ecological models, i.e. models differing in functional form or parameter values, allows us to identify key processes maintaining diversity.

A full model of evolution would consider both the ecological forces generating selection (i.e effects of *phenotypes*) together with genetic mechanisms responsible for generating diversity (i.e. production of *genotypes*). Here we are primarily interested in the former, and as such assume that the relevant types can be generated through appropriate mutations. The theory of adaptive dynamics offers a conceptual framework for understanding how selection on phenotypes, driven by competitive interactions and population dynamics, shapes biological communities, in particular the distribution of quantitative traits [@geritz_evolutionarily_1998; @dieckmann_dynamical_1996; @dieckmann_can_1997; @dieckmann_adaptive_2007].

Consider a community of  species potentially differing from one another in a number of quantitative traits. Let:

- $N$ be the number of species,
- $K$ be the number of quantitative traits,
- $x_{ij}$ be the value of the $j^{\textrm{th}}$ trait in species $i$, and
- $y_{i}$ be the abundance of individuals in species $i$ .

The (phenotypic) composition of the community is then described by the vector
of trait values
\begin{equation}x=\left(x_{11}, \ldots, x_{1K}, \ldots, x_{N1} , \ldots, x_{NK}\right)\end{equation}

and their corresponding abundances

\begin{equation}y=(y_{1},\ldots,y_{N}).\end{equation}

Define $f(x^\prime,x, y)$ to be a scalar-valued function giving
the **fitness (per capita rate of increase)** for individuals with trait values
$x^\prime$ growing in an environment shaped by resident community $(x,y)$.

Because fitness is calculated on the basis of physiology and under competition from residents, it captures all of the ecology encapsulated in different *ecological models*.

A number of different *community assembly models* can then be used to simulate the evolutionary dynamics of the system. The models are distinguished by details such as the frequency and size of mutations. Below we adopt the terminology of @dieckmann_adaptive_2007. For the most part, the formulation of fitness function and details of community assembly model are largely orthogonal: any given ecological model can be implemented under a range of assembly models and vice versa.

# Polymorphic and stochastic model

If mutations occur frequently, evolution proceeds via **polymorphic and stochastic** model of trait evolution. In this model, clouds of residents drift and diffuse through trait space via selection and mutation.

Some advantages of the polymorphic and stochastic model are:

- that it involves the fewest assumptions and is therefore the most realistic: mutations are not rare in real communities
- that it involves the least code
- that it may be faster in high dimensional systems, with many traits or
species, as solving of demographic attractors becomes more difficult.

## Numerical approach

1.  PRIMING THE SYSTEM
	1.  Set parameters of ecological model, i.e. fitness function
	1.  Set parameters of the assembly model
		* Step size for population dynamics $\Delta t$
		* Parameters for **birth function** - see below
		* Parameters for **death function** - extinction threshold $\epsilon$
	1.  **Choose initial resident community**: with traits $x(0)=\left(x_{1}, \ldots, x_{N}\right)$ and corresponding densities $y(0)=\left(y_{1}, \ldots, y_{N}\right)$.
2. MAIN ASSEMBLY LOOP
	1. **Deaths**: Remove any strategy that has dropped below extinction threshold; i.e. check $y_i >  \epsilon$
	2. **Births**: Generate new types and introduce these into the community, if viable
		* Check fitness is $\geq$ 0, i.e. can invade
		* Increase dimension of community: $N \rightarrow N+1$
		* Add mutant at minimum density: $y_N = \epsilon$
	3. **Population dynamics**: Step population by $\Delta t$
		* Using difference equation for discrete time models: $y_i (t+ 1) = y_i(t) (1+ \exp(\Delta t \, f(x_i, x,y)))$
		* Using an adaptive ODE solver for continuous time models: $y_i (t+ \Delta t) = y_i(t) +\int_t^{t+\Delta t} y_i \, f(x_i, x,y)) \, dt.$
	4. (Optional) Generate **fitness landscape**
	5. Repeat above steps until **exit condition** is satisfied. Could be
		* Number of steps completed
		* Nothing alive and no regions of positive fitness
		* Stability: rate of change in traits, population abundance, or fitness landscape dropped below some threshold (but what about cyclic dynamics)

The **initial resident community** may be chosen using a variety of methods:

* Empty community ($N=0$)
* Random sample within some bounds ($N\geq1$)
* Values maximising fitness in an empty community ($N=1$)
* Single-species attractor ($N=1$) (see below for details).

The **birth function** is used to introduce new types into the community. The following are some different options

1. *Local mutations (1)*:
	* New types generated via mutations of residents at rate $\mu \, y_i$, the product of the per capita  mutation rate ($\mu$) and  population size ($y_i$).
	* The actual number of mutations in any given step is drawn from a poisson distribution with parameter $\mu \, y_i \, \Delta t$.
	* For each mutant, traits are given by drawing a mutation step $\Delta x$ from assumed distribution and adding this to resident trait values.
2. *Local mutations (2)*:
	* As above, but  community level mutation rate is held constant at some rate $M_c$. Thus $\mu$ is set such that  $M_c = \mu \sum_{i=1}^{i=N} \, y_i$.
3. *Immigration*:
	* New types are generated through arrival of seeds from outside the community at some rate $M_I$.
	* The actual number of immigrants in any given step is drawn from a poisson distribution with parameter $\M_I \, \Delta t$.
	* For each immigrant, traits are given by drawing from an assumed  distribution. For example, values might be drawn from a uniform distribution within certain bounds or a normal distribution.
4. *Most fit*:
	* New types are generated in proportion to fitness landscape.


## Monomorphic systems

When mutations are rare, ecological dynamics proceed much faster than evolutionary dynamics. Thus when a successful mutation arises, we can assume it spreads before any additional mutations appear. Under this assumption we have a formal time-scale separation between evolutionary and ecological processes, leading to the **monomorphic** models of adaptive dynamics, so called because each species is represented by a single strategy. The term **oligomorhpic*** is also used to describe communities with several species (*oligo-*, containing a relatively small number of units).

Another assumption about the size of mutation then distinguishes between
two modes of adaptive evolution that can be handled using the
monomorphic model. If mutations are large, evolution proceeds as directed 'random walk' through trait space; the **monomorphic and stochastic** model. On the other hand, if mutations are small, then evolution can be
modelled deterministically via gradient ascent, as a standard non-linear, dynamical system. This is the **monomorphic and deterministic** model.

Running any of the monomorphic models requires that we solve for the demographic attractor of the resident community. Let $\bar{y}$ be the value of $y$ satisfying

\begin{equation}\label{ybar}f(x_{i},x,\bar{y})=0 \, \forall i.\end{equation}

Thus when $y=\bar{y}$ residents are at their demographic attractors, also called the *carrying capacity*.

We can then write

\begin{equation} \hat{f}(x^\prime,x) = f(x^\prime,x, y)\bigg|_{y=\bar{y}} \end{equation}

to represent the fitness a mutant in a demographically stable resident population.

### Numerical approach
The following techniques may be applied to find  $\bar{y}$:

1.  **Analytical**: ideal, but only possible in some models and then only with single resident.
2.  **Multi-dimensional root solving**: This approach works well only if you have good initial guess for the solution. In systems with many species, obtaining a good initial guess is problematic.
3.  **Iteration**: Using a difference equation, $y_i (t+ 1) = y_i(t) (1+ f(x_i, x,y))$, iterate the population until stable. This approach is fail proof, but a large number of iterations may be required to reach stability, especially when fitness is close to zero.
4.  **Solve system of ODEs**: With $\frac{\partial}{dt} y_{i} = y_{i} \, f(x_i,x, y)$, advance the system using an adaptive ODE solver. Depending on the system a stiff solver may perform better, this in turn requires calculation of jacobian.

For the root solving and ODE approach, it may be preferable to make $\log y$ the state variable, as this prevents negative population sizes being generated by the solver.

## Monomorphic and stochastic model

This model assumes infrequent but large mutations; evolution then follows a directed random walk through trait space [@dieckmann_dynamical_1996; @dieckmann_adaptive_2007].  Mathematically, the model is described by a master equation for the probability density $P(x,t)$ of realising a given phenotypic
distribution $x$ at time $t$ [see @dieckmann_dynamical_1996]. The shape of the fitness landscape influences the probability per unit time that a species transitions from its current trait value to another trait value in that direction. Evolution within each species then follows a Markovian *trait substitution sequence*. Mutational processes interact with the selective forces generated by ecological interaction to determine evolutionary trajectories.

### Numerical approach

The challenge when using this model is that there can large periods of time when nothing happens, while we wait for a new mutation to arise. This suggests an algorithm with an adaptive step size will be needed. One approach is using Gillespie’s minimal process method
[@gillespie_general_1976; @dieckmann_dynamical_1996; @ito_new_2007; @brannstrom_emergence_2010]. Rather than
stepping the system over a given time interval and considering any
transitions that may occur in that time period, the system is stepped
between successive events by drawing the next invading mutant trait
value and the time at which the invasion occurs from appropriate
distributions. The method is founded on the idea that with a given
mutation rate, the distribution of waiting times follows an exponential
distribution.

The following algorithm is adapted from @ito_new_2007..... modify from above.


The above recipe estimates the waiting time between successive mutants,
and then asks whether such mutants can invade. It is also possible to
estimate the waiting time between successful mutants, using an extension
of this model [see @brannstrom_emergence_2010].
Although superior in many ways, this revised technique requires an
integration over the fitness landscape, which is costly in systems with
more than a single trait and is to be avoided.

## Monomorphic and deterministic model

When mutations are small, evolutionary trajectories approach their deterministic limits and can thus be approximated by a differential equation, leading to a gradient ascent model of phenotypic change. Although the dynamics of the system under this model are entirely deterministic, most evolutionary problems are inherently non-linear, meaning they cannot be solved analytically. Thus, the  dynamics must be analysed numerically  using standard techniques from dynamic systems research.

The canonical equation of adaptive dynamics [@dieckmann_dynamical_1996]captures the combined influences of selection and mutation on phenotypic evolution. Assuming small, random mutational steps, @dieckmann_dynamical_1996 show that the temporal evolution of traits in a large population can be approximated by the following canonical equation:

\begin{equation} \label{eq:CanEq} \frac{d}{dt} x_{i}(x)  = \frac{1}{2} \, \bar{y_{i}}  \, \mu(x_{i}) \, \sigma^{2}(x_{i})  \, g_{i}(x) \end{equation}

for $i=1,\ldots,N$, where $\mu(x_{i})$ is the mutation rate for species with traits $x_{i}$ , $\bar{y_{i}}$ is the species population size, $\sigma^{2}(x_{i})$ is the covariance matrix of the symmetric mutation distribution around $x_{i}$, and  $g_{i}(x)$ is the selection gradient. Let us define a new function, $g_{ij}(x)$, to be the selection pressure acting on trait $j$ in species $i$, given by the first-order partial derivative of $f$ with respect to trait $j$, evaluated at $x_{i}^{\prime} = x_{i}$:

\begin{equation} \label{eq:SelecGrad_i}  g_{ij}(x)  = \frac{\partial}{\partial x_{ij}^{\prime}} f(x_i^{\prime},x) \mid _{x_i^{\prime} =x_{i}}.\end{equation}

The selection gradient on all traits of species $i$ is then
\begin{equation} \label{eq:SelecGrad_x} g_{i}(x)= \left(\frac{\partial}{\partial x_{i1}^{\prime}} f(x_i^{\prime},x) \mid _{x_i^{\prime}=x_i}, \ldots, \frac{\partial}{\partial x_{iK}^{\prime}} f(x_i^{\prime},x) \mid _{x_i^{\prime}=x_i}\right)^{\top},\end{equation}

while the selection pressure on the entire multi-species communities is
\begin{equation} \label{eq:SelecGrad} g(x) = \left(g_{1}(x)^{\top}, \ldots , g_{N}(x)^{\top} \right)^{\top} = \left(g_{11}(x), \ldots, g_{1K}(x), \dots, g_{N1}(x), \ldots, g_{NK}(x)\right)^{\top}.\end{equation}

Note that $\frac{d}{dt}x_{i}$ and $g_{i}(x)$ are each vectors of length $K$, $\bar{y_{i}}$ and $\mu(x_{i})$ are scalars, and $\sigma^{2}(x_{i})$ is a $K\times K$ matrix representing the mutation distribution. For example, if each species has 2 traits, and mutations follow a bi-variate normal mutation distribution,
$$\sigma^{2}(x_{i}) = \begin{pmatrix} \sigma_{i1}^2 & \sigma_{i1}\sigma_{i2}\rho_{12}\\
\sigma_{i1}\sigma_{i2}\rho_{12} & \sigma_{i2}^2 \end{pmatrix},$$
where $\rho_{12}$ is chance of drawing mutations influencing trait $1$ and $2$ simultaneously.

The structure of the equation \ref{eq:CanEq} means that evolution of the whole community can be written in the form
\begin{equation} \label{eq:CE}
\frac{d}{dt} x  = \mathbf {B}(x) g(x),\end{equation}
where $\mathbf {B}$ is a block, diagonal, symmetric, positive-definite matrix with blocks $\mathbf{B_{ii}}(x_i) = \frac{1}{2} \, \bar{y_{i}}  \, \mu(x_{i}) \, \sigma^{2}(x_{i})$, comprising the mutational-matrix for each species. An advantage of the formulation given in equation \ref{eq:CE} is that it separates the two elements influencing evolution: (i) the mutational process $\mathbf {B}(x)$, and (ii) the selective forces determined by ecological interactions, $g(x)$. The matrix  $\mathbf {B}$ is block diagonal because there may be covariance in mutations within species but not between. For example, imagine a system with 3 species and 2 traits per species. $\mathbf {B}$ would be of type:
\begin{equation} \begin{split}
\mathbf {B}(x) & = \left[ \begin{array}{ccc} B_{11}(x_1) & 0  & 0 \\  0 & B_{22}(x_2)  & 0 \\  0 & 0& B_{33}(x_3)  \\ \end{array} \right]\\
& = \left[ \begin{array}{cccccc} \sigma_{11}^2 & \sigma_{11}\sigma_{12} & 0 & 0 &0 & 0 \\
			\sigma_{11}\sigma_{12} & \sigma_{22}^2 & 0 & 0 &0 & 0 \\
			0 & 0 & \sigma_{21}^2 & \sigma_{21}\sigma_{22} &0 & 0 \\
			0 & 0 & \sigma_{21} \sigma_{22} & \sigma_{22}^2 &0 & 0 \\
			0&0&0 & 0 & \sigma_{31}^2 & \sigma_{31}\sigma_{32} \\
			0&0&0 & 0 & \sigma_{31} \sigma_{32} & \sigma_{32}^2 \\ \end{array} \right]. \\
\end{split} \end{equation}
(NB: To simplify presentation of this matrix, the terms $\frac{1}{2} \, \bar{y_{i}}  \, \mu(x_{i})$ are not included in the matrix, but should be.)


The key difference between single and multi-dimensional systems is that mutational covariance can influence the course of evolution for the latter. In the simplest 2D case, when there is no mutational covariance between traits (i.e. $\sigma^{2}(x_{i})$ is an identity matrix, and the 2D mutational distribution corresponds to delta sheets along the trait axes), mutation rates and variance simply carry forward as a multiplier, adjusting the rate of evolution in each trait. However, with pleiotropic effects (where a single mutation influences multiple traits), the selective force on each trait is given by cumulative influence of several traits:
	$$\frac{d}{dt} x_{ij}(x)  = \sum_{k=1}^K (\mathbf{B_{ii}}(x_i))_{jk} \, g_{ik}(x).$$
In single-dimensional systems, or in multi-dimensional systems without pleiotropy, the mutational matrix $\mathbf B$ is an identify matrix. Factors influencing mutation still carry over as multiplication factors on the selection gradient, but the change in each trait is independent of the selective pressure on other traits.

## Numerical approach - root solving

Assuming we are dealing with a single species and a simple mutation matrix, we identify the end point for equation \ref{eq:CE} using root finding techniques. We want find $x^*$ such that $g_{ij}(x^*) \approx 0$. This can be achieved by first finding an interval bracketing $x^*$, i.e. finding $x_1, x_2$ s.t. $g_{ij}(x_1) > 0$ and $g_{ij}(x_2) < 0$. Then $x_1 < x^* < x_2$.

While root -solving techniques could in principle be extended to more than one dimension, the algorithms tend to be unstable.

## Numerical approach - Canonical equation

The system of ODEs described by equation \ref{eq:CE} can be advanced within a specified accuracy using a suitable ODE solver. For simple systems, with a single trait and few species, introductory techniques such as a Runge-Kutta solver with adaptive time step work well. All that is required from the ecological model is ability to calculate the fitness gradient for any arbitrary trait combination (eq. \ref{eq:SelecGrad_x}).

An algorithm for assembly using the **monomorphic and deterministic model** is as follows:

1. PRIMING THE SYSTEM
	1.  Set parameters of ecological model
	2.  Set parameters of the assembly model
		* Variance covariance matrix
		* Parameters for **birth function**
		* Parameters for **death function** - extinction threshold $\epsilon$
2. MAIN ASSEMBLY LOOP
	1. **Births**: Generate new resident and introduce into the community
	2. **Step canonical equation until stable**: Using an adaptive ODE solver step equation \ref{eq:CE} until either stable or one species goes extinct
	3. Repeat above steps until **exit condition** is satisfied. Could be
		* Time point reached
		* Stability: rate of change in traits, population abundance, or fitness landscape dropped below some threshold (but what about cyclic dynamics)

In practice several difficulties may arise when running this model

* **Solving demographic attractor is slow**
* **Stiff systems:** In multi-trait systems, however, the possibility of a stiff set of equations arises. Stiffness occurs when there are two very different time scales on which the variables are changing. The solution space of a stiff ODE has a slow manifold, on which the state point moves slowly, and off which the state point moves rapidly toward the slow manifold.  The downside of using a stiff solver is that they (mostly) require the first derivative (Jacobian) of the selection gradient to be calculated, which requires many different resident populations to be simulated.

The Jacobian of the selection gradient is given by,
\begin{equation} \label{eq:SelecGrad-Jacobian1}
	J_{g_i}(x) = \frac{\partial}{\partial x_i} g_i(x) = \begin{pmatrix} \frac{\partial}{\partial x_{i1}} g_{i1}(x) \mid _{x_i^{\prime}=x_i} & \ldots & \frac{\partial}{\partial x_{iK}} g_{i1}(x)\mid _{x_i^{\prime}=x_i} \cr
		\ldots & \ldots &\ldots \cr
		\frac{\partial}{\partial x_{i1}} g_{iK}(x) \mid _{x_i^{\prime}=x_i}& \ldots & \frac{\partial}{\partial x_{iK}} g_{iK}(x)\mid _{x_i^{\prime}=x_i} \end{pmatrix}
\end{equation}
where the elements of the Jacobian are given by
\begin{equation} \label{eq:SelecGrad-Jacobian2}
	\frac{\partial}{\partial x_{ij}} g_{ik}(x)  = \frac{\partial^2}{\partial x^{\prime}_{ik} \partial x^{\prime}_{ij}} f(x_i^{\prime},x) \mid _{x_i^{\prime}=x_i} + \frac{\partial^2}{\partial x^{\prime}_{ik} \partial x_{ij}} f(x_i^{\prime},x) \mid _{x_i^{\prime}=x_i}.
\end{equation}
(Calculated using the chain rule: $\frac{\partial}{\partial x_{ij}} g_{ik}(x)  = \frac{\partial}{\partial x^{\prime}_{ij}} g_{ik}(x) \frac{\partial x^{\prime}_{ij}}{\partial x_{ij}} + \frac{\partial}{\partial x_{ij}} g_{ik}(x) \frac{\partial x_{ij}}{\partial x_{ij}}$).

In practice, the selection gradient (eq \ref{eq:SelecGrad}) and its Jacobian (eqs \ref{eq:SelecGrad-Jacobian2} \& \ref{eq:SelecGrad-Jacobian2}) are calculated numerically, for example using the finite difference method.


# References
