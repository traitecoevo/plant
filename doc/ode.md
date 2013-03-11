# ODE solver

The GSL ode solvers (etc) are going to be tedious to use, as our
problems regularly change size.  We really don't want to have to free
and reallocate all the memory at every step, as that will be very slow
(based on similar experience with same in the adaptive spline
calculations).

We could fake it, by having spare capacity but setting the derivatives
for this spare space as zero.  This is also going to be ugly and we'll
need to abstract that away anyway.

The basic ODE solver approach needs to do the following things:

- set a state (in our case this might be a partial set with just new
  individuals added).
- step: take one step with the problem (this might take several
  attempts if the step size is too large)

# First pass

The first pass is a non-templated, basic system.  It does not do any
clever allocation, nor any cleverl switching between steppers.  The
basic idea is just to get some integration working.

I have lumped the only stepper into Step.{cpp,h}, but this could
easily be generalised, and we could switch between different steppers
defined in a series of classes.  That is probably not something that I
will do without packaging this separately though.

I have not formalised the bits about being able to use `dydt_in` and
gives exact `dydt_out`, yet.

## Ownership

At the moment, `Control` is fully self contained, but `Step` contains
the ODE derivatives function; this is probably OK, though it makes it
a `Step<Problem>`, so `Evolve` will also become an `Evolve<Problem>`,
too.

However, if all Steps were templated, (say `StepRKCK<T>` and
`StepRK4<T>`):

```c++
	template <class Problem>
	class StepRKCK {
		StepRKCK(Problem *pr_);
		...
		Problem *pr_; // defining pr->derivs(), pr->size()
	};
```

Then we can do this:
```c++
	template <class Stepper>
	class OdeSolver {
		OdeSolver(Stepper st);
		...
		Stepper st;
	}
```
and then this

```c++
	RKCK<MyProblem> stepper(my_instance);
	OdeSolver< RKCK<MyProblem> > solver(stepper);
```

and if I look at odeint I can probably see a neater way of setting
that up.  Probably a nice intermediate amount of typing is.

``` c++
	Odesolver< RKCK<MyProblem> > solver(stepper);
```
which might require
```c++
	template<class Stepper, Problem>
	void OdeSolver<Stepper, Problem>::OdeSolver(Problem *p) {
		Stepper<Problem> s(p);
		OdeSolver(s);
	}
```

but at the least this might be possible:
```c++
	Odesolver< RKCK<MyProblem> > solver(stepper(stepper));
```

None of this is an issue for a while though.




## Divergence from GSL

I think that my Evolve will be where we diverge from GSL (especially
as I'm not that interested in the higher level driver object).  As a
result, I end up with something of a hybrid between `Evolve` and
`Driver`, which I will probably call `Ode` or `OdeSolver`.

Extra logic that is required:
  min and max step size
  current and end time
  

driver::apply (multiple step size between t0 and t1)
driver::fixed_step (one single step)

## Error checking

The GSL version has very heavy error checking (after basically every
operation, it checks to see if things worked OK).  We don't do that.
This means that if the derivative function returns NaN values, etc,
we're going to fail terribly, quite possibly.

## Clear modification rules

In the GSL, a ton of stuff is passed as pointers and modified.
Getting this into an OO style is giving me grief.  This is especially
true in `Evolve::apply()`, which I need to go through carefully at
some point.  h0, t0, t, etc.  All horribly named.

It seems that there is a 'state' object, containing:
  * `y` -- state
  * `dydt_in` -- derivative of `y` wrt time at beginning of step
  * `dydt_out` -- computed derivative of y wrt time
  * `yerr` -- error estimate
could be passed by reference much more happily.  

`Step::apply` would fill this up; call as step.apply(t, h, state) with
definition `Step::apply(double t, double h, State &y);` The (a)
problem with this is that if the step fails we do `state.y = y0` and
we are relying on `dydt_in` not changing.  Perhaps just references
(const and not)?
  
## Iterators for the derivatives function

This is probably going to seem really weird to people, but it exists
here because we are going to need it in our model.  Not sure how to
make this optional, but iterators are at least fairly common, and the
intention is fairly clear.
