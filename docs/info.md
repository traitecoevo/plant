
# The Plant and Strategy objects

## Naming of variables

We use the followig conventions for name of varaibles:

- Specific variables are decsirbed by the measure then tissue, separated by an underscore,
e.g. `mass_leaf`, `area_sapwood`, `diameter_stem`
- Partial derivatives are written as `dX_dY' where `X,Y` are varaible names, e.g.
`dmass_leaf_darea_stem`
- time derivatives (rates) are written as `X_dt` where `X` is the varaible whsoe
rate of change we are calculating, e.g. `height_dt`, `area_leaf_dt`
