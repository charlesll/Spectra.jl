# Baseline subtraction

Baseline subtraction can be made with using the baseline function, which is the main API using other dedicated baseline functions. The documentation of all those functions is available below.

```@docs
baseline
```

The `baseline` function actually calls more dedicated functions, in particular for the automatic baselines. 

```@docs
als_baseline
arPLS_baseline
drPLS_baseline
```

The `rubberband` method does not work well for now but here is its docstrings for reference:

```@docs
rubberband_baseline
```


