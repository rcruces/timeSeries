Logistic map: discrete dynamical systems
========================================

Most realistic systems with chaotic behavior, such as fluid flow around
an obstacle, are described by a non-linear partial differential
equation, which determines how the velocity changes in space and time.
We will see later in the course that chaos also occurs in many ordinary
(non-linear) differential equations. Here we will consider something
even simpler, an “iterative map”. This means that starting from an
initial value of a variable, *x*<sub>0</sub> say, we generate a sequence
of values, *x*<sub>1</sub>, *x*<sub>2</sub>, etc. from the map
(i.e. function), *x*<sub>*n* + 1</sub> = *f*(*x*<sub>*n*</sub>) where we
here make a simple choice
*f*(*x*) = 4*λ**x*<sub>*n*</sub>(1 − *x*<sub>*n*</sub>)
where *λ* is parameter. In other words,
*x*<sub>1</sub> = 4*λ**x*<sub>0</sub>(1 − *x*<sub>0</sub>),
*x*<sub>2</sub> = 4*λ**x*<sub>1</sub>(1 − *x*<sub>1</sub>), etc. We will
be interested in the behavior of successive iterations of this map, as a
function of the parameter *λ*. In particular we will study the behavior
of the *x*<sub>*n*</sub> for large n.  
We consider *λ* in the range from 0 to 1, so, if *x*<sub>0</sub> is
between 0 and 1, it is easy to see that all subsequent values of x also
lie in this range. In fact the largest value of *x*<sub>*n* + 1</sub>
(which occurs for *x*<sub>*n*</sub> = 1/2) is equal to *λ*.  
This so-called “logistic map” has been used as model for population
dynamics, but here we just treat it as a toy model which has a
transition to chaos.

The logistic map function for population growing is define following:  
*X*<sub>*n*</sub> = *r**x*(1 − *x*)
Where *x* ∈ {0, 1}, represents the ratio of existing population to the
maximum possible population (1).  
*X*<sub>*n*</sub> is the new population after *n* generations.  
*r* is the combined rate between reproduction and mortality.

Population for different *r* values
-----------------------------------

![](logistic_map_files/figure-markdown_github/population-1.png)

![](logistic_map_files/figure-markdown_github/logisticmap-var-1.png)![](logistic_map_files/figure-markdown_github/logisticmap-var-2.png)

Poincaré phase diagram (self-similarity recurrence plot)
--------------------------------------------------------

![](logistic_map_files/figure-markdown_github/phase-1.png)![](logistic_map_files/figure-markdown_github/phase-2.png)![](logistic_map_files/figure-markdown_github/phase-3.png)

Recurrence plot of the time series
----------------------------------

<a href="http://www.recurrence-plot.tk/glance.php" class="uri">http://www.recurrence-plot.tk/glance.php</a>
![](logistic_map_files/figure-markdown_github/recur-1.png)
![](logistic_map_files/figure-markdown_github/recur-2.png)
![](logistic_map_files/figure-markdown_github/recur-4.png)

Further lectures:
-----------------

<a href="https://geoffboeing.com/2015/03/chaos-theory-logistic-map/" class="uri">https://geoffboeing.com/2015/03/chaos-theory-logistic-map/</a>  
<a href="http://www.kierandkelly.com/from-chaos-to-creativity/" class="uri">http://www.kierandkelly.com/from-chaos-to-creativity/</a>  
<a href="http://www.kierandkelly.com/what-is-chaos/logistic-map/" class="uri">http://www.kierandkelly.com/what-is-chaos/logistic-map/</a>

Extra: r for negative values
----------------------------

![](logistic_map_files/figure-markdown_github/logisticmapLong-1.png)

![](logistic_map_files/figure-markdown_github/logisticmapNeg-1.png)![](logistic_map_files/figure-markdown_github/logisticmapNeg-2.png)

![](logistic_map_files/figure-markdown_github/logisticmapMirror-1.png)
