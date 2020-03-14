Implementation
===============

This section will introduce some algorithms use in spatialtis.

Determine cell shape
--------------------

There are two way to determine the shape of a cell, "convex" or "concave". The "concave" is slower and need an extra
parameter, an alpha value. This will influence the `cell_shape` and `area` but no other geometry information.

.. image:: ../src/convex_concave.png
    :align: center
    :width: 60%

Cell co-occurrence
------------------

The occurrence or absence of a cell type in a ROI (Region of interest) is determined by a user-defined threshold value,
for example if we set the threshold at 50, the cell type is determined as occurrence if it's counts exceed 50 otherwise
it's absence. Afterwards, if two types of cell are both determine as occurrence, one-way :math:`\chi^2` test are
conducted to determine the significance, if one or both are absence would be determined as non co-occurrence.

Find cell neighbors
-------------------

For each ROI, every cell will be stored in a R-tree structure allow for fast neighbor search. In spatialtis we use
shapely's `STRtree <https://shapely.readthedocs.io/en/latest/manual.html#str-packed-r-tree>`_.

Neighborhood analysis
----------------------

This analysis is from `histocat <https://www.nature.com/articles/nmeth.4391>`_, it's a Monte Carlo method. After `Find cell neighbors <#find-cell-neighbors>`_,
we know the neighbor cells for each cell. Therefore, we know the counts of cell type `X` at the neighborhood of cell type `Y`. Then we randomize the cell in a ROI, the position of cell and the sum of every
cell type will not change, only cell type will be reassigned. Then we count the neighbors again. User can define how many times randomization (resample) will happen. And then a pseudo `P` value is calculated:

.. math::
    P_{association} = \frac{\text{Numbers of }(\overline{perm}\geq\overline{real})}{\text{Number of resample} + 1}

    P_{avoidance} = \frac{\text{Numbers of }(\overline{perm}\leq\overline{real})}{\text{Number of resample} + 1}

With user defined p-value threshold, we can determine the significance.


Spatial enrichment analysis
---------------------------

This analysis is from `MIBI paper <https://www.nature.com/articles/nm.3488>`_, it's basically the same as the neighborhood analysis, the major
difference is they used Z-score instead of a `P` value. Positive Z-score means association and negative Z-score means avoidance. You might notice that
in that paper they do it based on markers, because they defined different cell types based on marker expression and assign the type name based
on marker name. This is also possible in spatialtis simply by storing these info in a new key in `anndata.obs` field, and then tell
spatialtis the `type_col` has changed to a new key.

Spatial distribution
---------------------

There are three point distribution patterns in general, random, regular and cluster. Random means the point pattern follows the poisson process,
the regular means evenly distributed and cluster means the points tend to aggregate. (Cells are represented by their centroid)

.. image:: ../src/distribution_pattern.png
    :align: center
    :width: 50%

To determine the cell distribution patterns in each ROI, spatialtis provided three methods.

     - Index of Dispersion (ID)
     - Morisita’s index of dispersion (MID)
     - Clark and Evans aggregation index (CE)

+--------------------------------------+--------+---------+---------+
|                                      | Random | Regular | Clumped |
+======================================+========+=========+=========+
| Index of dispersion: ID              | ID = 1 | ID < 1  | ID > 1  |
+--------------------------------------+--------+---------+---------+
| Morisita’s index of dispersion: I    | I = 1  |  I < 1  |  I > 1  |
+--------------------------------------+--------+---------+---------+
| Clark and Evans aggregation index: R | R = 1  |  R < 1  |  R > 1  |
+--------------------------------------+--------+---------+---------+

Index of dispersion
###################

.. figure:: ../src/index_of_dispersion.png
    :width: 50%
    :align: center
    :figclass: align-center

    Sampling process, the orange circle is the sampling windows, the number is the count of points

First we store all the point in a ROI in KD tree. A random sample window with diameter `r` is generated,
the count of points in this window is `x`, a number of counts are generated after sampling many times. The null
hypothesis is that the points are randomly distributed. :math:`s^2` is the variance of all samples,
:math:`\overline{x}` is the average of all samples. Index of dispersion is calculated
as follow:

.. math:: ID = \frac{s^2}{\overline{x}}


Morisita’s index of dispersion
##############################

This is a quadratic statistic method, user need to define how to rasterize the ROI.

.. figure:: ../src/quadratic_statistic.png
    :width: 50%
    :align: center
    :figclass: align-center

    In this example, the ROI is divided into :math:`3\times3` grids, the number is the count of points

Morisita’s index of dispersion is calculated as follow:

.. math:: I_d = n[\frac{\sum x^2 - \sum x}{(\sum x)^2 - \sum x}]

:math:`\sum x` sum of the quadrat counts :math:`\sum x = x_1+x_2+x_3+...`

:math:`\sum x^2` sum of quadrat counts squared :math:`\sum x = x_1^2+x_2^2+x_3^2+...`

:math:`\chi^2 = I_d(\sum x - 1)+n-\sum x`  (:math:`df = n-1`)

Clark and Evans aggregation index
##################################

This method evaluate the distribution pattern base on distance between points. The points are stored in KD tree
at the first place.

Index of aggregation is calculated as follow:

.. math::
    R = \frac{\overline{r}_A}{\overline{r}_E}

:math:`\overline{r}_A` Mean distance to nearest neighbor: :math:`\overline{r}_A = \frac{\sum r_i}{n}`

:math:`r_i` Distance to nearest neighbor for individual :math:`i` (here we use euclidean distance)

:math:`n` number of individuals

:math:`\overline{r}_E` Expected distance to nearest neighbor: :math:`\overline{r}_E = \frac{1}{2\sqrt{\rho}}`

:math:`\rho` density of individuals: :math:`\rho = \frac{n}{\text{area size}}`

:math:`z = \frac{\overline{r}_A - \overline{r}_E}{S_r}`

:math:`S_r` Standard error of the expected distance to nearest neighbor: :math:`S_r = \frac{0.26136}{\sqrt{n\rho}}`


Spatial heterogeneity
----------------------

In spatialtis, we use Shannon entropy to quantify the heterogeneity in a ROI

.. math::
    H(X) = -\sum P_i log_2(P_i)

To compare the difference within a group (eg. different samples from same tumor),
Kullback–Leibler divergences for each sample within the group are computed, smaller value indicates less difference within group.

.. math::
    D = \sum P_i log_2(\frac{P_i}{Q_i})

Hotspot detection
------------------

Hotspot detection is used to find the cells that form clumps. Here we use
Getis–Ord hotspot analysis. First we rasterize the ROI into grids, for each small
square, we will compare it to its neighbor cells. User can define the level of neighbors
to search.

.. image:: ../src/hotspot_search.png
    :width: 50%
    :align: center

z score for a region :math:`i`:

.. math::
    z_i=\sum_{j=1}^n W_{i,j} C_j - \frac{\overline{C}\sum_{j=1}^n W_{i,j}}{SU}

.. math::
    S=\sqrt{\frac{\sum_{j=1}^n C_j^2}{n} - (\overline{c})^2}

.. math::
    U=\sqrt{\frac{[n\sum_{j=1}^n W_{i,j}^2 - (\sum_{j=1}^n W_{i,j})^2]}{n-1}}

.. math::
    W_{i,j} = \left\{\begin{equation}\begin{array}{lr}
                 \text{1 if j is a neighbor of i}\\
                 \text{0 if j is not a neighbor of i}
                 \end{array}
    \end{equation}\right.

:math:`n` total number of grid regions

:math:`C_j` Cell count for region j

:math:`\overline{C}` mean of cell count in all region

`A more illustrative example <https://www.nature.com/articles/modpathol201537>`_

Communities detection
----------------------

This is used to find communities in a ROI.
The neighbors relationships are convert to graph, each cell is a node, two nodes are connected if
they are neighbors, edge weight is represented by distance. Using leidenalg algorithm, we can detect
the communities within a ROI.


