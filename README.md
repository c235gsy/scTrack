

# scTrack
 scTrack is a single-cell trajectory inference tool developed by [Jin_lab@SUSTech](http://jinwlab.net/).

### Before using

#### Modules scTrack depend on 

> 1. Python3
> 2. numpy
> 3. pandas
> 4. sk-learn
> 5. matplotblib
> 6. seaborn
> 7. colour
> 8. rpy2
> 9. R

#### Install those Modules

[Anaconda](https://www.anaconda.com/products/individual) is recommended to used to install the modules which the scTrack depends on. 

### Test

After successfully installing anaconda, you can run the following code.

```
conda install scikit-learn seaborn rpy2 matplotblib r-base numpy pandas colour 
```

After cloning this repository, you can run the following code to test. (Before running code, you need to decompressing the file **test.zip**)

```
python get_paths.py -e ./test/exp_matrix.txt -c ./test/ident.txt -o ./test/output/Test -b "[1]" -d "[3,4]" -ul -pbf  
python analyze_paths.py -e ./test/exp_matrix.txt -c ./test/ident.txt -o ./test/output/Test -i ./test/cell_id.txt -p ./test/output/Test_paths.txt -p2D -p3D -cd -cb
```

### Description of the parameters

-------

#### get_paths.py

Firstly, run **get_paths.py** to get the general structure of the cell lineage, and the user will get **a path file**. The output will contain **xxx_paths.txt, xxx_all_paths.txt, xxx_BF.txt, xxx_BF.pdf.**

* xxx_paths.txt
    - **_This file contains the filtered inferred paths._** The format of lines in **xxx_paths.txt** is **[path_name, path_score, cluster1_in_path, cluster2_in_path, cluster3_in_pat, ...... ]**.
* xxx_all_paths.txt
    - **_This file contains all paths inferred._** The format of lines in **xxx_paths.txt** is **[path_name, path_score, cluster1_in_path, cluster2_in_path, cluster3_in_pat, ...... ]**.
* xxx_BF.txt
    - **_This file contains the Bayes Factor Matrix between cell clusters._**
* xxx_BF.pdf
    - Only output if the user use **--plotBayesFactorMatrix(-pbf)**
    - **_This is the heatmap of the Bayes Factor Matrix between cell clusters._**

-------
> **--expMatrix(-e) :**  *file path* **(required)**

> File path of expression matrix.

-------
> **--clusterArray(-c) :**  *file path* **(required)**

> File path of cell cluster marker array, _the number of clusters should not be less than 3._

-------
> **--output(-o) :**  *string* 

> Paths and names of the output files, **default='./YOUR_OUTPUT_NAME'**.

-------
> **--nComponentsPCA(-np) :**  *integer*

> The data dimension retained by PCA dimensionality reduction. Note that this parameter should be greater than **[cell_cluster_number-1].** By default it is equal to **max[5*[[cell_cluster_number]-1], int[0.01*[gene_number]].** 

-------
> **--nComponentsLDA(-nl) :**  *integer*

> the data dimension retained by LDA dimensionality reduction. Note that this parameter should not be greater than **[cell_cluster_number]-1**. By default it is equal to **min[[cell_cluster_number]-1, --nComponentsPCA]**.

-------
> **--sizeSubcluster(-s) :**  *float*

>  The approximate average size of the subclusters, **default=50**.

-------
> **--kNeighbor(-k) :**  *integer*

> A cell cluster will establish connections with the surrounding **-k** clusters, **default=3**.

-------
> **--nCalculationBayesFactor(-nbf) :**  *integer*

> The number of cell pairs in calculation to get the distribution of distances between cells from 2 cell clusters, **default=10^5**.

-------
> **--beginning(-b) :**  *string* **(required)**

> The format is **"[clusterA,clusterB,...]"**. Select one or more clusters as the beginning.

-------
> **--destination(-d) :**  *string*

> The format is **"[clusterA,clusterB,...]"**. Select any number of clusters as the destination. (type **"[]"** if there is no destination), **default="[]"**.

-------
> **--newDestination(-nd) :**  *no input*

> Only work while **--destination != []**. If **--destination== []**, it will be used automaticallyIf it is used, the program will output the paths whose end points are not contained in **--destination**. Otherwise, those paths will be ignored.

-------
> **--mandatoryLink(-ml) :**  *string*

> The format is **"{clusterA:[clusterB,clusterC,...], clusterL:[clusterM,clusterN,...],...}"**. Default value is **"{}"**. The command in the example means to force the establishment of links from **clusterA** to **clusterB** and **clusterC**, and links from **clusterL** to **clusterM** and **clusterN**.

-------
> **--onlyMandatoryLink(-oml) :**  *no input*

> If used, the links form clusters included in **[--mandatoryLink]** to those not included in **[--mandatoryLink]** will be canceled.

-------
> **--cancelLink(-cl) :**  *string*

> The format is **"{clusterA:[clusterB,clusterC,...],clusterL:[clusterM,clusterN,...],...}"**. Default value is **"{}"** The command in the example means to force the cancel of links from **clusterA** to **clusterB** and **clusterC**, and links from **clusterL** to **clusterM** and **clusterN**.

-------
> **--undirectedLink(-ul) :**  *no input*

> If used, the links between clusters will be transformed into bidirectional links.

-------
> **--plotBayesFactorMatrix(-pbf) :**  *no input*

> If used, plot the heat map of BF values.

-------

#### analyze_paths.py

You can view and modify the **xxx_paths.txt** or **xxx_all_paths.txt**, and run the second step of the program **analyze_paths.py** to complete the trajectory inference, where **xxx_paths.txt** or **xxx_all_paths.txt** is one of the input files. 

-------
> **--expMatrix(-e) :**  *file path* **(required)**

> File path of expression matrix.

-------
> **--clusterArray(-c) :**  *file path* **(required)**

> File path of cell cluster marker array, _the number of clusters should not be less than 3._

-------
> **--output(-o) :**  *string* 

> Paths and names of the output files, **default='./YOUR_OUTPUT_NAME'**.

-------
> **--idArray(-i) :**  *file path* **(required)**

> File path of cell-id array. If missing, cell-id defaults to **Cell1, Cell2, Cell3, ...**

-------
> **--nComponentsPCA(-np) :**  *integer*

> The data dimension retained by PCA dimensionality reduction. Note that this parameter should be greater than **[cell_cluster_number-1].** By default it is equal to **max[5*[[cell_cluster_number]-1], int[0.01*[gene_number]].**

-------
> **--nComponentsLDA(-nl) :**  *integer*

> the data dimension retained by LDA dimensionality reduction. Note that this parameter should not be greater than **[cell_cluster_number]-1**. By default it is equal to **min[[cell_cluster_number]-1, --nComponentsPCA]**.

-------
> **--maxIter(-mi) :**  *integer*

> The maximum number of iterations in the principal curve fitting process, **default=10**.
 
 -------
> **--cutoffIter(-ci) :**  *float*

> The iterations in the principal curve fitting process will stop if the rate of change between the point and its corresponding point on the curve is less than this value, **default=0.001**.
 
 -------
> **--pathFile(-p) :**  *file path* **(required)**

> The path of the path-file inferred in the previous step.

 -------
> **--combineBeginnings(-cb) :**  *no input*

>  If used, combining the shared parts of paths in the beginning.

 -------
> **--combineDestinations(-cd) :**  *no input*

> If used, combining the shared parts of paths in the destination.

-------
> **--plot2D(-p2D) :**  *no input*

>  If used, plotting 2D data figures.
 
 -------
> **--plot3D(-p3D) :**  *no input*

>  If used, plotting  a lot 3D data figures.
 
 -------
> **--colorMap(-cm) :**  *string* **(Only can be chosen form {R,r,Python,py})**

> If **R** or **r**, you also need to install **Python** package **rpy2** in order touse the function **hue_pal()** in **R** to get colormap. But the **hue_pal()** do not use standard **HUE space**. If **Python** or **py**, program will use **Python** package **Colour**, where the standard **HUE space** is used, **default=R**.

 -------
> **--smoothness(-s) :**  *float*

> Positive smoothing factor used to choose the number of knots. **-s = smoothness\*[point_number_in_curve]**, **default=1.0**. More information: [https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html)

### Author
Siyuan Guo 11611118@mail.sustech.edu.cn





