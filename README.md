# SIMCCDA

## Title:

Prediction of circRNA-disease associations based on inductive matrix completion.

## Developers: 

Menglu Li (mengluli@foxmail.com) and Mengya Liu ([woshiliangming@outlook.com](mailto:woshiliangming@outlook.com) ) from Institutes of Physical Science and Information Technology and School of Computer Science and Technology, Anhui University.



## Related Files

### Dataset:

```
========================================================================================
| FILE NAME            | DESCRIPTION                                                   |
========================================================================================
|chr_diseasematrix.csv | known circRNA-disease associations information, 1 indicates that circRNA and disease are associated, 0 indicates that their associations is currently in an unknown state.                                                         |
|dissimilarity.csv     | disease semantic similarity (DOSim (DOSE R package)).         |
|seqsimilarity.csv     | circRNA sequence similarity (Levenshtein distance).           |
```

### Code:

```
========================================================================================
| FILE NAME       | DESCRIPTION                                                        |
========================================================================================
| main_LOOCV.m    | leave-one-out cross validation.                                    |
| SIMCCDA_demo.m  | function predicting potential circRNA-disease associations.        |
| gkl.m           | function computing Gaussian interaction profile kernel.            |
| PCA.m           | function extracting primary feature vectors via PCA.               |
| SIMC            | inductive matrix completion.                                       |
| roc.m           | function computing area under ROC curve.                           |

```



## **Contact**

Please feel free to contact us if you need any help: mengluli@foxmail.com