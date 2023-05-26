# mrMediator

Quick review on mediator MR.

Please have a look at the following github page [https://github.com/CESP-ExpHer/MVMR/blob/main/README.md](https://github.com/CESP-ExpHer/MVMR) in order to understand the concept of Mediator MR is and how to prepare the dataset correctly. 

**Make sure the column names correspond exactly as speicified below**
```
head(DTC)
  chromosome  position EA OA       EAF    Beta          SE             P Phenotype        SNP
1          5  76515824  A  G 0.5714411 -0.1475 0.004671074 4.794161e-219       TSH  rs9687206
2          6 166043959  C  G 0.3071071 -0.1155 0.004975709 1.386952e-120       TSH  rs2983511
3          4 149665602  T  C 0.8041840  0.1176 0.005788069  6.765580e-93       TSH rs11732089
4          6  43806609  A  G 0.6938409  0.0992 0.005077254  1.525257e-86       TSH   rs881858
5          2 217628430  A  G 0.2792577 -0.0970 0.005178799  4.109376e-79       TSH  rs1861628
6         16  79747855  A  G 0.6673459  0.0905 0.004874164  2.616128e-78       TSH rs73575083
```

The first step is to calculate both the direct and indirect effect using **mrMediator** function.

```
fit = mrMediator(data=DTC,
                 exposure = list(name='TSH', type = 'continuous'),
                 mediator = list(name ='SHBG', type = 'continuous'),
                 outcome = list(name='DTC', type ='binary'))
```
Once fitted, we can draw the mediation diagram. We can either put **method = 'prod'** to have indirect effect using product method
```
mrDiagram(fit =fit, method = 'prod')
```
![mr-prod](https://github.com/CESP-ExpHer/mrMediator/assets/24691084/b4d6e560-933a-4151-9df1-9d314ba0f3a2)


or **method = 'diff'** to have indirect effect using difference method
```
mrDiagram(fit =fit, method = 'diff')
```
![mr-diff](https://github.com/CESP-ExpHer/mrMediator/assets/24691084/b654230d-620f-469b-8ef2-b6e467c65d32)
