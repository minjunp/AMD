# Inferring genomic signature in Age Macular Degeneration (AMD)
We developed machine learning models (Multinomial Lasso regression, XGboost, Random Forest) to predict distinct stages of age-related macular degeneration (AMD).
You can find our final report:
```
Final-Report.pdf
```
and proposal for this project:
```
Proposal.pdf
```

### Due to the massive size of the data, datasets needed to run the code are available at:
https://rice.box.com/s/lvuwn3g8e4az5kua2y1qz3drmkgnhq9t
Please download the data in folder: D2K_BCM_DATASET

### You can import the below packages using the conda command:
`conda env create -f environment.yml`

This is our Envrionment and Packages that we need

packages:  
pandas==0.24.2   
seaborn==0.9.0  
umap_learn==0.3.10  
matplotlib==3.1.0          
patsy==0.5.1      
numpy==1.18.1  
scipy==1.4.1    
scikit_learn==0.22.2.post1
py-xgboost-mutex=2.0
