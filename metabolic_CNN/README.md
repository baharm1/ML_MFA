# Metabolic CNN: Quantification of relative fluxes of metabolite sources


## Requirements
We recommend installing [Anaconda](https://www.anaconda.com/download) and create a new conda environment to manage the versioning of dependencies:
```
conda create -n mCNN python=3.8.15
conda activate mCNN
```
To run Python codes, the following packages have to be installed:
```
pip install -r requirements.txt
```


Add mouse and patient data 
Run CNN_pred_mouse_patient

R version 4.2.2

ggplot2
3.4.2

ggbeeswarm
0.7.2

ggpubr
0.4.0

Seurat
4.2.0

tidyr
1.3.0

dplyr
1.1.2

stringr
1.5.0



Python 3.8.15

GPU-enabled workstation with NVIDIA Quadro P2200

matplotlib
3.7.2

optuna
3.1.0

plotly
5.15.0

joblib
1.3.0

scipy
1.10.1

scikit-learn
1.3.0

tensorflow
'2.10.1'

tensorflow.keras
'2.10.0'

pandas
'2.0.3'

numpy
'1.24.4'

installing tensorflow on windows:

8.1.0.77
11.2.2
conda install -c conda-forge cudatoolkit=11.2 cudnn=8.1.0
# Anything above 2.10 is not supported on the GPU on Windows Native
python -m pip install "tensorflow<2.11"
# Verify the installation:
python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"