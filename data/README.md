# Glucose data set

Section 7.3 of the paper uses the public "Diabetes prediction dataset" by Mohammed Mustafa:
https://www.kaggle.com/datasets/iammustafatz/diabetes-prediction-dataset

The file `diabetes_prediction_dataset.csv` (100,000 records, 9 columns, 3,810,356 bytes) is not redistributed here.
Download it from Kaggle, place it in this folder and run

    cd ../code
    python kaggle_summary.py ../data/diabetes_prediction_dataset.csv

The script prints the SHA-256 checksum of the file. The copy used for the paper has

    f233fc470e904201f35310d061cb46ca0a6bc7ca5c3e428b769fb4f87f69929b

Only the column `blood_glucose_level` (mg/dL) is used, together with the `diabetes` label for display in Fig. 12(a).
The aggregate that enters the paper, the 18 distinct glucose levels and the number of records at each level, is stored
in `../out/kaggle_levels.json`. The data set contains no optical measurements; it supplies test concentrations and
their weights for the accuracy test against ISO 15197:2013.
