import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
import shap
import matplotlib.pyplot as plt
from openpyxl import Workbook

# 1. Reading data from an Excel file with multiple sheets
file_path = "G:\\variation_and_attribution_for_runoff_and_sediment_flux\\codes2\\input_data\\relative_changes_in_pre_ta_ndvi_es_q_for_random_forest.xlsx"
excel_file = pd.ExcelFile(file_path)
sheet_names = excel_file.sheet_names

# Create a dictionary to store results for each sheet
results = {}

for sheet_name in sheet_names:
    df = excel_file.parse(sheet_name)

    # Extract predictor variables (X1, X2, X3) and response variable (Y)
    X = df[['pre', 'ta', 'ndvi','es']]
    Y = df['q']

    # Split the data into training and testing sets
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3, random_state=42)

    # 2. Building a random forest model for each sheet
    rf = RandomForestRegressor(n_estimators=100, random_state=42)
    rf.fit(X_train, Y_train)

    Y_train_fitted=rf.predict(X_train)
    Y_test_fitted=rf.predict(X_test)
    # 3. Calculating SHAP values
    explainer = shap.TreeExplainer(rf)
    shap_values_train = explainer.shap_values(X_train)
    shap_values_test = explainer.shap_values(X_test)

    shap_values=np.append(shap_values_train,shap_values_test,axis=0)
    # 4. Exporting results to Excel files
    # Create a workbook for this sheet
    workbook = Workbook()
    shap_sheet = workbook.active
    shap_sheet.title = 'SHAP Values'

    # Add SHAP values to the sheet
    for row in shap_values:
        shap_sheet.append(list(row))

    # Create a sheet for fitted and test values
    values_sheet = workbook.create_sheet(title='Fitted and Test Values')
    values_sheet.append(['q_train_obs', 'q_train_fitted','q_test_obs','q_test_fitted'])
 
    print(Y_train)
    Y_train_obs_fitted=np.transpose([Y_train,Y_train_fitted])
    Y_test_obs_fitted=np.transpose([Y_test,Y_test_fitted])

        
    Y_test_obs_fitted1=np.full((len(Y_test_obs_fitted),len(Y_test_obs_fitted[0])+2),None,dtype=object)
    Y_test_obs_fitted1[:,2:4]=Y_test_obs_fitted
  
    for row in Y_train_obs_fitted:
        values_sheet.append(list(row))
  

    for row in Y_test_obs_fitted1:
        
        values_sheet.append(list(row))

    # Save the workbook to a file
    workbook.save(f"G:\\variation_and_attribution_for_runoff_and_sediment_flux\\codes2\\output\\{sheet_name}_results.xlsx")

    # Store SHAP values for this sheet in the results dictionary
    

# 5. Plotting bar graphs for SHAP values of all sheets
#plt.figure(figsize=(15, 10))
#for sheet_name, shap_valuesx in results.items():
#    plt.bar(range(len(shap_values)), np.mean(np.abs(shap_values), axis=0), label=sheet_name)
#plt.xlabel('Features')
#plt.ylabel('Mean SHAP Value Magnitude')
#plt.xticks(range(len(X.columns)))
##plt.legend()
##plt.title('SHAP Values for Each Sheet')
#plt.show()

