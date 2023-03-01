import pandas as pd

# Load pressure data from CSV file. 
csv = pd.read_csv('slice_h_10m.csv')


# Only keep spatial and pressure data. 
csv = csv[['Points_0', 'Points_1', 'Points_2','p_Ar']]

# Drop x-coordinates.
csv = csv.drop(['Points_0'], axis = 1)

# Convert coordinates from [m] to [mm]
csv['Points_1'] = 1000 * csv['Points_1'].convert_dtypes(convert_integer = True)
csv['Points_2'] = 1000 * csv['Points_2'].convert_dtypes(convert_integer = True)

# Print some statistics
print('average values')
print(csv.mean(axis = 0))
print()
print('max values')
print(csv.max(axis = 0))
print()
print('mininum values')
print(csv.min(axis = 0))
print()
# Save data to a new csv file. 
csv.to_csv('pressure_distribution_h_10m.csv')

print(csv.columns)
print(csv)