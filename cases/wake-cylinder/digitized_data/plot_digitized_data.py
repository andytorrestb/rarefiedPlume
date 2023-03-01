import matplotlib.pyplot as plt
import pandas as pd

pressure_ratio = pd.read_csv('pressure_ratio.csv')
# print(pressure_ratio)
print(pressure_ratio.columns)

# Break up curve data in the csv file into seperate dataframes according each plotted case.
curves = {}
for index, column in enumerate(pressure_ratio.columns):
    # Find every column that contains x data. 
    if 'x' in column:
        # Group columns of data acording to the current and proceeding index.
        # print(index, column)
        # print(pressure_ratio.columns[index], pressure_ratio.columns[index+1])
        x_data = pressure_ratio[pressure_ratio.columns[index]]
        p_ratio = pressure_ratio[pressure_ratio.columns[index+1]]
        curve = [x_data, p_ratio]
        curves[pressure_ratio.columns[index+1]] = curve

# Plot data on the same plot. 
import matplotlib
print(matplotlib.get_configdir())
input()
print(plt.style.available)
# plt.style.use(['science', 'ieee'])
plt.style.use('seaborn-paper')


for curve in curves:
    x = curves[curve][0]
    y = curves[curve][1]
    plt.plot(x,y, label = curve)
    print(type(curves[curve]))

font_dict  = {"fontname":"Times New Roman", "color":"black", "size":20}

plt.yscale('log')
plt.xscale('log')
plt.legend(title='Plate Distance')
plt.ylabel('Pressure ratio (torr)')
plt.xlabel('Reservoir Pressure (torr)')
plt.title('Cylinder Windward to Leeward Pressure Ratio')
plt.savefig('digitized_pressure_ratio.png', dpi=800)
