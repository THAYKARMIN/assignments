#%%
import csv
import matplotlib.pyplot as plt
#%%
#%%
def calculate_annual_mean(row):
    try:
        monthly_rainfall = [float(row[month]) for month in ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']]
    except ValueError:
        return None  # Return None if there is invalid data

    annual_mean = sum(monthly_rainfall) / len(monthly_rainfall)
    return annual_mean

#%%
def create_subdivisions_dict(data):
    subdivisions = {}
    for row in data:
        subdivision = row['SUBDIVISION']
        year = row['YEAR']
        
        annual_mean = calculate_annual_mean(row)
        if annual_mean is None:
            continue  # Skip rows with invalid data

        if subdivision not in subdivisions:
            subdivisions[subdivision] = {'years': [], 'means': []}
        subdivisions[subdivision]['years'].append(year)
        subdivisions[subdivision]['means'].append(annual_mean)
        
    return subdivisions

#%%
def plot_subdivision_means(subdivisions):
    for subdivision, data in subdivisions.items():
        years = data['years']
        means = data['means']
        plt.figure(figsize=(20, 6))  
        plt.plot(years, means, marker='o')
        plt.title(f'Mean Annual Rainfall for {subdivision}')
        plt.xlabel('Year')
        plt.ylabel('Mean Annual Rainfall (mm)')
        plt.grid(True)
        plt.xticks(rotation=90)  
        plt.tight_layout()
        plt.xlim(min(years), max(years))
        plt.show()         
#%%

def main():
    filename = "rainfall_in_india_1901-2015.csv"
    data = []
    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            data.append(row)

    subdivisions = create_subdivisions_dict(data)
    plot_subdivision_means(subdivisions)


main()
#%%