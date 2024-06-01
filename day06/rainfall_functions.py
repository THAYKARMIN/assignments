
def calculate_annual_mean(row):
    try:
        monthly_rainfall = [float(row[month]) for month in ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']]
        print(monthly_rainfall)
    except ValueError:
        return None  # Return None if there is invalid data

    annual_mean = sum(monthly_rainfall) / len(monthly_rainfall)
    print(annual_mean)
    return annual_mean


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