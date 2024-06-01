import csv
import pytest
import rainfall_functions

@pytest.fixture
def data():
    filename = "rainfall_in_india_1901-2015.csv"

    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        data_samp = next(reader)
        return data_samp

def test_annual_mean(data):

    result = rainfall_functions.calculate_annual_mean(data)
    assert(result == 281.1)

def test_subdivisions_dict(data):
    result = rainfall_functions.create_subdivisions_dict([data])  
    assert(result == {'ANDAMAN & NICOBAR ISLANDS': {'years': ['1901'], 'means': [281.1]}})
