import importlib.util
import pathlib
import pandas as pd
from unittest import mock

# Load module from file path with patched pandas and print to avoid heavy output
module_name = 'Fitness_Athletics'
module_path = pathlib.Path(__file__).resolve().parents[1] / f'{module_name}.py'
spec = importlib.util.spec_from_file_location(module_name, module_path)

with mock.patch('pandas.read_csv') as read_csv, mock.patch('builtins.print'):
    read_csv.return_value = pd.DataFrame({'rsid': [], 'chromosome': [], 'position': [], 'genotype': []})
    Fitness_Athletics = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(Fitness_Athletics)

calculate_probability = Fitness_Athletics.calculate_probability


def test_calculate_probability_zero():
    assert calculate_probability(0) == 50.0


def test_calculate_probability_negative_one():
    assert calculate_probability(-1) == 0.0


def test_calculate_probability_one():
    assert calculate_probability(1) == 100.0
