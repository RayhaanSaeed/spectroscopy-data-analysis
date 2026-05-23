# spectroscopy-data-analysis
Python-based spectroscopy database and analysis pipeline for processing, standardising, and analysing fluorescence spectra of photofluorescent proteins.

# Spectroscopy Database and Analysis Pipeline

## Overview

This project focuses on the collection, processing, organisation, and analysis of spectroscopy data for photofluorescent proteins. The aim of the project is to create a structured and reusable database of excitation and emission spectra collected from journal papers, repositories, and experimental datasets.

Using Python, the project automates the cleaning, restructuring, standardisation, and export of spectroscopy datasets to allow meaningful comparison between different fluorescent proteins and their spectral behaviour.

---

## Project Structure

The project currently consists of two main Python workflows:

### 1. Data Import and Cleaning Pipeline

This script collates spectroscopy datasets from multiple unstructured sources, including CSV files, repositories, and experimental data exports. The workflow cleans and restructures inconsistent data formats, removes duplicates, handles missing information, and standardises the output into organised CSV datasets suitable for analysis.

Key tasks include:
- Importing spectroscopy datasets from multiple sources
- Cleaning inconsistent formatting and messy data
- Removing duplicate or incomplete entries
- Standardising wavelength and intensity data
- Generating structured CSV outputs for downstream analysis

---

### 2. Data Export and Analysis Preparation Pipeline

This script exports selected spectroscopy datasets for groups of proteins into OriginPro for graphical comparison and further data analysis. The workflow allows spectral datasets to be compared visually to identify trends and differences in fluorescence behaviour between proteins.

Key tasks include:
- Selecting datasets for groups of proteins
- Exporting processed spectroscopy data into analysis-ready formats
- Preparing datasets for OriginPro visualisation
- Supporting graphical comparison of excitation and emission spectra
- Enabling trend analysis and spectral interpretation

---

## Objectives

- Collect spectroscopy datasets from scientific literature and repositories
- Standardise inconsistent and unstructured spectroscopy data
- Automate data cleaning and restructuring workflows
- Extract meaningful spectral information from raw datasets
- Enable comparison and visualisation of fluorescence behaviour between proteins
- Create reusable datasets for future computational and statistical analysis

---

## Techniques Used

### Programming and Libraries
- Python
- pandas
- NumPy
- matplotlib
- SciPy
- pyodbc

### Data Science Techniques
- Data cleaning and preprocessing
- Handling unstructured datasets
- Data standardisation
- Automated data pipelines
- Feature extraction
- Scientific data visualisation
- Statistical analysis

---

## Example Features Analysed

- Peak excitation wavelength
- Peak emission wavelength
- Relative fluorescence intensity
- Spectral trends between proteins
- Comparative fluorescence behaviour

---

## Skills Demonstrated

- Scientific programming
- Data pipeline development
- Handling large unstructured datasets
- Data cleaning and preprocessing
- Computational analysis
- Data visualisation
- Workflow automation
- Independent project management

---

## Future Improvements

- Integration of machine learning models for spectral prediction
- Automated feature classification
- Interactive data visualisation dashboards
- Expansion of dataset coverage and metadata
- Expanded dataset coverage
