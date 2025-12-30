EVT Climate Analysis Tool

EVT Climate Analysis Tool is a professional, open-source Python package for analyzing extreme values in climate and environmental datasets.
It applies the principles of Extreme Value Theory (EVT) to identify, model, and report rare but high-impact events — such as floods, heatwaves, and extreme precipitation.

This tool automates the entire workflow — from data import to generating publication-ready PDF reports with figures, statistical tables, and trend analyses.

🧩 Key Features

📈 Annual Block Maxima (GEV): Fits Generalized Extreme Value distribution to yearly maxima.

🎚 Peaks Over Threshold (POT): Detects and models exceedances above a chosen threshold.

🔁 Rolling POT analysis: Computes time-varying shape parameter ξ(t) using a moving window.

🎲 Bootstrap uncertainty estimation: Confidence intervals via repeated resampling.

📉 Non-stationary σ(t) modeling: Detects trends in scale parameter over time.

📊 Return levels estimation: Calculates expected magnitudes of 50-year and 100-year extremes.

🧠 Trend and change-point tests: Mann–Kendall (trend) and Pettitt (abrupt change) tests.

📑 Automated PDF reporting: Generates a complete, publication-ready report with all results and figures.


📊 Input Data Format

The input dataset should be a CSV file containing daily observations with date columns (year, month, day) and at least one numeric variable to analyze.

Example:

year,month,day,value
2000,1,1,12.4
2000,1,2,13.2
2000,1,3,10.8
...
2020,12,31,14.5


Each numeric column (like “value”, “precipitation”, or “temperature”) will be analyzed separately.

✅ You can include multiple numeric columns — the script automatically detects and processes each one.

⚙️ How to Use
1️⃣ Install Dependencies

Make sure you have Python 3.8+ installed, then install required packages:

pip install -r requirements.txt

2️⃣ Prepare Your Data

Place your input CSV file (e.g., btmax.csv) inside the sample_data/ folder.
Ensure it has year, month, and day columns and one or more numeric variables.

3️⃣ Run the Analysis

Execute the script from the project root:

python evt_analysis.py


The program will:

Read and process your data,

Fit GEV and POT models,

Estimate uncertainty and return levels,

Perform Mann–Kendall and Pettitt tests,

Generate plots and tables,

Export a comprehensive PDF report.

4️⃣ View Results

After execution, the final report will be saved automatically at:

output/EVT_Report.pdf


The PDF includes:

Fitted GEV histogram

Time series of annual maxima

POT threshold stability plots

Rolling ξ(t) trend visualization

Return levels table

Statistical interpretation text

📘 Example Output

File: EVT_Report.pdf

Contents:

Title page with analyzed variable name

Tables for GEV parameters (μ, σ, ξ)

Graphs for annual maxima and fitted curves

POT threshold sensitivity analysis

Rolling shape parameter ξ(t) plot with bootstrap confidence intervals

Trend test results and return level estimates

This report is publication-ready, suitable for scientific papers, theses, or climate risk assessments.

🧮 Example Use Cases

Climate extremes (rainfall, temperature, wind speed)

Hydrological modeling (river discharge peaks, floods)

Environmental risk assessment

Infrastructure resilience studies

Insurance and catastrophe modeling

📦 Requirements
Library	Purpose
numpy	Numerical operations
pandas	Data handling and time series management
matplotlib	Plotting and visualization
scipy	Statistical modeling (GEV, GPD)
reportlab	PDF report generation

Install all dependencies via:

pip install -r requirements.txt

🧠 Methodology Summary

This tool implements key components of Extreme Value Theory:

Method	Description
GEV (Generalized Extreme Value)	Models block maxima (e.g., yearly max precipitation).
POT (Peaks Over Threshold)	Uses all data exceeding a defined threshold.
Bootstrap	Resamples data to estimate uncertainty in ξ (shape parameter).
Rolling window	Tracks temporal change in extreme behavior.
Mann–Kendall	Non-parametric test for monotonic trends.
Pettitt test	Detects change-points in time series.
📈 Return Levels

Return levels (e.g., 50-year, 100-year) represent the magnitude of an event expected once every T years, on average.

This helps planners and engineers assess design safety margins for rare, extreme events.

🧩 Customization

In the script (evt_analysis.py), you can modify:

FILE_PATH = "sample_data/btmax.csv"
OUTPUT_DIR = "output"
POT_QUANTILES = [0.9, 0.95, 0.99, 0.995]
ROLLING_WINDOW = 5
RETURN_YEARS = [50, 100]


These parameters control:

Input file location

Output folder

POT quantiles tested

Rolling window size

Return periods for calculation

🧑‍💻 Example Terminal Output
✅ PDF generated successfully:
output/EVT_Report.pdf

🧾 License

This project is released under the MIT License — you’re free to use, modify, and distribute it for academic or professional purposes.

👨‍🔬 Citation (for research use)

If you use this tool in your research, please cite:

Dashtbozorgi,A; (2025). EVT Climate Analysis Tool: An open-source Python framework for modeling extreme climate events.
