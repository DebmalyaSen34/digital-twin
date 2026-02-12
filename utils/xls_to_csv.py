import pandas as pd

def xls_to_csv(xls_file, csv_file):
    # Read the Excel file
    df = pd.read_excel(xls_file, sheet_name=0)
    
    # Save as CSV
    df.to_csv(csv_file, index=False)
    
    print(f"Converted '{xls_file}' to '{csv_file}' successfully.")

if __name__ == "__main__":
    xls_file = "Table_S2.xls"
    csv_file = "biomass.csv"
    xls_to_csv(xls_file, csv_file)
    print(pd.read_csv(csv_file).head())