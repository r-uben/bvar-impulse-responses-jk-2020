from bvar_impulse_responses_jk_2020.aws_manager.bucket_manager import BucketManager
from bvar_impulse_responses_jk_2020.date_format import DateFormatUtils

import calendar
import pandas as pd

def get_gdp_data(bucket: BucketManager):
    """
    Retrieve GDP data from S3 bucket and format the date column.
    
    Args:
        bucket: BucketManager instance to access S3 data
        
    Returns:
        DataFrame containing GDP data with properly formatted date column
    """
    gdp_data = bucket.read_document("usa/national_accounting/gdp.csv", format="csv")
    gdp_data.columns = [x.lower() for x in gdp_data.columns]
    
    if 'date' not in gdp_data.columns:
        gdp_data = DateFormatUtils.add_date_column(gdp_data, frequency='monthly')
    return gdp_data


def get_sp500():
    sp500_data = pd.read_csv("data/fred/SP500/data.csv")
    sp500_data.columns = [x.lower() for x in sp500_data.columns]
    sp500_data = DateFormatUtils.add_date_column(sp500_data, frequency='daily')
    sp500_data.columns = ['sp500']
    return sp500_data

def fed_funds():
    fed_funds_data = pd.read_csv("data/fred/FEDFUNDS/data.csv")
    fed_funds_data.columns = [x.lower() for x in fed_funds_data.columns]
    fed_funds_data = DateFormatUtils.add_date_column(fed_funds_data, frequency='daily')
    fed_funds_data = fed_funds_data.resample('ME').last()
    fed_funds_data.columns = ['fed_funds']
    return fed_funds_data



def main():
    bucket = BucketManager("macroeconomic-data")
    gdp_data = get_gdp_data(bucket)
    sp500_data = get_sp500()
    fed_funds_data = fed_funds()
    print(gdp_data)
    print(sp500_data)
    print(fed_funds_data)




if __name__ == "__main__":
    main()
