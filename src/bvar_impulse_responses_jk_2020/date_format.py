import pandas as pd
import calendar

class DateFormatUtils:

    @staticmethod
    def add_date_column(df, frequency='monthly'):
        """
        Add a date column to the dataframe based on time columns and frequency.
        If date column already exists, return DataFrame with date as index.
        
        Args:
            df: DataFrame containing time columns (year/month/day or year/quarter)
                or DataFrame with existing date column
            frequency: String indicating data frequency ('monthly', 'daily', or 'quarterly')
            
        Returns:
            DataFrame with 'date' as index in YYYY-MM-DD format
        """
        df = df.copy()
        
        if 'date' in df.columns:
            df = df.set_index('date')
            df.index = pd.to_datetime(df.index)
            return df
            
        if frequency == 'monthly':
            year_col, month_col = DateFormatUtils._get_year_month_columns(df)
            df = DateFormatUtils._convert_year_month_to_int(df, year_col, month_col)
            DateFormatUtils._validate_month_values(df, month_col)
            df = DateFormatUtils._add_date_with_last_day(df, year_col, month_col)
        
        elif frequency == 'daily':
            year_col, month_col, day_col = DateFormatUtils._get_year_month_day_columns(df)
            df = DateFormatUtils._convert_year_month_day_to_int(df, year_col, month_col, day_col)
            DateFormatUtils._validate_month_values(df, month_col)
            DateFormatUtils._validate_day_values(df, year_col, month_col, day_col)
            df = DateFormatUtils._add_date_with_day(df, year_col, month_col, day_col)
            
        elif frequency == 'quarterly':
            year_col, quarter_col = DateFormatUtils._get_year_quarter_columns(df)
            df = DateFormatUtils._convert_year_quarter_to_int(df, year_col, quarter_col)
            DateFormatUtils._validate_quarter_values(df, quarter_col)
            df = DateFormatUtils._add_date_for_quarter(df, year_col, quarter_col)
            
        else:
            raise ValueError("Frequency must be 'monthly', 'daily', or 'quarterly'")
            
        df.index = pd.to_datetime(df.index)
        return df

    @staticmethod
    def _get_year_month_columns(df):
        """Get year and month column names from dataframe."""
        year_col = next((col for col in df.columns if col.lower() == 'year'), None)
        month_col = next((col for col in df.columns if col.lower() == 'month'), None)
        
        if not year_col or not month_col:
            raise ValueError("DataFrame must contain 'year' and 'month' columns (case-insensitive)")
            
        return year_col, month_col

    @staticmethod
    def _get_year_month_day_columns(df):
        """Get year, month and day column names from dataframe."""
        year_col = next((col for col in df.columns if col.lower() == 'year'), None)
        month_col = next((col for col in df.columns if col.lower() == 'month'), None)
        day_col = next((col for col in df.columns if col.lower() == 'day'), None)
        
        if not all([year_col, month_col, day_col]):
            raise ValueError("DataFrame must contain 'year', 'month' and 'day' columns (case-insensitive)")
            
        return year_col, month_col, day_col

    @staticmethod
    def _get_year_quarter_columns(df):
        """Get year and quarter column names from dataframe."""
        year_col = next((col for col in df.columns if col.lower() == 'year'), None)
        quarter_col = next((col for col in df.columns if col.lower() == 'quarter'), None)
        
        if not year_col or not quarter_col:
            raise ValueError("DataFrame must contain 'year' and 'quarter' columns (case-insensitive)")
            
        return year_col, quarter_col

    @staticmethod
    def _convert_year_month_to_int(df, year_col, month_col):
        """Convert year and month columns to integers."""
        df[year_col] = pd.to_numeric(df[year_col], errors='coerce').fillna(2000).astype(int)
        df[month_col] = pd.to_numeric(df[month_col], errors='coerce').fillna(1).astype(int)
        return df

    @staticmethod
    def _convert_year_month_day_to_int(df, year_col, month_col, day_col):
        """Convert year, month and day columns to integers."""
        df[year_col] = pd.to_numeric(df[year_col], errors='coerce').fillna(2000).astype(int)
        df[month_col] = pd.to_numeric(df[month_col], errors='coerce').fillna(1).astype(int)
        df[day_col] = pd.to_numeric(df[day_col], errors='coerce').fillna(1).astype(int)
        return df

    @staticmethod
    def _convert_year_quarter_to_int(df, year_col, quarter_col):
        """Convert year and quarter columns to integers."""
        df[year_col] = pd.to_numeric(df[year_col], errors='coerce').fillna(2000).astype(int)
        df[quarter_col] = pd.to_numeric(df[quarter_col], errors='coerce').fillna(1).astype(int)
        return df

    @staticmethod
    def _validate_month_values(df, month_col):
        """Validate that month values are between 1 and 12."""
        if (df[month_col] < 1).any() or (df[month_col] > 12).any():
            raise ValueError("Month values must be between 1 and 12")

    @staticmethod
    def _validate_day_values(df, year_col, month_col, day_col):
        """Validate that day values are valid for given year and month."""
        def is_valid_day(row):
            year, month, day = row[year_col], row[month_col], row[day_col]
            _, last_day = calendar.monthrange(year, month)
            return 1 <= day <= last_day
            
        if not df.apply(is_valid_day, axis=1).all():
            raise ValueError("Invalid day values for given year and month")

    @staticmethod
    def _validate_quarter_values(df, quarter_col):
        """Validate that quarter values are between 1 and 4."""
        if (df[quarter_col] < 1).any() or (df[quarter_col] > 4).any():
            raise ValueError("Quarter values must be between 1 and 4")

    @staticmethod
    def _add_date_with_last_day(df, year_col, month_col):
        """Add date column with last day of each month and set as index."""
        df['date'] = df.apply(
            lambda row: pd.Timestamp(
                year=int(row[year_col]),
                month=int(row[month_col]),
                day=calendar.monthrange(int(row[year_col]), int(row[month_col]))[1]
            ),
            axis=1
        )
        return df.set_index('date')

    @staticmethod
    def _add_date_with_day(df, year_col, month_col, day_col):
        """Add date column with specified day and set as index."""
        df['date'] = pd.to_datetime(df[[year_col, month_col, day_col]])
        return df.set_index('date')

    @staticmethod
    def _add_date_for_quarter(df, year_col, quarter_col):
        """Add date column with last day of each quarter and set as index."""
        df['date'] = df.apply(
            lambda row: pd.Timestamp(
                year=int(row[year_col]),
                month=int(row[quarter_col] * 3),
                day=calendar.monthrange(int(row[year_col]), int(row[quarter_col] * 3))[1]
            ),
            axis=1
        )
        return df.set_index('date')