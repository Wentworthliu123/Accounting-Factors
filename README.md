## Study Log
### WRDS file list for compustat
https://wrds-web.wharton.upenn.edu/wrds/tools/variable.cfm?library_id=129&file_id=65811
### Must find the direct path for the dataset you wish you use, meaning on the same page
https://wrds-web.wharton.upenn.edu/wrds/ds/compd/funda/index.cfm?navId=83

STILL MISSING
### cik,substr(sic,1,2) as sic2, sic,naics
Standard Industrial Classification Code
Variable Name = SIC

North American Industrial Classification System
Variable Name = NAICS

### Compustat: COMPANY Dataset Vs. NAMES Dataset
https://wrds-www.wharton.upenn.edu/pages/support/support-articles/compustat/north-america/compustat-company-dataset-vs-names-dataset/?_ga=2.149312118.1772785986.1574277426-1924758175.1568907974

### Direct access to SIC NAICS
#### current segments
libname comp '/wrds/comp/sasdata/d_na/segments_current';
#### something that record company wise info but using historical fomrat
libname comp '/wrds/comp/sasdata/d_na/company';
#### The most close one that I found so far (11.23)
libname comp '/wrds/comp/sasdata/d_na';
comp.company can be simply replaced by comp.names, and cusip in select sql needs to be specified as f.cusip

### adsprate Company S&P Credit Ratings gvkey datadate 

## Notes from the paper "Machine Valuation"
- Guetal.(2018) assess the performance of a wide range of machine learning approaches in the context of predicting future stock returns. They show that machine learning techniques, in particular gradient boosted trees and neural networks, substantially outperform classic predictive methods based on least squares regressions.




