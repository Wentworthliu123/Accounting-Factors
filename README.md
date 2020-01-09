## Before everything, technical reference for github
https://scotch.io/tutorials/how-to-push-an-existing-project-to-github


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

## Notes from the paper "Machine Valuation"
- Intuitively, investment predicts returns because given expected cash ﬂows, high costs of capital imply low net present values of new capital and low investment, and low costs of capital imply high net present values of new capitalandhighinvestment.ROEpredictsreturnsbecausehighexpectedROE relative to low investment must imply high discount rates. The high discount rates are necessary to offset the high expected ROE to induce low net present values of new capital and low investment. If the discount rates were not high enough, ﬁrms would instead observe high net present values of new capital and invest more. Conversely, low expected ROE relative to high investment must imply low discount rates. If the discount rates were not low enough to counteractthelowexpectedROE,ﬁrmswouldinsteadobservelownetpresent values of new capital and invest less.
- 折现率 = 盈利率/投资成本。相对于盈利率，投资越多
的公司，折现率越低，股票预期收益率也越低。相当于投
资，盈利率越高的公司，折现率越高，股票预期收益率也
越高。

## Translating SAS into Python
https://github.com/asnr/sas-to-python






