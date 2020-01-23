#coding=utf-8
#the first line is necessary to run this code on server

##########################################
# Accounting Factors translated from SAS 
# December 03 2019
# Created by Xinyu LIU
##########################################

import pandas as pd
import numpy as np
import datetime as dt
import wrds
import psycopg2 
import matplotlib.pyplot as plt
from dateutil.relativedelta import *
from pandas.tseries.offsets import *
from pandas.core.frame import DataFrame
from scipy import stats
import datetime
from matplotlib.backends.backend_pdf import PdfPages

###################
# Connect to WRDS #
###################
conn = wrds.Connection(wrds_username='dachxiu')
#make it a constant portal by creating ppass

#%%
"""
proc sql;
	create table data_0 as select 		/*header info*/
	substr(compress(f.cusip), 1, 6) as cnum, c.gvkey, datadate, datadate as 
		datadate_a, fyear, c.cik, cat(sic, 1, 2) as sic2, sic, naics, 
		/*firm variables*/
		/*income statement*/
		sale, revt, cogs, xsga, dp, xrd, xad, ib, ebitda, ebit, nopi, spi, pi, txp, 
		ni, txfed, txfo, txt, xint, /*CF statement and others*/
		capx, oancf, dvt, ob, gdwlia, gdwlip, gwo, /*assets*/
		rect, act, che, ppegt, invt, at, aco, intan, ao, ppent, gdwl, fatb, fatl, 
		/*liabilities*/
		lct, dlc, dltt, lt, dm, dcvt, cshrc, dcpstk, pstk, ap, lco, lo, drc, drlt, 
		txdi, /*equity and other*/
		ceq, scstkc, emp, csho, /*addition*/
		pstkrv, pstkl, txditc, year(datadate) as year, /*market*/
		abs(prcc_f) as prcc_f, csho*calculated prcc_f as mve_f, /*HXZ*/
		am, ajex, txdb, seq, dvc, dvp, dp, dvpsx_f, mib, ivao, ivst, sstk, prstkc, 
		dv, dltis, dltr, dlcch, oibdp, dvpa, tstkp, oiadp, xpp, xacc, re, ppenb, 
		ppenls, capxv, fopt, wcap from comp.names as c, comp.funda as f where 
		f.gvkey=c.gvkey
		/*get consolidated, standardized, industrial format statements*/
		and not missing(at) and not missing(prcc_f) and not missing(ni) and 
		datadate>='01Jan2015'd and f.indfmt='INDL' and f.datafmt='STD' and 
		f.popsrc='D' and f.consol='C';
quit;
"""

comp = conn.raw_sql("""
                    select 
                    f.cusip as cnum, c.gvkey, datadate, datadate as
                    datadate_a, fyear, c.cik, sic as sic2, sic, naics, 
                    sale, revt, cogs, xsga, xrd, xad, ib, ebitda, ebit, nopi, spi, pi, txp, 
                    ni, txfed, txfo, txt, xint, 
                    capx, oancf, dvt, ob, gdwlia, gdwlip, gwo, 
                    rect, act, che, ppegt, invt, at, aco, intan, ao, ppent, gdwl, fatb, fatl, 
                    lct, dlc, dltt, lt, dm, dcvt, cshrc, dcpstk, pstk, ap, lco, lo, drc, drlt, txdi,
                    ceq, scstkc, emp, csho, /*addition*/
                    pstkrv, pstkl, txditc, datadate as year, /*market*/
                    abs(prcc_f) as prcc_f, csho*prcc_f as mve_f, /*HXZ*/
                    am, ajex, txdb, seq, dvc, dvp, dp, dvpsx_f, mib, ivao, ivst, sstk, prstkc, 
                    dv, dltis, dltr, dlcch, oibdp, dvpa, tstkp, oiadp, xpp, xacc, re, ppenb, 
                    ppenls, capxv, fopt, wcap
                    from comp.names as c, comp.funda as f
                    where 
                    f.gvkey=c.gvkey
                    /*get consolidated, standardized, industrial format statements*/
                    and f.indfmt='INDL' 
                    and f.datafmt='STD'
                    and f.popsrc='D'
                    and f.consol='C'
                    and datadate >= '01/01/2015'
                    """)

comp.cnum=comp.cnum.replace(' ','').str.slice(0, 6)
comp.sic2=comp.sic2+'12'
comp.datadate=pd.to_datetime(comp.datadate)
comp.year = comp.datadate.dt.year
comp=comp.dropna(subset=['at','prcc_f','ni'])
#About same as SAS (line 75) (31304 rows in Python, 31318 rows in SAS)

#%%
"""
data data_1;
	set data_0;
	retain BE;
	ps_beme=coalesce(pstkrv, pstkl, pstk, 0);

	if missing(txditc) then
		txditc=0;
	BE=ceq + txditc - ps_beme;

	if BE<0 then
		BE=.;
run;

proc sort data=data_1 nodupkey;
	by gvkey datadate;
run;
"""
comp['ps_beme']=np.where(comp['pstkrv'].isnull(), comp['pstkl'], comp['pstkrv'])
comp['ps_beme']=np.where(comp['ps_beme'].isnull(),comp['pstk'], comp['ps_beme'])
comp['ps_beme']=np.where(comp['ps_beme'].isnull(),0,comp['ps_beme'])
comp['txditc']=comp['txditc'].fillna(0)
comp['be']=comp['ceq']+comp['txditc']-comp['ps_beme']
comp['be']=np.where(comp['be']>0,comp['be'],None)

comp=comp.sort_values(by=['gvkey','datadate']).drop_duplicates(['gvkey', 'datadate'])
#Sams as SAS line 93

#%%

# number of years in Compustat
comp['count']=comp.groupby(['gvkey']).cumcount()
#Sort DataFrame by column gvkey and datadate
#Mark cumulative number of each gvkey as of that row, starting from 0

"""
%crspmerge(s=m, start=01Jan2015, end=30jun2018, sfvars=&dsfvars, 
	sevars=&dsevars, filters=exchcd in (1, 2, 3) and shrcd in (10, 11));
"""

crsp_m = conn.raw_sql("""
                      select a.permno, a.permco, a.date, b.ticker, b.ncusip, b.shrcd, b.exchcd, b.siccd,
                      a.prc, a.ret, a.retx, a.shrout, a.vol
                      from crsp.msf as a
                      left join crsp.msenames as b
                      on a.permno=b.permno
                      and b.namedt<=a.date
                      and a.date<=b.nameendt
                      where a.date between '01/01/2015' and '06/30/2018'
                      and b.exchcd between 1 and 3
                      and b.shrcd between 10 and 11
                      """) 
#Note: b.dlprc does not exist
                      
crsp_m['me']=crsp_m['prc'].abs()*crsp_m['shrout']
crsp_m['prca']=crsp_m['prc'].abs()
crsp_m['lprc']=crsp_m.groupby(['permno','permco'])['prca'].shift(1)
crsp_m['lme']=crsp_m.groupby(['permno','permco'])['me'].shift(1)
#About as same as CRSP_M (Line 317) (154,792 rows in Python, 154,0112 in SAS)

#%%
"""
proc sql;
	create table crsp_m2 as select a.*, b.dlret, sum(1, ret)*sum(1, dlret)-1 as 
		retadj, abs(a.prc)*a.shrout as meq from CRSP_M as a left join 
		crsp.dsedelist (where=(missing(dlret)=0) ) as b on a.permno=b.permno and 
		intnx('day', a.date, 0, 'E')=intnx('day', b.DLSTDT, 0, 'E') order by permno, 
		date;
quit;
"""
# change variable format to int
crsp_m[['permco','permno','shrcd','exchcd']]=crsp_m[['permco','permno','shrcd','exchcd']].astype(int)

# Line up date to be end of month
crsp_m['date']=pd.to_datetime(crsp_m['date'])
crsp_m['jdate']=crsp_m['date']+MonthEnd(0)
#The 1 in MonthEnd just specifies to move one step forward to the next date that's a month end.

# add delisting return
dlret = conn.raw_sql("""
                     select permno, dlret, dlstdt 
                     from crsp.msedelist
                     """)
#MSEDELIST		CRSP Monthly Stock Event - Delisting
#DLRET 	Num	8	Delisting Return,DLRET is the return of the security after it is delisted. 
#It is calculated by comparing a value after delisting against the price on the security's last trading date. 
#The value after delisting can include a delisting price or the amount from a final distribution.
#DLSTDT 	Num	8	Delisting Date,DLSTDT contains the date (in YYMMDD format) of a security's last price on the current exchange.

#process dlret
dlret.permno=dlret.permno.astype(int)
dlret['dlstdt']=pd.to_datetime(dlret['dlstdt'])
dlret['jdate']=dlret['dlstdt']+MonthEnd(0)

#merge dlret and crsp_m
crsp = pd.merge(crsp_m, dlret, how='left',on=['permno','jdate'])
#crsp and dlret share the same column names: permno and jdate

#process crsp
crsp['dlret']=crsp['dlret'].fillna(0)
crsp['ret']=crsp['ret'].fillna(0)
crsp['retadj']=(1+crsp['ret'])*(1+crsp['dlret'])-1

crsp = crsp.sort_values(['permno', 'date'])

#process crsp
crsp=crsp.drop(['dlret','dlstdt'], axis=1)
crsp=crsp.sort_values(by=['jdate','permco','me']).drop_duplicates()

### Aggregate Market Cap ###
# sum of me across different permno belonging to same permco a given date
crsp_summe = crsp.groupby(['jdate','permco'])['me'].sum().reset_index()
# largest mktcap within a permco/date
crsp_maxme = crsp.groupby(['jdate','permco'])['me'].max().reset_index()
# join by jdate/maxme to find the permno
crsp1=pd.merge(crsp, crsp_maxme, how='inner', on=['jdate','permco','me'])
# drop me column and replace with the sum me
crsp1=crsp1.drop(['me'], axis=1)
# join with sum of me to get the correct market cap info
crsp2=pd.merge(crsp1, crsp_summe, how='inner', on=['jdate','permco'])
# sort by permno and date and also drop duplicates
crsp2=crsp2.sort_values(by=['permno','jdate']).drop_duplicates()
# important to have a duplicate check

#crsp same as CRSP_M (line 327)
#%%

# keep December market cap
crsp2['year']=crsp2['jdate'].dt.year
crsp2['month']=crsp2['jdate'].dt.month
decme=crsp2[crsp2['month']==12]
decme=decme[['permno','date','jdate','me','year']].rename(columns={'me':'dec_me'})

### July to June dates
crsp2['ffdate']=crsp2['jdate']+MonthEnd(-6)
crsp2['ffyear']=crsp2['ffdate'].dt.year
crsp2['ffmonth']=crsp2['ffdate'].dt.month
crsp2['1+retx']=1+crsp2['retx']
crsp2=crsp2.sort_values(by=['permno','date'])

# cumret by stock
crsp2['cumretx']=crsp2.groupby(['permno','ffyear'])['1+retx'].cumprod()
#cumprod returns the product of the year in this case, which is the cumulative return as time goes by

# lag cumret
crsp2['lcumretx']=crsp2.groupby(['permno'])['cumretx'].shift(1)

# lag market cap
crsp2['lme']=crsp2.groupby(['permno'])['me'].shift(1)

# if first permno then use me/(1+retx) to replace the missing value
crsp2['count']=crsp2.groupby(['permno']).cumcount()
crsp2['lme']=np.where(crsp2['count']==0, crsp2['me']/crsp2['1+retx'], crsp2['lme'])

# baseline me
mebase=crsp2[crsp2['ffmonth']==1][['permno','ffyear', 'lme']].rename(columns={'lme':'mebase'})

# merge result back together
crsp3=pd.merge(crsp2, mebase, how='left', on=['permno','ffyear'])
crsp3['wt']=np.where(crsp3['ffmonth']==1, crsp3['lme'], crsp3['mebase']*crsp3['lcumretx'])

decme['year']=decme['year']+1
decme=decme[['permno','year','dec_me']]

# Info as of June
crsp3_jun = crsp3[crsp3['month']==6]

crsp_jun = pd.merge(crsp3_jun, decme, how='inner', on=['permno','year'])

# Because I haven't reach the end of the code so will temporarily leave this subslicing open
# crsp_jun=crsp_jun[['permno','date', 'jdate', 'shrcd','exchcd','retadj','me','wt','cumretx','mebase','lme','dec_me']]
crsp_jun=crsp_jun.sort_values(by=['permno','jdate']).drop_duplicates()
crsp_jun=crsp_jun.drop(columns=['count'])

#decme about equal to DECME2 (10,861 in Python, 10,979 rows in SAS)
#crsp_jun about equal to CRSP_JUNE (10,496 rows in Python, 10,611 rows in SAS)
#crsp3 about equal to CRSP_M3 (152,687 rows in Python, 154,012 rows in SAS)
#%%

#######################
# CCM Block           #
#######################
ccm=conn.raw_sql("""
                  select gvkey, lpermno as permno, linktype, linkprim, 
                  linkdt, linkenddt
                  from crsp.ccmxpf_linktable
                  where substr(linktype,1,1)='L'
                  and linkprim in ('P', 'C')
                  """)

ccm['linkdt']=pd.to_datetime(ccm['linkdt'])
ccm['linkenddt']=pd.to_datetime(ccm['linkenddt'])
# if linkenddt is missing then set to next month's june (Changed)
ccm['linkenddt']=ccm['linkenddt'].fillna(pd.to_datetime('today')+YearEnd(0)+MonthEnd(6)) 
ccm['linkenddt']=ccm['linkenddt'].dt.date
ccm['linkenddt']=pd.to_datetime(ccm['linkenddt'])

ccm1=pd.merge(comp,ccm,how='left',on=['gvkey'])
ccm1['yearend']=ccm1['datadate']+YearEnd(0)
ccm1['jdate']=ccm1['yearend']+MonthEnd(6)

ccm2=ccm1[(ccm1['jdate']>=ccm1['linkdt'])&(ccm1['jdate']<=(ccm1['linkenddt']))]
ccm2=ccm2.drop(columns=['datadate_a','linktype','linkdt','linkenddt'])

ccm_data=pd.merge(ccm2,crsp_jun, how='left', on=['permno', 'jdate'])

"""
#Why did Xinyu do June? 

# # link comp and crsp
# Note: Different from SAS code, I left merge CCM2 to CRSP_JUN
# It could be exactly the same if using the form below
# ccm_jun = pd.merge( ccm2,crsp_jun, how='left', on=['permno', 'jdate'])
ccm_jun=pd.merge(crsp_jun, ccm2, how='inner', on=['permno', 'jdate'])

#filtering out prc==nan and dec_me==0
ccm_jun=ccm_jun[ccm_jun.dec_me!=0]
ccm_jun['beme']=ccm_jun['be']*1000/ccm_jun['dec_me']

# drop duplicates
ccm_jun=ccm_jun.sort_values(by=['permno','date']).drop_duplicates()
ccm_jun=ccm_jun.sort_values(by=['gvkey','date']).drop_duplicates()

# Note: Different from SAS, Python count start from zero, will see if I need to add 1 to better serve the need
ccm_jun['count']=ccm_jun.groupby(['gvkey']).cumcount()

# Parallel to the cleaning step for 'dr'
ccm_jun['dr']=np.where(ccm_jun.drc.notna() & ccm_jun.drlt.notna(),ccm_jun.drc+ccm_jun.drlt,None)
ccm_jun['dr']=np.where(ccm_jun.drc.notna() & ccm_jun.drlt.isna(),ccm_jun.drc,ccm_jun['dr'])
ccm_jun['dr']=np.where(ccm_jun.drc.isna() & ccm_jun.drlt.notna(),ccm_jun.drlt,ccm_jun['dr'])
# Parallel to the cleaning step for 'dc'
ccm_jun['dc']=np.where(ccm_jun.dcvt.isna() & ccm_jun.dcpstk.notna() & ccm_jun.pstk.notna() & (ccm_jun.dcpstk>ccm_jun.pstk),\
                       ccm_jun.dcpstk-ccm_jun.pstk,None)
ccm_jun['dc']=np.where(ccm_jun.dcvt.isna() & ccm_jun.dcpstk.notna() & ccm_jun.pstk.isna(),\
                       ccm_jun.dcpstk,ccm_jun['dr'])
ccm_jun['dc']=np.where(ccm_jun.dc.isna(), ccm_jun.dcvt, ccm_jun['dr'])
ccm_jun['xint']=ccm_jun['xint'].fillna(0)
ccm_jun['xsga']=ccm_jun['xsga'].fillna(0)

ccm_jun=ccm_jun.sort_values(by=['permno','date']).drop_duplicates()
"""

ccm_data=ccm_data[ccm_data.dec_me!=0]
ccm_data['beme']=ccm_data['be']*1000/ccm_data['dec_me']

# drop duplicates
ccm_data=ccm_data.sort_values(by=['permno','date']).drop_duplicates()
ccm_data=ccm_data.sort_values(by=['gvkey','date']).drop_duplicates()

# Note: Different from SAS, Python count start from zero, will see if I need to add 1 to better serve the need
ccm_data['count']=ccm_data.groupby(['gvkey']).cumcount()

# Parallel to the cleaning step for 'dr'
ccm_data['dr']=np.where(ccm_data.drc.notna() & ccm_data.drlt.notna(),ccm_data.drc+ccm_data.drlt,None)
ccm_data['dr']=np.where(ccm_data.drc.notna() & ccm_data.drlt.isna(),ccm_data.drc,ccm_data['dr'])
ccm_data['dr']=np.where(ccm_data.drc.isna() & ccm_data.drlt.notna(),ccm_data.drlt,ccm_data['dr'])
# Parallel to the cleaning step for 'dc'
ccm_data['dc']=np.where(ccm_data.dcvt.isna() & ccm_data.dcpstk.notna() & ccm_data.pstk.notna() & (ccm_data.dcpstk>ccm_data.pstk),\
                       ccm_data.dcpstk-ccm_data.pstk,None)
ccm_data['dc']=np.where(ccm_data['dcvt'].isna() & ccm_data['dcpstk'].notna() & ccm_data['pstk'].isna(),\
                       ccm_data['dcpstk'],ccm_data['dr'])
ccm_data['dc']=np.where(ccm_data.dc.isna(), ccm_data.dcvt, ccm_data['dr'])
ccm_data['xint']=ccm_data['xint'].fillna(0)
ccm_data['xsga']=ccm_data['xsga'].fillna(0)

ccm_data=ccm_data.sort_values(by=['permno','date']).drop_duplicates()
ccm_jun = ccm_data.copy()
#ccm_data about equal to DATA (20,139 rows Python vs 20,130 SAS)
#%%

#######################
# more clean-up and create first pass of variables           #
#######################
#create simple-just annual Compustat variables
#Line 532

ccm_jun['ep']=ccm_jun.ib/ccm_jun.mve_f
ccm_jun['cashpr']=(ccm_jun.mve_f+ccm_jun.dltt-ccm_jun['at'])/ccm_jun.che
ccm_jun['dy']=ccm_jun.dvt/ccm_jun.mve_f
ccm_jun['lev']=ccm_jun['lt']/ccm_jun.mve_f
ccm_jun['sp']=ccm_jun.sale/ccm_jun.mve_f
ccm_jun['roic']=(ccm_jun.ebit-ccm_jun.nopi)/(ccm_jun.ceq+ccm_jun['lt']-ccm_jun.che)
ccm_jun['rd_sale']=ccm_jun.xrd/ccm_jun.sale
ccm_jun['sp']=ccm_jun.sale/ccm_jun.mve_f #duplicate?

#Deleting duplicated columns
#ccm_jun = ccm_jun.loc[:,~ccm_jun.columns.duplicated()]

# treatment for lagged terms
ccm_jun['lagat']=ccm_jun.groupby(['permno'])['at'].shift(1)
ccm_jun['lagcsho']=ccm_jun.groupby(['permno'])['csho'].shift(1)
ccm_jun['laglt']=ccm_jun.groupby(['permno'])['lt'].shift(1)
ccm_jun['lagact']=ccm_jun.groupby(['permno'])['act'].shift(1)
ccm_jun['lagche']=ccm_jun.groupby(['permno'])['che'].shift(1)
ccm_jun['lagdlc']=ccm_jun.groupby(['permno'])['dlc'].shift(1)
ccm_jun['lagtxp']=ccm_jun.groupby(['permno'])['txp'].shift(1)
ccm_jun['laglct']=ccm_jun.groupby(['permno'])['lct'].shift(1)
ccm_jun['laginvt']=ccm_jun.groupby(['permno'])['invt'].shift(1)
ccm_jun['lagemp']=ccm_jun.groupby(['permno'])['emp'].shift(1)
ccm_jun['lagsale']=ccm_jun.groupby(['permno'])['sale'].shift(1)
ccm_jun['lagib']=ccm_jun.groupby(['permno'])['ib'].shift(1)
ccm_jun['lag2at']=ccm_jun.groupby(['permno'])['at'].shift(2)
ccm_jun['lagrect']=ccm_jun.groupby(['permno'])['rect'].shift(1)
ccm_jun['lagcogs']=ccm_jun.groupby(['permno'])['cogs'].shift(1)
ccm_jun['lagxsga']=ccm_jun.groupby(['permno'])['xsga'].shift(1)
ccm_jun['lagppent']=ccm_jun.groupby(['permno'])['ppent'].shift(1)
ccm_jun['lagdp']=ccm_jun.groupby(['permno'])['dp'].shift(1)
ccm_jun['lagxad']=ccm_jun.groupby(['permno'])['xad'].shift(1)
ccm_jun['lagppegt']=ccm_jun.groupby(['permno'])['ppegt'].shift(1)
ccm_jun['lagceq']=ccm_jun.groupby(['permno'])['ceq'].shift(1)
ccm_jun['lagcapx']=ccm_jun.groupby(['permno'])['capx'].shift(1)
ccm_jun['lag2capx']=ccm_jun.groupby(['permno'])['capx'].shift(2)
ccm_jun['laggdwl']=ccm_jun.groupby(['permno'])['gdwl'].shift(1)
ccm_jun['lagdvt']=ccm_jun.groupby(['permno'])['dvt'].shift(1)
ccm_jun['lagob']=ccm_jun.groupby(['permno'])['ob'].shift(1)
ccm_jun['lagaco']=ccm_jun.groupby(['permno'])['aco'].shift(1)
ccm_jun['lagintan']=ccm_jun.groupby(['permno'])['intan'].shift(1)
ccm_jun['lagao']=ccm_jun.groupby(['permno'])['ao'].shift(1)
ccm_jun['lagap']=ccm_jun.groupby(['permno'])['ap'].shift(1)
ccm_jun['laglco']=ccm_jun.groupby(['permno'])['lco'].shift(1)
ccm_jun['laglo']=ccm_jun.groupby(['permno'])['lo'].shift(1)
ccm_jun['lagdr']=ccm_jun.groupby(['permno'])['dr'].shift(1)
ccm_jun['lagxrd']=ccm_jun.groupby(['permno'])['xrd'].shift(1)
ccm_jun['lagni']=ccm_jun.groupby(['permno'])['ni'].shift(1)
ccm_jun['lagdltt']=ccm_jun.groupby(['permno'])['dltt'].shift(1)



ccm_jun['agr']=np.where(ccm_jun['at'].isna() | ccm_jun.lagat.isna(), np.NaN, (ccm_jun.lagat-ccm_jun['at'])/ccm_jun.lagat)
ccm_jun['gma']=(ccm_jun['revt']-ccm_jun['cogs'])/ccm_jun['lagat']
ccm_jun['chcsho']=ccm_jun.csho/ccm_jun.lagcsho -1
ccm_jun['lgr']=ccm_jun['lt']/ccm_jun.laglt -1
ccm_jun['acc']=(ccm_jun.ib-ccm_jun.oancf)/(ccm_jun['at']+ccm_jun.lagat) * 2

ccm_jun['pctacc']=np.where(ccm_jun['ib']==0,(ccm_jun['ib']-ccm_jun['oancf'])/0.01, np.NaN)
ccm_jun['pctacc']=np.where(ccm_jun['oancf'].isna(),(ccm_jun['act']-ccm_jun['lagact']-(ccm_jun['che']-ccm_jun['lagche']))\
                           -(ccm_jun['lct']-ccm_jun['laglct']-(ccm_jun['dlc']-ccm_jun['lagdlc'])\
                            -(ccm_jun['txp']-ccm_jun['lagtxp'])-ccm_jun['dp'])/ccm_jun['ib'].abs(), ccm_jun['pctacc'])
ccm_jun['pctacc']=np.where(ccm_jun['oancf'].isna() & ccm_jun['ib']==0, (ccm_jun['act']-ccm_jun['lagact']-(ccm_jun['che']-ccm_jun['lagche']))\
                           -(ccm_jun['lct']-ccm_jun['laglct']-(ccm_jun['dlc']-ccm_jun['lagdlc'])\
                            -(ccm_jun['txp']-ccm_jun['lagtxp'])-ccm_jun['dp'])/0.01, ccm_jun['pctacc'])

ccm_jun['cfp']= (ccm_jun['ib']-(ccm_jun['act']-ccm_jun['lagact']-(ccm_jun['che']-ccm_jun['lagche'])))\
                -(ccm_jun['lct']-ccm_jun['laglct']-(ccm_jun['dlc']-ccm_jun['lagdlc'])\
                  -(ccm_jun['txp']-ccm_jun['lagtxp'])-ccm_jun['dp'])/ccm_jun['mve_f']
ccm_jun['cfp']=np.where(ccm_jun['oancf'].notna(),ccm_jun['oancf']/ccm_jun['mve_f'], ccm_jun['cfp'])
ccm_jun['absacc']=ccm_jun['acc'].abs()
ccm_jun['chinv']=2*(ccm_jun['invt']-ccm_jun['laginvt'])/(ccm_jun['at']+ccm_jun['lagat'])
ccm_jun['spii']=np.where((ccm_jun['spi']!=0)&ccm_jun['spi'].notna(), 1, 0)

ccm_jun['spi']=2*ccm_jun['spi']/(ccm_jun['at']+ccm_jun['lagat'])
ccm_jun['cf']=2*ccm_jun['oancf']/(ccm_jun['at']+ccm_jun['lagat'])

ccm_jun['cf']=np.where(ccm_jun['oancf'].isna(), (ccm_jun['ib']-(ccm_jun['act']-ccm_jun['lagact']-(ccm_jun['che']-ccm_jun['lagche'])))\
                -(ccm_jun['lct']-ccm_jun['laglct']-(ccm_jun['dlc']-ccm_jun['lagdlc'])\
                  -(ccm_jun['txp']-ccm_jun['lagtxp'])-ccm_jun['dp'])/((ccm_jun['at']+ccm_jun['lagat'])/2),ccm_jun['cf'])  
ccm_jun['hire']=(ccm_jun['emp']-ccm_jun['lagemp'])/ccm_jun['lagemp']
ccm_jun['hire']=np.where(ccm_jun['emp'].isna() | ccm_jun['lagemp'].isna(), 0, ccm_jun['hire'])

ccm_jun['sgr']=ccm_jun['sale']/ccm_jun['lagsale'] -1
ccm_jun['chpm']=ccm_jun['ib']/ccm_jun['sale']-ccm_jun['lagib']/ccm_jun['lagsale']
ccm_jun['chato']=(ccm_jun['sale']/((ccm_jun['at']+ccm_jun['lagat'])/2)) - (ccm_jun['lagsale']/((ccm_jun['lagat']+ccm_jun['lag2at'])/2))
ccm_jun['pchsale_pchinvt']=((ccm_jun['sale']-(ccm_jun['lagsale']))/(ccm_jun['lagsale']))-((ccm_jun['invt']-(ccm_jun['laginvt']))/(ccm_jun['laginvt']))
ccm_jun['pchsale_pchrect']=((ccm_jun['sale']-(ccm_jun['lagsale']))/(ccm_jun['lagsale']))-((ccm_jun['rect']-(ccm_jun['lagrect']))/(ccm_jun['lagrect']))
ccm_jun['pchgm_pchsale']=(((ccm_jun['sale']-ccm_jun['cogs'])-((ccm_jun['lagsale'])-(ccm_jun['lagcogs'])))/((ccm_jun['lagsale'])-(ccm_jun['lagcogs'])))-((ccm_jun['sale']-(ccm_jun['lagsale']))/(ccm_jun['lagsale']))
ccm_jun['pchsale_pchxsga']=((ccm_jun['sale']-(ccm_jun['lagsale']))/(ccm_jun['lagsale']) )-((ccm_jun['xsga']\
        -(ccm_jun['lagxsga'])) /(ccm_jun['lagxsga']) )
ccm_jun['depr']=ccm_jun['dp']/ccm_jun['ppent']
ccm_jun['pchdepr']=((ccm_jun['dp']/ccm_jun['ppent'])-((ccm_jun['lagdp'])/(ccm_jun['lagppent'])))/((ccm_jun['lagdp'])/(ccm_jun['lagppent']))
ccm_jun['chadv']=np.log(1+ccm_jun['xad'])-np.log((1+(ccm_jun['lagxad'])))
ccm_jun['invest']=((ccm_jun['ppegt']-(ccm_jun['lagppegt'])) +  (ccm_jun['invt']-(ccm_jun['laginvt'])) ) / (ccm_jun['lagat'])
ccm_jun['invest']=np.where(ccm_jun['ppegt'].isna(), ((ccm_jun['ppent']-(ccm_jun['lagppent'])) +  (ccm_jun['invt']-(ccm_jun['laginvt'])) ) / (ccm_jun['lagat']), ccm_jun['invest'])
ccm_jun['egr']=((ccm_jun['ceq']-(ccm_jun['lagceq']))/(ccm_jun['lagceq']) )

# Note here instead of using count>=2, I use 1 to stay in line iwth python cumcount
# Also starting from here I'll keep the SAS in the notes for comparison and debug
    # 	if missing(capx) and count>=2 then
    # 		capx=ppent-lag(ppent);
    # 	pchcapx=(capx-lag(capx))/lag(capx);
    # 	grcapx=(capx-lag2(capx))/lag2(capx);
    # 	grGW=(gdwl-lag(gdwl))/lag(gdwl);
ccm_jun['capx']=np.where(ccm_jun['capx'].isna() & ccm_jun['count']>=1,ccm_jun['ppent']-(ccm_jun['lagppent']), ccm_jun['capx'])
ccm_jun['pchcapx']=(ccm_jun['capx']-ccm_jun['lagcapx'])/ccm_jun['lagcapx']
ccm_jun['grcapx']=(ccm_jun['capx']-ccm_jun['lag2capx'])/ccm_jun['lag2capx']
ccm_jun['grGW']=(ccm_jun['gdwl']-ccm_jun['laggdwl'])/ccm_jun['laggdwl']
    # 	if missing(gdwl) or gdwl=0 then
    # 		grGW=0;
    # 	if gdwl ne 0 and not missing(gdwl) and missing(grGW) then
    # 		grGW=1;
ccm_jun['grGW']=np.where(ccm_jun['gdwl'].isna() | ccm_jun['gdwl']==0, 0, ccm_jun['grGW'])  
ccm_jun['grGW']=np.where(ccm_jun['gdwl'].notna() & ccm_jun['gdwl']!=0 & ccm_jun['grGW'].isna(), 1, ccm_jun['grGW']) 
    # 	if (not missing(gdwlia) and gdwlia ne 0) or (not missing(gdwlip) and gdwlip ne 
    # 		0) or (not missing(gwo) and gwo ne 0) then
    # 			woGW=1;
ccm_jun['woGW']=np.where((ccm_jun['gdwlia'].notna()&ccm_jun['gdwlia']!=0)|(ccm_jun['gdwlip'].notna()&(ccm_jun['gdwlip']!=0))|\
                                                                          (ccm_jun['gwo'].notna()&ccm_jun['gwo']!=0) , 1, 0)
    #	tang=(che+rect*0.715+invt*0.547+ppent*0.535)/at;
ccm_jun['tang']=(ccm_jun['che']+ccm_jun['rect']*0.715+ccm_jun['invt']*0.547+ccm_jun['ppent']*0.535)/ccm_jun['at']

# 	if (2100<=sic<=2199) or (2080<=sic<=2085) or (naics in ('7132', '71312', 
# 		'713210', '71329', '713290', '72112', '721120')) then
# 			sin=1;
# 	else
# 		sin=0;

# 	if missing(act) then
# 		act=che+rect+invt;

# 	if missing(lct) then
# 		lct=ap;
ccm_jun['sic']=ccm_jun['sic'].astype(int)
ccm_jun['sin']=np.where(ccm_jun['sic'].between(2100,2199) | ccm_jun['sic'].between(2080,2085) | (ccm_jun['naics'].isin(['7132', '71312', \
                                                                             '713210', '71329', '713290', '72112', '721120'])), 1, 0)
ccm_jun['act']=np.where(ccm_jun['act'].isna(), ccm_jun['che']+ccm_jun['rect']+ccm_jun['invt'],ccm_jun['act'])
ccm_jun['lct']=np.where(ccm_jun['lct'].isna(), ccm_jun['ap'], ccm_jun['lct'])

# 	currat=act/lct;
# 	pchcurrat=((act/lct)-(lag(act)/lag(lct)))/(lag(act)/lag(lct));
# 	quick=(act-invt)/lct;
# 	pchquick=((act-invt)/lct - (lag(act)-lag(invt))/lag(lct) )/ 
# 		((lag(act)-lag(invt) )/lag(lct) );
# 	salecash=sale/che;
# 	salerec=sale/rect;
# 	saleinv=sale/invt;
# 	pchsaleinv=((sale/invt)-(lag(sale)/lag(invt)) ) / (lag(sale)/lag(invt));
# 	cashdebt=(ib+dp)/((lt+lag(lt))/2);
# 	realestate=(fatb+fatl)/ppegt;
    
ccm_jun['currat']=ccm_jun['act']/ccm_jun['lct']
ccm_jun['pchcurrat']=((ccm_jun['act']/ccm_jun['lct'])-((ccm_jun['lagact'])/(ccm_jun['laglct'])))/((ccm_jun['lagact'])/(ccm_jun['lct']))
ccm_jun['quick']=(ccm_jun['act']-ccm_jun['invt'])/ccm_jun['lct']
ccm_jun['pchquick']=((ccm_jun['act']-ccm_jun['invt'])/ccm_jun['lct'] - ((ccm_jun['lagact'])-(ccm_jun['laginvt']))/(ccm_jun['laglct']) )/ (((ccm_jun['lagact'])-(ccm_jun['laginvt']))/(ccm_jun['laglct']))
ccm_jun['salecash']=ccm_jun['sale']/ccm_jun['che']
ccm_jun['salerec']=ccm_jun['sale']/ccm_jun['rect']
ccm_jun['saleinv']=ccm_jun['sale']/ccm_jun['invt']
ccm_jun['pchsaleinv']=((ccm_jun['sale']/ccm_jun['invt'])-((ccm_jun['lagsale'])/(ccm_jun['laginvt'])) ) / ((ccm_jun['lagsale'])/(ccm_jun['laginvt']))
ccm_jun['cashdebt']=(ccm_jun['ib']+ccm_jun['dp'])/((ccm_jun['lt']+(ccm_jun['laglt']))/2)
ccm_jun['realestate']=(ccm_jun['fatb']+ccm_jun['fatl'])/ccm_jun['ppegt']

# 	if missing(ppegt) then
# 		realestate=(fatb+fatl)/ppent;
# 	if (not missing(dvt) and dvt>0) and (lag(dvt)=0 or missing(lag(dvt))) then
# 		divi=1;
# 	else
# 		divi=0;
# 	if (missing(dvt) or dvt=0) and (lag(dvt)>0 and not missing(lag(dvt))) then
# 		divo=1;
# 	else
# 		divo=0;
ccm_jun['realestate']=np.where(ccm_jun['ppegt'].isna(), (ccm_jun['fatb']+ccm_jun['fatl'])/ccm_jun['ppent'], ccm_jun['realestate'])
ccm_jun['divi']=np.where((ccm_jun['dvt'].notna() & ccm_jun['dvt']>0) & ((ccm_jun['lagdvt'])==0 | (ccm_jun['lagdvt'].isna())),1,0)
ccm_jun['divo']=np.where((ccm_jun['dvt'].isna() | ccm_jun['dvt']==0) & ((ccm_jun['lagdvt'])>0 & (ccm_jun['lagdvt'].notna())),1,0)

# 	obklg=ob/((at+lag(at))/2);
# 	chobklg=(ob-lag(ob))/((at+lag(at))/2);
ccm_jun['obklg']=ccm_jun['ob']/((ccm_jun['at']+(ccm_jun['lagat']))/2)
ccm_jun['chobklg']=(ccm_jun['ob']-(ccm_jun['lagob']))/((ccm_jun['at']+(ccm_jun['lagat']))/2)

# 	if not missing(dm) and dm ne 0 then
# 		securedind=1;
# 	else
# 		securedind=0;
# 	secured=dm/dltt;
# 	if not missing(dc) and dc ne 0 or (not missing(cshrc) and CSHRC ne 0) then
# 		convind=1;
# 	else
# 		convind=0;
# 	conv=dc/dltt;
ccm_jun['securedind']=np.where(ccm_jun['dm'].notna() &ccm_jun['dm']!=0, 1, 0)
ccm_jun['secured']=ccm_jun['dm']/ccm_jun['dltt']
ccm_jun['convind']=np.where((ccm_jun['dc'].notna() & ccm_jun['dc']!=0) | (ccm_jun['cshrc'].notna() & ccm_jun['cshrc']!=0) , 1, 0)
ccm_jun['dc']=ccm_jun['dc'].astype(float)
## There will be inf in the result
ccm_jun['conv']=ccm_jun['dc']/ccm_jun['dltt']

# 	grltnoa=((rect+invt+ppent+aco+intan+ao-ap-lco-lo)-(lag(rect)+lag(invt)+lag(ppent)+lag(aco)+lag(intan)+lag(ao)-lag(ap)-lag(lco)-lag(lo)) 
# 		-(rect-lag(rect)+invt-lag(invt)+aco-lag(aco)-(ap-lag(ap)+lco-lag(lco)) 
# 		-dp))/((at+lag(at))/2);
# 	chdrc=(dr-lag(dr))/((at+lag(at))/2);
ccm_jun['grltnoa']=((ccm_jun['rect']+ccm_jun['invt']+ccm_jun['ppent']+ccm_jun['aco']+ccm_jun['intan']+ccm_jun['ao']-ccm_jun['ap']-ccm_jun['lco']-ccm_jun['lo'])\
                    -((ccm_jun['lagrect'])+(ccm_jun['laginvt'])+(ccm_jun['lagppent'])+(ccm_jun['lagaco'])+(ccm_jun['lagintan'])+(ccm_jun['lagao'])-(ccm_jun['lagap'])\
                      -(ccm_jun['laglco'])-(ccm_jun['laglo'])) -(ccm_jun['rect']-(ccm_jun['lagrect'])+ccm_jun['invt']-(ccm_jun['laginvt'])+ccm_jun['aco']-(ccm_jun['lagaco'])\
                        -(ccm_jun['ap']-(ccm_jun['lagap'])+ccm_jun['lco']-(ccm_jun['laglco'])) -ccm_jun['dp']))/((ccm_jun['at']+(ccm_jun['lagat']))/2)
ccm_jun['chdrc']=(ccm_jun['dr']-(ccm_jun['lagdr']))/((ccm_jun['at']+(ccm_jun['lagat']))/2)

# 	if ((xrd/at)-(lag(xrd/lag(at))))/(lag(xrd/lag(at))) >.05 then
# 		rd=1;
# 	else
# 		rd=0;
# 	rdbias=(xrd/lag(xrd))-1-ib/lag(ceq);
# 	roe=ib/lag(ceq);
# additional process to carry on the next calculation
ccm_jun['xrd/lagat']=ccm_jun['xrd']/(ccm_jun['lagat'])
ccm_jun['lag(xrd/lagat)']=ccm_jun.groupby(['permno'])['xrd/lagat'].shift(1)
ccm_jun['rd']=np.where(((ccm_jun['xrd']/ccm_jun['at'])-ccm_jun['lag(xrd/lagat)'])/ccm_jun['lag(xrd/lagat)']>0.05, 1, 0)
ccm_jun['rdbias']=(ccm_jun['xrd']/(ccm_jun['lagxrd']))-1-ccm_jun['ib']/(ccm_jun['lagceq'])
ccm_jun['roe']=ccm_jun['ib']/(ccm_jun['lagceq'])

# 	ps_beme=coalesce(pstkrv, pstkl, pstk, 0);

# 	if missing(txditc) then
# 		txditc=0;
# 	BE=ceq + txditc - ps_beme;

# 	if BE<0 then
# 		BE=.;
# 	operprof=(revt-cogs-xsga-xint)/BE;

# 	if missing(revt) then
# 		operprof=.;

# 	if missing(cogs)=1 and missing(xsga)=1 and missing(xint)=1 then
# 		operprof=.;

# 	if missing(BE) then
# 		operprof=.;
ccm_jun['ps_beme']=np.where(ccm_jun['pstkrv'].isnull(), ccm_jun['pstkl'], ccm_jun['pstkrv'])
ccm_jun['ps_beme']=np.where(ccm_jun['ps_beme'].isnull(),ccm_jun['pstk'], ccm_jun['ps_beme'])
ccm_jun['ps_beme']=np.where(ccm_jun['ps_beme'].isnull(),0,ccm_jun['ps_beme'])
ccm_jun['txditc']=ccm_jun['txditc'].fillna(0)
ccm_jun['be']=ccm_jun['ceq']+ccm_jun['txditc']-ccm_jun['ps_beme']
ccm_jun['be']=np.where(ccm_jun['be']>0,ccm_jun['be'],np.NaN)
ccm_jun['operprof']=np.where(ccm_jun['be'].notna() & ccm_jun['revt'].notna() & (ccm_jun['cogs'].notna() | ccm_jun['xsga'].notna() | ccm_jun['xint'].notna()),\
                          (ccm_jun['revt']-ccm_jun['cogs']-ccm_jun['xsga']-ccm_jun['xint'])/ccm_jun['be'], np.nan)

#	ps=(ni>0)+(oancf>0)+(ni/at > lag(ni)/lag(at))+(oancf>ni)+(dltt/at < lag(dltt)/lag(at))+(act/lct > lag(act)/lag(lct)) 
#+((sale-cogs)/sale > (lag(sale)-lag(cogs))/lag(sale))+ (sale/at > lag(sale)/lag(at))+ (scstkc=0)
#!!! Not sure if the boolean adding is the same 
ccm_jun['ps']=(ccm_jun['ni']>0)+(ccm_jun['oancf']>0)+(ccm_jun['ni']/ccm_jun['at'] > (ccm_jun['lagni'])/(ccm_jun['lagat']))+(ccm_jun['oancf']>ccm_jun['ni'])+(ccm_jun['dltt']/ccm_jun['at'] < (ccm_jun['lagdltt'])/(ccm_jun['lagat']))\
                +(ccm_jun['act']/ccm_jun['lct'] > (ccm_jun['lagact'])/(ccm_jun['laglct'])) +((ccm_jun['sale']-ccm_jun['cogs'])/ccm_jun['sale'] > ((ccm_jun['lagsale'])-(ccm_jun['lagcogs']))/(ccm_jun['lagsale']))\
                + (ccm_jun['sale']/ccm_jun['at'] > (ccm_jun['lagsale'])/(ccm_jun['lagat']))+ (ccm_jun['scstkc']==0)

# 	if fyear<=1978 then
# 		tr=.48;

# 	if 1979<=fyear<=1986 then
# 		tr=.46;

# 	if fyear=1987 then
# 		tr=.4;

# 	if 1988<=fyear<=1992 then
# 		tr=.34;

# 	if 1993<=fyear then
# 		tr=.35;
# 	tb_1=((txfo+txfed)/tr)/ib;

# 	if missing(txfo) or missing(txfed) then
# 		tb_1=((txt-txdi)/tr)/ib;
# 	*they rank within industries;

def tr_fyear(row):
    if row['fyear']<=1978:
        value = 0.48
    elif row['fyear']<=1986:
        value = 0.46
    elif row['fyear']==1987:
        value = 0.4
    elif row['fyear']>=1988 and row['fyear']<=1992:
        value = 0.34
    elif row['fyear']>=1993:
        value = 0.35
    else:
        value=''
    return value
ccm_jun['tr']=ccm_jun.apply(tr_fyear, axis=1)
ccm_jun['tb_1']=((ccm_jun['txfo']+ccm_jun['txfed'])/ccm_jun['tr'])/ccm_jun['ib']
ccm_jun['tb_1']=np.where(ccm_jun['txfo'].isna() | ccm_jun['txfed'].isna(),((ccm_jun['txt']+ccm_jun['txdi'])/ccm_jun['tr'])/ccm_jun['ib'], ccm_jun['tb_1'])

# 	if (txfo+txfed>0 or txt>txdi) and ib<=0 then
# 		tb_1=1;
#!!! Caution that for condition, when using | and &, one must apply parenthesis 
ccm_jun['tb_1']=np.where(((ccm_jun['txfo']+ccm_jun['txfed'])>0 | (ccm_jun['txt']>ccm_jun['txdi'])) & (ccm_jun['ib']<=0), 1, ccm_jun['tb_1'])

# 	*variables that will be used in subsequent steps to get to final RPS;
# 	*--prep for for Mohanram (2005) score;
# 	roa=ni/((at+lag(at))/2);
# 	cfroa=oancf/((at+lag(at))/2);
ccm_jun['roa']=ccm_jun['ni']/((ccm_jun['at']+(ccm_jun['lagat']))/2)
ccm_jun['cfroa']=ccm_jun['oancf']/((ccm_jun['at']+(ccm_jun['lagat']))/2)

# 	if missing(oancf) then
# 		cfroa=(ib+dp)/((at+lag(at))/2);
# 	xrdint=xrd/((at+lag(at))/2);
# 	capxint=capx/((at+lag(at))/2);
# 	xadint=xad/((at+lag(at))/2);
ccm_jun['cfroa']=np.where(ccm_jun['oancf'].isna(),ccm_jun['ib']+ccm_jun['dp'] /((ccm_jun['at']+(ccm_jun['lagat']))/2), ccm_jun['cfroa'])
ccm_jun['xrdint']=ccm_jun['xrd']/((ccm_jun['at']+(ccm_jun['lagat']))/2)
ccm_jun['capxint']=ccm_jun['capx']/((ccm_jun['at']+(ccm_jun['lagat']))/2)
ccm_jun['xadint']=ccm_jun['xad']/((ccm_jun['at']+(ccm_jun['lagat']))/2)

# 	/*HXZ*/
# 	adm=xad/mve6b;
# 	gad=(xad-lag(xad))/lag(xad);
# 	rdm=xrd/mve6b;
# 	rds=xrd/sale;
# 	ol=(cogs+xsga)/at;
# 	rc_1=xrd+0.8*lag(xrd)+0.6*lag2(xrd)+0.4*lag3(xrd)+0.2*lag4(xrd);
    
#New lag terms for this section
ccm_jun['lag2xrd']=ccm_jun.groupby(['permno'])['lagxrd'].shift(1)
ccm_jun['lag3xrd']=ccm_jun.groupby(['permno'])['lag2xrd'].shift(1)
ccm_jun['lag4xrd']=ccm_jun.groupby(['permno'])['lag3xrd'].shift(1)

# Here I follow previous naming of mve6b as dec_me
ccm_jun['adm']=ccm_jun['xad']/ccm_jun['dec_me']  
ccm_jun['gad']=(ccm_jun['xad']-(ccm_jun['lagxad']))/(ccm_jun['lagxad'])
ccm_jun['rdm']=ccm_jun['xrd']/ccm_jun['dec_me']
ccm_jun['rds']=ccm_jun['xrd']/ccm_jun['sale']
ccm_jun['ol']=(ccm_jun['cogs']+ccm_jun['xsga'])/ccm_jun['at']
ccm_jun['rc_1']=ccm_jun['xrd']+0.8*(ccm_jun['lagxrd'])+0.6*(ccm_jun['lag2xrd'])+0.4*(ccm_jun['lag3xrd'])+0.2*(ccm_jun['lag4xrd'])

# 	cdd=dcvt/(dlc+dltt);
# 	roaq_a=ib/lag(at);
ccm_jun['cdd']=ccm_jun['dcvt']/(ccm_jun['dlc']+ccm_jun['dltt'])
ccm_jun['roaq_a']=ccm_jun['ib']/(ccm_jun['lagat'])
# roavol_1=std(roaq_a, lag(roaq_a), lag2(roaq_a), lag3(roaq_a), lag4(roaq_a)), lag5(roaq_a), lag6(roaq_a), lag7(roaq_a), lag8(roaq_a), lag9(roaq_a))
#!!! not pretty sure about this, need wider time span to test
#!!! I am using min_periods to temporarily bypass this issue (and it feels like the way SAS does this)
#!!! Pertaining the 10th lag term, whether or not using the same min_periods depends on the SAS array process condition,
#!!! The std function itself will not set nan if there is a missing value in the arguments, it will just ignore it
#!!! Now I take out "min_periods=1" from rolling() condition to match what SAS does, but you'd always like to remain cautious

# must reindex in order to set back the value
ccm_jun['roavol_1']=ccm_jun.groupby(['permno'])['roaq_a'].rolling(10).std().reset_index()['roaq_a']

# cs_1=(ib-(act-lag(act)-(lct-lag(lct))-(che-lag(che))+dlc-lag(dlc)))/lag(at)
ccm_jun['cs_1']=(ccm_jun['ib']-(ccm_jun['act']-(ccm_jun['lagact'])-(ccm_jun['lct']-(ccm_jun['laglct']))-(ccm_jun['che']-(ccm_jun['lagche']))+ccm_jun['dlc']-(ccm_jun['lagdlc'])))/(ccm_jun['lagat'])
# 	roavol_2=std(cs_1, lag(cs_1), lag2(cs_1), lag3(cs_1), lag4(cs_1), lag5(cs_1), 
# 		lag6(cs_1), lag7(cs_1), lag8(cs_1), lag9(cs_1));
# 	roavol_a=roavol_1/roavol_2;
ccm_jun['roavol_2']=ccm_jun.groupby(['permno'])['cs_1'].rolling(10).std().reset_index()['cs_1']
ccm_jun['roavol_a']=ccm_jun['roavol_1']/ccm_jun['roavol_2']

# 	if missing(gdwl) then
# 		gdwl=0;

# 	if missing(intan) then
# 		intan=0;
# 	ala=che+0.75*(act-che)-0.5*(at-act-gdwl-intan);
# 	alm=ala/(at+prcc_f*csho-ceq);
# 	ob_a=ob/(.5*at+.5*lag(at));
# 	x_3=capx/sale;
# 	cinvest_a=x_3/((lag(x_3)+lag2(x_3)+lag3(x_3))/3)-1;
# 	dpia=(ppegt-lag(ppegt)+invt-lag(invt))/lag(at);

# 	if missing(dlc) then
# 		dlc=0;

# 	if missing(dltt) then
# 		dltt=0;

# 	if missing(mib) then
# 		mib=0;

# 	if missing(pstk) then
# 		pstk=0;

# 	if missing(ceq) then
# 		ceq=0;

ccm_jun['gdwl']=ccm_jun['gdwl'].fillna(0)
ccm_jun['intan']=ccm_jun['intan'].fillna(0)
ccm_jun['ala']=ccm_jun['che']+0.75*(ccm_jun['act']-ccm_jun['che'])-0.5*(ccm_jun['at']-ccm_jun['act']-ccm_jun['gdwl']-ccm_jun['intan'])
ccm_jun['alm']=ccm_jun['ala']/(ccm_jun['at']+ccm_jun['prcc_f']*ccm_jun['csho']-ccm_jun['ceq'])
ccm_jun['ob_a']=ccm_jun['ob']/(0.5*ccm_jun['at']+0.5*(ccm_jun['lagat']))
ccm_jun['x_3']=ccm_jun['capx']/ccm_jun['sale']
#New lag terms for this section
ccm_jun['lagx_3']=ccm_jun.groupby(['permno'])['x_3'].shift(1)
ccm_jun['lag2x_3']=ccm_jun.groupby(['permno'])['lagx_3'].shift(1)
ccm_jun['lag3x_3']=ccm_jun.groupby(['permno'])['lag2x_3'].shift(1)

ccm_jun['cinvest_a']=ccm_jun['x_3']/(((ccm_jun['lagx_3'])+(ccm_jun['lag2x_3'])+(ccm_jun['lag3x_3']))/3)-1

ccm_jun['dlc']=ccm_jun['dlc'].fillna(0)
ccm_jun['dltt']=ccm_jun['dltt'].fillna(0)
ccm_jun['mib']=ccm_jun['mib'].fillna(0)
ccm_jun['pstk']=ccm_jun['pstk'].fillna(0)
ccm_jun['ceq']=ccm_jun['ceq'].fillna(0)

# 	noa=((at-che)-(at-dlc-dltt-mib-pstk-ceq))/lag(at);
# 	dnoa=noa-lag(noa);
# 	pchcapx3=capx/lag3(capx)-1;
# 	x_4=dlc+dltt;
# 	cdi=log(x_4/lag5(x_4));
# 	ivg=invt/lag(invt)-1;
# 	dcoa=(act-lag(act)-(che-lag(che)))/lag(at);
# 	dcol=(lct-lag(lct)-(dlc-lag(dlc)))/lag(at);
# 	dwc=(dcoa-dcol)/lag(at);
# 	dnca=(at-act-ivao-(lag(at)-lag(act)-lag(ivao)))/lag(at);
# 	dncl=(lt-lct-dltt-(lag(lt)-lag(lct)-lag(dltt)))/lag(at);
# 	dnco=(dnca-dncl)/lag(at);
# 	dfin=(ivst+ivao-dltt-dlc-pstk-(lag(ivst)+lag(ivao)-lag(dltt)-lag(dlc)-lag(pstk)))/lag(at);
# 	ta=(dwc+dnco+dfin)/lag(at);
# 	dsti=(ivst-lag(ivst))/lag(at);
# 	dfnl=(dltt+dlc+pstk-(lag(dltt)+lag(dlc)+lag(pstk)))/lag(at);
# 	egr_hxz=(ceq-lag(ceq))/lag(at);

ccm_jun['noa']=((ccm_jun['at']-ccm_jun['che'])-(ccm_jun['at']-ccm_jun['dlc']-ccm_jun['dltt']-ccm_jun['mib']-ccm_jun['pstk']-ccm_jun['ceq']))/(ccm_jun['lagat'])
ccm_jun['lagnoa']=ccm_jun.groupby(['permno'])['noa'].shift(1)
ccm_jun['dnoa']=ccm_jun['noa']-(ccm_jun['lagnoa'])
ccm_jun['lag3capx']=ccm_jun.groupby(['permno'])['capx'].shift(3)
ccm_jun['pchcapx3']=ccm_jun['capx']/(ccm_jun['lag3capx'])-1
ccm_jun['x_4']=ccm_jun['dlc']+ccm_jun['dltt']
ccm_jun['lag5x_4']=ccm_jun.groupby(['permno'])['capx'].shift(5)
ccm_jun['cdi']=np.log(ccm_jun['x_4']/(ccm_jun['lag5x_4']))

ccm_jun['ivg']=ccm_jun['invt']/(ccm_jun['laginvt'])-1
ccm_jun['dcoa']=(ccm_jun['act']-(ccm_jun['lagact'])-(ccm_jun['che']-(ccm_jun['lagche'])))/(ccm_jun['lagat'])
ccm_jun['dcol']=(ccm_jun['lct']-(ccm_jun['laglct'])-(ccm_jun['dlc']-(ccm_jun['lagdlc'])))/(ccm_jun['lagat'])
ccm_jun['dwc']=(ccm_jun['dcoa']-ccm_jun['dcol'])/(ccm_jun['lagat'])
ccm_jun['lagivao']=ccm_jun.groupby(['permno'])['ivao'].shift(1)
ccm_jun['dnca']=(ccm_jun['at']-ccm_jun['act']-ccm_jun['ivao']-((ccm_jun['lagat'])-(ccm_jun['lagact'])-(ccm_jun['lagivao'])))/(ccm_jun['lagat'])
ccm_jun['dncl']=(ccm_jun['lt']-ccm_jun['lct']-ccm_jun['dltt']-((ccm_jun['laglt'])-(ccm_jun['laglct'])-(ccm_jun['lagdltt'])))/(ccm_jun['lagat'])
ccm_jun['dnco']=(ccm_jun['dnca']-ccm_jun['dncl'])/(ccm_jun['lagat'])
ccm_jun['lagivst']=ccm_jun.groupby(['permno'])['ivst'].shift(1)
ccm_jun['lagpstk']=ccm_jun.groupby(['permno'])['pstk'].shift(1)
ccm_jun['dfin']=(ccm_jun['ivst']+ccm_jun['ivao']-ccm_jun['dltt']-ccm_jun['dlc']-ccm_jun['pstk']-((ccm_jun['lagivst'])+(ccm_jun['lagivao'])-(ccm_jun['lagdltt'])-(ccm_jun['lagdlc'])-(ccm_jun['lagpstk'])))/(ccm_jun['lagat'])
ccm_jun['ta']=(ccm_jun['dwc']+ccm_jun['dnco']+ccm_jun['dfin'])/(ccm_jun['lagat'])
ccm_jun['dsti']=(ccm_jun['ivst']-(ccm_jun['lagivst']))/(ccm_jun['lagat'])
ccm_jun['dfnl']=(ccm_jun['dltt']+ccm_jun['dlc']+ccm_jun['pstk']-((ccm_jun['lagdltt'])+(ccm_jun['lagdlc'])+(ccm_jun['lagpstk'])))/(ccm_jun['lagat'])
ccm_jun['egr_hxz']=(ccm_jun['ceq']-(ccm_jun['lagceq']))/(ccm_jun['lagat'])

# 	if missing(txp) then
# 		txp=0;
# 	poa=(act-lag(act)-(che-lag(che))-(lct-lag(lct)-(dlc-lag(dlc))-(txp-lag(txp)))-dp)/abs(ni);
# 	nef=(sstk-prstkc-dv)/((at+lag(at))/2);

# 	if missing(dlcch) then
# 		dlcch=0;
# 	ndf=(dltis-dltr+dlcch)/((at+lag(at))/2);
# 	nxf=nef+ndf;
# 	atm=at/mve6b;
# 	cp=(ib+dp)/mve6b;
# 	op=(dvc+prstkc-(pstkrv-lag(pstkrv)))/mve6b;
# 	nop=(dvc+prstkc-(pstkrv-lag(pstkrv))-sstk+pstkrv-lag(pstkrv))/mve6b;
# 	em=(mve6b+dlc+dltt+pstkrv-che)/oibdp;

ccm_jun['lagpstkrv']=ccm_jun.groupby(['permno'])['pstkrv'].shift(1)

ccm_jun['txp']=ccm_jun['txp'].fillna(0)
ccm_jun['poa']=(ccm_jun['act']-(ccm_jun['lagact'])-(ccm_jun['che']-(ccm_jun['lagche']))-(ccm_jun['lct']-(ccm_jun['laglct'])-(ccm_jun['dlc']-(ccm_jun['lagdlc']))-(ccm_jun['txp']-(ccm_jun['lagtxp'])))-ccm_jun['dp'])/(ccm_jun['ni'].abs())
ccm_jun['nef']=(ccm_jun['sstk']-ccm_jun['prstkc']-ccm_jun['dv'])/((ccm_jun['at']+(ccm_jun['lagat']))/2)
ccm_jun['dlcch']=ccm_jun['dlcch'].fillna(0)
ccm_jun['ndf']=(ccm_jun['dltis']-ccm_jun['dltr']+ccm_jun['dlcch'])/((ccm_jun['at']+(ccm_jun['lagat']))/2)
ccm_jun['atm']=ccm_jun['at']/ccm_jun['dec_me']
ccm_jun['cp']=(ccm_jun['ib']+ccm_jun['dp'])/ccm_jun['dec_me']
ccm_jun['op']=(ccm_jun['dvc']+ccm_jun['prstkc']-(ccm_jun['pstkrv']-(ccm_jun['lagpstkrv'])))/ccm_jun['dec_me']
ccm_jun['nop']=(ccm_jun['dvc']+ccm_jun['prstkc']-(ccm_jun['pstkrv']-(ccm_jun['lagpstkrv']))-ccm_jun['sstk']+ccm_jun['pstkrv']-(ccm_jun['lagpstkrv']))/ccm_jun['dec_me']

# 	if missing(dvpa) then
# 		dvpa=0;

# 	if missing(tstkp) then
# 		tstkp=0;
# 	ndp=dltt+dlc+pstk+dvpa-tstkp-che;
# 	ebp=(ndp+ceq+tstkp-dvpa)/(ndp+mve6b);
# 	rna=oiadp/lag(noa);
# 	pm=rna/sale;
# 	ato=sale/lag(noa);
# 	cto=sale/lag(at);
# 	gpa=(revt-cogs)/at;
    
# 	rmw=(revt-cogs-xsga-xint)/(ceq+pstk);
# 	ole=(revt-cogs-xsga-xint)/(lag(ceq)+lag(pstk));
# 	opa=(revt-cogs-xsga+xrd)/at;
# 	ola=(revt-cogs-xsga+xrd)/lag(at);
# 	cop=(revt-cogs-xsga+xrd-(rect-lag(rect))-(invt-lag(invt))-(xpp-lag(xpp))+drc-lag(drc)+drlt-lag(drlt)+ap-lag(ap)+xacc-lag(xacc))/at;
# 	cla=(revt-cogs-xsga+xrd-(rect-lag(rect))-(invt-lag(invt))-(xpp-lag(xpp))+drc-lag(drc)+drlt-lag(drlt)+ap-lag(ap)+xacc-lag(xacc))/lag(at);

ccm_jun['dvpa']=ccm_jun['dvpa'].fillna(0)
ccm_jun['tstkp']=ccm_jun['tstkp'].fillna(0)
ccm_jun['ndp']=ccm_jun['dltt']+ccm_jun['dlc']+ccm_jun['pstk']+ccm_jun['dvpa']-ccm_jun['tstkp']-ccm_jun['che']
ccm_jun['ebp']=(ccm_jun['ndp']+ccm_jun['ceq']+ccm_jun['tstkp']-ccm_jun['dvpa'])/(ccm_jun['ndp']+ccm_jun['dec_me'])
ccm_jun['rna']=ccm_jun['oiadp']/(ccm_jun['lagnoa'])
ccm_jun['pm']=ccm_jun['rna']/ccm_jun['sale']
ccm_jun['ato']=ccm_jun['sale']/(ccm_jun['lagnoa'])
ccm_jun['cto']=ccm_jun['sale']/(ccm_jun['lagat'])
ccm_jun['gpa']=(ccm_jun['revt']-ccm_jun['cogs'])/ccm_jun['at']
ccm_jun['rmw']=(ccm_jun['revt']-ccm_jun['cogs']-ccm_jun['xsga']-ccm_jun['xint'])/(ccm_jun['ceq']+ccm_jun['pstk'])
ccm_jun['ole']=(ccm_jun['revt']-ccm_jun['cogs']-ccm_jun['xsga']-ccm_jun['xint'])/((ccm_jun['lagceq'])+(ccm_jun['lagpstk']))
ccm_jun['opa']=(ccm_jun['revt']-ccm_jun['cogs']-ccm_jun['xsga']+ccm_jun['xrd'])/ccm_jun['at']
ccm_jun['ola']=(ccm_jun['revt']-ccm_jun['cogs']-ccm_jun['xsga']+ccm_jun['xrd'])/(ccm_jun['lagat'])

ccm_jun['lagxpp']=ccm_jun.groupby(['permno'])['xpp'].shift(1)
ccm_jun['lagdrc']=ccm_jun.groupby(['permno'])['drc'].shift(1)
ccm_jun['lagdrlt']=ccm_jun.groupby(['permno'])['drlt'].shift(1)
ccm_jun['lagxacc']=ccm_jun.groupby(['permno'])['xacc'].shift(1)

ccm_jun['cop']=(ccm_jun['revt']-ccm_jun['cogs']-ccm_jun['xsga']+ccm_jun['xrd']-(ccm_jun['rect']-(ccm_jun['lagrect']))-(ccm_jun['invt']-(ccm_jun['laginvt']))-\
         (ccm_jun['xpp']-(ccm_jun['lagxpp']))+ccm_jun['drc']-(ccm_jun['lagdrc'])+ccm_jun['drlt']-(ccm_jun['lagdrlt'])+ccm_jun['ap']-(ccm_jun['lagap'])+ccm_jun['xacc']-(ccm_jun['lagxacc']))/ccm_jun['at']
ccm_jun['cla']=(ccm_jun['revt']-ccm_jun['cogs']-ccm_jun['xsga']+ccm_jun['xrd']-(ccm_jun['rect']-(ccm_jun['lagrect']))-(ccm_jun['invt']-(ccm_jun['laginvt']))-\
         (ccm_jun['xpp']-(ccm_jun['lagxpp']))+ccm_jun['drc']-(ccm_jun['lagdrc'])+ccm_jun['drlt']-(ccm_jun['lagdrlt'])+ccm_jun['ap']-(ccm_jun['lagap'])+ccm_jun['xacc']-(ccm_jun['lagxacc']))/ccm_jun['lagat']
    
#%%
    
    # 	if lt>at then
# 		i_1=1;
# 	else
# 		i_1=0;

# 	if ni<0 and lag(ni)<0 then
# 		i_2=1;
# 	else
# 		i_2=0;
# 	os=-1.32-0.407*log(at)+6.03*(dlc+dltt)/at-1.43*(act-lct)/at+0.076*(lct/act)-1.72*i_1-2.37*ni/at-1.83*(pi+dp)/lt+0.285*i_2-0.521*(ni+lag(ni))/(abs(ni)+abs(lag(ni)));
# 	zs=1.2*(act-lct)/at+1.4*re/at+3.3*oiadp/at+0.6*mve6b/lt+sale/at;
ccm_jun['i_1']=np.where(ccm_jun['lt']>ccm_jun['at'],1,0)
ccm_jun['i_2']=np.where((ccm_jun['ni']<0) & (ccm_jun['lagni']<0),1,0)
ccm_jun['os']=-1.32-0.407*np.log(ccm_jun['at'])+6.03*(ccm_jun['dlc']+ccm_jun['dltt'])/ccm_jun['at']-1.43*(ccm_jun['act']-ccm_jun['lct'])/ccm_jun['at']+0.076*(ccm_jun['lct']/ccm_jun['act'])-1.72*ccm_jun['i_1']-2.37*ccm_jun['ni']/ccm_jun['at']-1.83*(ccm_jun['pi']+ccm_jun['dp'])/ccm_jun['lt']+0.285*ccm_jun['i_2']-0.521*(ccm_jun['ni']+(ccm_jun['lagni']))/((ccm_jun['ni'].abs())+ccm_jun['lagni'].abs())
ccm_jun['zs']=1.2*(ccm_jun['act']-ccm_jun['lct'])/ccm_jun['at']+1.4*ccm_jun['re']/ccm_jun['at']+3.3*ccm_jun['oiadp']/ccm_jun['at']+0.6*ccm_jun['dec_me']/ccm_jun['lt']+ccm_jun['sale']/ccm_jun['at']

# 	if BE>0 then
# 		bi=at/BE;

# 	if missing(bi) then
# 		bi=.;

ccm_jun['bi']=np.where(ccm_jun['be']>0,ccm_jun['at']/ccm_jun['be'],np.NaN)
# 	pchsale_pchinvt_hxz=(sale-lag(sale))/(0.5*sale+0.5*lag(sale))-(invt-lag(invt))/(0.5*invt+0.5*lag(invt));
# 	pchsale_pchrect_hxz=(sale-lag(sale))/(0.5*sale+0.5*lag(sale))-(rect-lag(rect))/(0.5*rect+0.5*lag(rect));
# 	gm_1=sale-cogs;
# 	pchgm_pchsale_hxz=(gm_1-lag(gm_1))/(0.5*(gm_1+lag(gm_1)))-(sale-lag(sale))/(0.5*sale+0.5*lag(sale));
# 	pchsale_pchxsga_hxz=(sale-lag(sale))/(0.5*sale+0.5*lag(sale))-(xsga-lag(xsga))/(0.5*xsga+0.5*lag(xsga));
# 	realestate_hxz=(fatb+fatl)/ppegt;
# 	if missing(fatb) then
# 		realestate_hxz=(ppenb+ppenls)/ppent;

ccm_jun['pchsale_pchinvt_hxz']=(ccm_jun['sale']-(ccm_jun['lagsale']))/(0.5*ccm_jun['sale']+0.5*(ccm_jun['lagsale']))-(ccm_jun['invt']-(ccm_jun['laginvt']))/(0.5*ccm_jun['invt']+0.5*(ccm_jun['laginvt']))
ccm_jun['pchsale_pchrect_hxz']=(ccm_jun['sale']-(ccm_jun['lagsale']))/(0.5*ccm_jun['sale']+0.5*(ccm_jun['lagsale']))-(ccm_jun['rect']-(ccm_jun['lagrect']))/(0.5*ccm_jun['rect']+0.5*(ccm_jun['lagrect']))
ccm_jun['gm_1']=ccm_jun['sale']-ccm_jun['cogs']
ccm_jun['laggm_1']=ccm_jun.groupby(['permno'])['gm_1'].shift(1)
ccm_jun['pchgm_pchsale_hxz']=(ccm_jun['gm_1']-(ccm_jun['laggm_1']))/(0.5*(ccm_jun['gm_1']+(ccm_jun['laggm_1'])))-(ccm_jun['sale']-(ccm_jun['lagsale']))/(0.5*ccm_jun['sale']+0.5*(ccm_jun['lagsale']))
ccm_jun['pchsale_pchxsga_hxz']=(ccm_jun['sale']-(ccm_jun['lagsale']))/(0.5*ccm_jun['sale']+0.5*(ccm_jun['lagsale']))-(ccm_jun['xsga']-(ccm_jun['lagxsga']))/(0.5*ccm_jun['xsga']+0.5*(ccm_jun['lagxsga']))
ccm_jun['realestate_hxz']=(ccm_jun['fatb']+ccm_jun['fatl'])/ccm_jun['ppegt']
#ccm_jun['realestate_hxz']=np.where(ccm_jun['fatb'].isna(), (ccm_jun['ppenb']+ccm_jun['ppenls'])/ccm_jun['ppent'], ccm_jun['realestate_hxz'])    
#TODO: ppenb and ppenls are all nan's

# 	secured_hxz=dm/(dltt+dlc);
# 	agr_hxz=at/lag(at)-1;
# 	x_1=ppent+intan+ao-lo+dp;
# 	grltnoa_hxz=(x_1-lag(x_1))/((at+lag(at))/2);
# 	chcsho_hxz=log(csho*ajex)-log(lag(csho)*lag(ajex));
# 	pchcapx_hxz=(capxv-0.5*lag(capxv)-0.5*lag2(capxv))/(0.5*lag(capxv)+0.5*lag2(capxv));
ccm_jun['secured_hxz']=ccm_jun['dm']/(ccm_jun['dltt']+ccm_jun['dlc'])
ccm_jun['agr_hxz']=ccm_jun['at']/(ccm_jun['lagat'])-1
ccm_jun['x_1']=ccm_jun['ppent']+ccm_jun['intan']+ccm_jun['ao']-ccm_jun['lo']+ccm_jun['dp']
ccm_jun['lagx_1']=ccm_jun.groupby(['permno'])['x_1'].shift(1)
ccm_jun['grltnoa_hxz']=(ccm_jun['x_1']-(ccm_jun['lagx_1']))/((ccm_jun['at']+(ccm_jun['lagat']))/2)
ccm_jun['lagajex']=ccm_jun.groupby(['permno'])['ajex'].shift(1)
ccm_jun['chcsho_hxz']=np.log(ccm_jun['csho']*ccm_jun['ajex'])-np.log((ccm_jun['lagcsho'])*(ccm_jun['lagajex']))
ccm_jun['lagcapxv']=ccm_jun.groupby(['permno'])['capxv'].shift(1)
ccm_jun['lag2capxv']=ccm_jun.groupby(['permno'])['capxv'].shift(2)
ccm_jun['pchcapx_hxz']=(ccm_jun['capxv']-0.5*(ccm_jun['lagcapxv'])-0.5*(ccm_jun['lag2capxv']))/(0.5*(ccm_jun['lagcapxv'])+0.5*(ccm_jun['lag2capxv']))

# 	if missing(txp) then
# 		txp=0;
# 	acc_hxz=(act-lag(act)-(che-lag(che))-(lct-lag(lct)-(dlc-lag(dlc))-(txp-lag(txp)))-dp)/lag(at);
# 	pctacc_hxz=(dwc+dnco+dfin)/abs(ni);
# 	lev_hxz=(dlc+dltt)/mve6b;
# 	ep_hxz=ib/mve6b;
# 	cfp_hxz=(fopt-(wcap-lag(wcap)))/mve6b;
ccm_jun['txp']=ccm_jun['txp'].fillna(0)
ccm_jun['acc_hxz']=(ccm_jun['act']-(ccm_jun['lagact'])-(ccm_jun['che']-(ccm_jun['lagche']))-(ccm_jun['lct']-(ccm_jun['laglct'])-(ccm_jun['dlc']-(ccm_jun['lagdlc']))-(ccm_jun['txp']-(ccm_jun['lagtxp'])))-ccm_jun['dp'])/(ccm_jun['lagat'])
ccm_jun['pctacc_hxz']=(ccm_jun['dwc']+ccm_jun['dnco']+ccm_jun['dfin'])/(ccm_jun['ni'].abs())
ccm_jun['lev_hxz']=(ccm_jun['dlc']+ccm_jun['dltt'])/ccm_jun['dec_me']
ccm_jun['ep_hxz']=ccm_jun['ib']/ccm_jun['dec_me']
ccm_jun['lagwcap']=ccm_jun.groupby(['permno'])['wcap'].shift(1)
ccm_jun['cfp_hxz']=(ccm_jun['fopt']-(ccm_jun['wcap']-(ccm_jun['lagwcap'])))/ccm_jun['dec_me']

# 	if not missing(oancf) then
# 		cfp_hxz=(fopt-oancf)/mve6b;
# 	tb_hxz=pi/ni;
ccm_jun['cfp_hxz']=np.where(ccm_jun['oancf'].notna(),(ccm_jun['fopt']-ccm_jun['oancf'])/ccm_jun['dec_me'],ccm_jun['cfp_hxz'])
ccm_jun['tb_hxz']=ccm_jun['pi']/ccm_jun['ni']


#%%

# /*other preparation steps for annual variables: industry adjustments*/
# proc sql;
# 	create table data2 as select *, chpm-mean(chpm) as chpmia, chato-mean(chato) 
# 		as chatoia, sum(sale) as indsale, hire-mean(hire) as chempia, bm-mean(bm) as 
# 		bm_ia, pchcapx-mean(pchcapx) as pchcapx_ia, tb_1-mean(tb_1) as tb, 
# 		cfp-mean(cfp) as cfp_ia, mve_f-mean(mve_f) as mve_ia, /*HXZ*/
# 		sum(at) as indat, sum(BE) as indbe, pchcapx_hxz-mean(pchcapx_hxz) as 
# 		pchcapx_ia_hxz from data2 group by sic2, fyear;
# quit;
ccm_jun = ccm_jun.sort_values(['sic2', 'fyear'])
a=ccm_jun.groupby(['sic2','fyear'])['chpm'].mean()
a=a.rename('meanchpm')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])
ccm_jun['chpmia']=ccm_jun['chpm']-ccm_jun['meanchpm']

a=ccm_jun.groupby(['sic2','fyear'])['chato'].mean()
a=a.rename('meanchato')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])
ccm_jun['chatoia']=ccm_jun['chato']-ccm_jun['meanchato']

a=ccm_jun.groupby(['sic2','fyear'])['sale'].sum()
a=a.rename('indsale')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])

a=ccm_jun.groupby(['sic2','fyear'])['hire'].mean()
a=a.rename('meanhire')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])
ccm_jun['chempia']=ccm_jun['hire']-ccm_jun['meanhire']

ccm_jun['beme']=ccm_jun['beme'].astype(float)
a=ccm_jun.groupby(['sic2','fyear'])['beme'].mean()
a=a.rename('meanbeme')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])
ccm_jun['beme_ia']=ccm_jun['beme']-ccm_jun['meanbeme']

a=ccm_jun.groupby(['sic2','fyear'])['pchcapx'].mean()
a=a.rename('meanpchcapx')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])
ccm_jun['pchcapx_ia']=ccm_jun['pchcapx']-ccm_jun['meanpchcapx']

a=ccm_jun.groupby(['sic2','fyear'])['tb_1'].mean()
a=a.rename('meantb_1')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])
ccm_jun['tb']=ccm_jun['tb_1']-ccm_jun['meantb_1']

a=ccm_jun.groupby(['sic2','fyear'])['cfp'].mean()
a=a.rename('meancfp')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])
ccm_jun['cfp_ia']=ccm_jun['cfp']-ccm_jun['meancfp']

a=ccm_jun.groupby(['sic2','fyear'])['mve_f'].mean()
a=a.rename('meanmve_f')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])
ccm_jun['mve_ia']=ccm_jun['mve_f']-ccm_jun['meanmve_f']

a=ccm_jun.groupby(['sic2','fyear'])['at'].sum()
a=a.rename('indat')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])

a=ccm_jun.groupby(['sic2','fyear'])['be'].sum()
a=a.rename('indbe')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])

a=ccm_jun.groupby(['sic2','fyear'])['pchcapx_hxz'].mean()
a=a.rename('meanpchcapx_hxz')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])
ccm_jun['pchcapx_ia_hxz']=ccm_jun['pchcapx_hxz']-ccm_jun['meanpchcapx_hxz']

# proc sql;
# 	create table data2 as select *, sum((sale/indsale)*(sale/indsale) ) as herf, 
# 		/*HXZ*/
# 		sum((at/indat)*(at/indat)) as ha, sum((BE/indbe)*(BE/indbe)) as he from data2 
# 		group by sic2, fyear;
# quit;
ccm_jun['herfraw']=(ccm_jun['sale']/ccm_jun['indsale'])*(ccm_jun['sale']/ccm_jun['indsale'])
a=ccm_jun.groupby(['sic2','fyear'])['herfraw'].sum()
a=a.rename('herf')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])

ccm_jun['haraw']=(ccm_jun['at']/ccm_jun['indat'])*(ccm_jun['at']/ccm_jun['indat'])
a=ccm_jun.groupby(['sic2','fyear'])['haraw'].sum()
a=a.rename('ha')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])

ccm_jun['heraw']=(ccm_jun['be']/ccm_jun['indbe'])*(ccm_jun['be']/ccm_jun['indbe'])
a=ccm_jun.groupby(['sic2','fyear'])['heraw'].sum()
a=a.rename('he')
ccm_jun=pd.merge(ccm_jun, a, how='left', on=['sic2','fyear'])


# proc sort data=data2;
# 	by fyear sic2;
# run;
ccm_jun = ccm_jun.sort_values(['sic2', 'fyear'])

# proc univariate data=data2 noprint;
# 	by fyear sic2;
# 	var roa cfroa xrdint capxint xadint;
# 	output out=indmd median=md_roa md_cfroa md_xrdint md_capxint md_xadint;
# run;

# proc sql;
# 	create table data2 as select * from data2 a left join indmd b on 
# 		a.fyear=b.fyear and a.sic2=b.sic2;
# quit;

indmd=pd.DataFrame()

indmd = ccm_jun.groupby(['sic2', 'fyear'], as_index=False)['roa'].median()
indmd = indmd.rename(columns={'roa':'md_roa'})
indmd['md_cfroa'] = ccm_jun.groupby(['sic2', 'fyear'], as_index=False)['cfroa'].median()['cfroa']
indmd['md_xrdint'] = ccm_jun.groupby(['sic2', 'fyear'], as_index=False)['xrdint'].median()['xrdint']
indmd['md_capxint'] = ccm_jun.groupby(['sic2', 'fyear'], as_index=False)['capxint'].median()['capxint']
indmd['md_xadint'] = ccm_jun.groupby(['sic2', 'fyear'], as_index=False)['xadint'].median()['xadint']
ccm_jun=pd.merge(ccm_jun, indmd, how='left', on=['sic2','fyear'])

# proc sort data=data2 nodupkey;
# 	by gvkey datadate;
# run;

# data data2;
# 	set data2;
# 	*more for Mohanram score;

# 	if roa>md_roa then
# 		m1=1;
# 	else
# 		m1=0;

# 	if cfroa>md_cfroa then
# 		m2=1;
# 	else
# 		m2=0;

# 	if oancf>ni then
# 		m3=1;
# 	else
# 		m3=0;

# 	if xrdint>md_xrdint then
# 		m4=1;
# 	else
# 		m4=0;

# 	if capxint>md_capxint then
# 		m5=1;
# 	else
# 		m5=0;

# 	if xadint>md_xadint then
# 		m6=1;
# 	else
# 		m6=0;
# 	*still need to add another thing for Mohanram (2005) score;
# run;

ccm_jun = ccm_jun.sort_values(['gvkey', 'datadate'])
ccm_jun['m1']=np.where(ccm_jun['roa']>ccm_jun['md_roa'], 1, 0)
ccm_jun['m2']=np.where(ccm_jun['cfroa']>ccm_jun['md_cfroa'], 1, 0)
ccm_jun['m3']=np.where(ccm_jun['oancf']>ccm_jun['ni'], 1, 0)
ccm_jun['m4']=np.where(ccm_jun['xrdint']>ccm_jun['md_xrdint'], 1, 0)
ccm_jun['m5']=np.where(ccm_jun['capxint']>ccm_jun['md_capxint'], 1, 0)
ccm_jun['m6']=np.where(ccm_jun['xadint']>ccm_jun['md_xadint'], 1, 0)


#%%
# *----add credit rating--------;
# proc sql;
# 	create table data2 as select a.*, b.splticrm from data2 a left join 
# 		compr.adsprate b on a.gvkey=b.gvkey and year(a.datadate)=year(b.datadate);
# quit;
###!!! It is now different from SAS when querying credit rating, instead of using compr, we shoud use comp to access data
compr=conn.raw_sql("""
                  select splticrm, gvkey, datadate
                  from comp.adsprate
                  """)
compr.datadate=pd.to_datetime(compr.datadate)

ccm_jun['year']=ccm_jun['datadate'].dt.year
compr['year']=compr['datadate'].dt.year
ccm_jun=pd.merge(ccm_jun, compr[['splticrm','gvkey','year']], how='left', on=['gvkey','year'])

# proc sort data=data2 nodupkey;
# 	by gvkey datadate;
# run;
ccm_jun = ccm_jun.sort_values(['gvkey', 'datadate'])


"""
data cpi;
	infile datalines;
	input yr 4.0 cpi 10.3;
	datalines;
    [manually inputted data]
;
run;

proc sql;
	create table data2 as select a.*, b.cpi from data2 a left join cpi b on 
		a.fyear=b.yr;
quit;

proc sort data=data2 nodupkey;
	by gvkey datadate;
run;

"""
cpi=pd.DataFrame()
cpi['fyear'] = list(reversed(range(1924, 2017+1)))
cpi['cpi'] = [246.19,242.23,236.53,229.91,229.17,229.594,
              224.939,218.056,214.537,215.303,207.342,201.6,195.3,188.9,
              183.96,179.88,177.1,172.2,166.6,163,160.5,156.9,
              152.4,148.2,144.5,140.3,136.2,130.7,124,118.3,
              113.6,109.6,107.6,103.9,99.6,96.5,90.9,82.4,
              72.6,65.2,60.6,56.9,53.8,49.3,44.2,41.8,
              40.6,38.9,36.7,34.8,33.4,32.5,31.6,31,
              30.7,30.2,29.9,29.6,29.2,28.9,28.2,27.3,
              26.8,26.9,26.8,26.7,25.9,24,23.7,24.4,
              22.2,19.7,18.1,17.6,17.3,16.3,14.7,14,13.8,14.1,
              14.4,13.9,13.6,13.3,13.1,13.6,15.1,16.6,17.2,17.1,17.2,17.5,17.7,17]
ccm_jun = pd.merge(ccm_jun, cpi, how='left', on='fyear')
ccm_jun = ccm_jun.sort_values(['gvkey', 'datadate'])
ccm_jun = ccm_jun.drop_duplicates(['gvkey','datadate'])

"""
data data2;
	set data2;
	*an attempt to coding credit ratings into numerical format;
	by gvkey datadate;

	if splticrm='D' then
		credrat=1;

	if splticrm='C' then
		credrat=2;

	if splticrm='CC' then
		credrat=3;

	if splticrm='CCC-' then
		credrat=4;

	if splticrm='CCC' then
		credrat=5;

	if splticrm='CCC+' then
		credrat=6;

	if splticrm='B-' then
		credrat=7;

	if splticrm='B' then
		credrat=8;

	if splticrm='B+' then
		credrat=9;

	if splticrm='BB-' then
		credrat=10;

	if splticrm='BB' then
		credrat=11;

	if splticrm='BB+' then
		credrat=12;

	if splticrm='BBB-' then
		credrat=13;

	if splticrm='BBB' then
		credrat=14;

	if splticrm='BBB+' then
		credrat=15;

	if splticrm='A-' then
		credrat=16;

	if splticrm='A' then
		credrat=17;

	if splticrm='A+' then
		credrat=18;

	if splticrm='AA-' then
		credrat=19;

	if splticrm='AA' then
		credrat=20;

	if splticrm='AA+' then
		credrat=21;

	if splticrm='AAA' then
		credrat=22;
	*if missing(credrat) then credrat=0;

	if credrat<lag(credrat) then
		credrat_dwn=1;
	else
		credrat_dwn=0;

	if count=1 then
		credrat_dwn=0;
run;

proc sort data=data2 nodupkey;
	by gvkey datadate;
run;

"""
def lag(df, col, n=1, on='gvkey'):
    return df.groupby(on)[col].shift(n)
credit_mapping = {'D': 1, 'C': 2, 'CC': 3, 'CCC-':4, 'CCC':5, 'CCC+':6, 'B-': 7, 'B': 8, 'B+':9, 'BB-':10, 'BB':11, 'BB+': 12, 'BBB-':13, 
                  'BBB':14, 'BBB+': 15, 'A-':16, 'A':17, 'A+':18, 'AA-':19,'AA':20,'AA+':21,'AAA':22}
ccm_jun['credrat'] = ccm_jun['splticrm'].map(credit_mapping)
ccm_jun['credrat'] = ccm_jun['credrat']
ccm_jun['credrat'] = 0
ccm_jun.loc[ccm_jun['credrat'] < lag(ccm_jun, 'credrat'), 'credrat_dwn'] = 1
ccm_jun.loc[ccm_jun['credrat'] >= lag(ccm_jun, 'credrat'), 'credrat_dwn'] = 0

#note: credit ratings are different (database different)
#%%
"""

*finish orgcap measure;

data data(drop=permno);
	set data2;
	by gvkey datadate;
	retain orgcap_1;
	avgat=((at+lag(at))/2);

	if first.gvkey then
		orgcap_1=(xsga/cpi)/(.1+.15);
	else
		orgcap_1=orgcap_1*(1-.15)+xsga/cpi;
	orgcap=orgcap_1/avgat;

	if count=1 then
		orgcap=.;

	if first.gvkey then
		oc_1=xsga/(.1+.15);
	else
		oc_1=(1-.15)*lag(oc_1)+xsga/cpi;
	oca=oc_1/at;
	retain bc_1;

	if first.gvkey then
		bc_1=xad/(0.1+0.5);
	else
		do;

			if missing(xad) then
				xad=0;
			bc_1=0.5*bc_1+xad;
		end;
	bca=bc_1/at;
run;

"""


ccm_jun = ccm_jun.sort_values(['gvkey','datadate'])
ccm_jun['avgat'] = (ccm_jun['at']+ccm_jun['lagat'])/2
ccm_jun.loc[ccm_jun['count']==0, 'orgcap_1'] = (ccm_jun['xsga']/ccm_jun['cpi'])/(.1+.15)
orgcap_1 = ccm_jun[['orgcap_1','xsga','cpi']]
prev_row = None
for i, row in orgcap_1.iterrows():
    if(np.isnan(row['orgcap_1'])):
        row['orgcap_1'] = prev_row['orgcap_1']*(1-0.15)+row['xsga']/row['cpi']
    prev_row = row
ccm_jun['orgcap_1'] = orgcap_1['orgcap_1']

ccm_jun['orgcap'] = ccm_jun['orgcap_1']/ccm_jun['avgat']
ccm_jun.loc[ccm_jun['count'] == 0, 'orgcap'] = np.nan

ccm_jun.loc[ccm_jun['count']==0, 'oc_1'] = ccm_jun['xsga']/(.1+.15)
oc_1 = ccm_jun[['oc_1','xsga','cpi']]
prev_row = None
for i, row in oc_1.iterrows():
    if(np.isnan(row['oc_1'])):
        row['oc_1'] = (1-0.15)*prev_row['oc_1'] + row['xsga']/row['cpi']
    prev_row = row
ccm_jun['oc_1'] = oc_1['oc_1']

ccm_jun['oca'] = ccm_jun['oc_1']/ccm_jun['at']

""" #TODO
ccm_jun.loc[ccm_jun['count']==0, 'bc_1'] = ccm_jun['xad']/(0.1+0.5)
bc_1 = ccm_jun[['bc_1','xad']]
prev_row = None
for i, row in bc_1.iterrows():
    if(np.isnan(row['bc_1'])):
        if(np.isnan(row['xad'])):
            row['xad'] = 0
        row['bc_1'] = 0.5 * prev_row['bc_1'] + row['xad']
    prev_row = row
ccm_jun['bc_1'] = bc_1['bc_1']
ccm_jun['bca'] = ccm_jun['bc_1']/ccm_jun['at']
"""

"""
/*HXZ*/
proc sql;
	create table data as select *, (orgcap-mean(orgcap))/std(orgcap) as orgcap_ia, 
		(oca-mean(oca))/std(oca) as oca_ia from data group by sic2, fyear;
quit;

data data;
	set data;

	if dvpsx_f>0 then
		dvpsx_1=1;
	else
		dvpsx_1=0;
	ww=-0.091*(ib+dp)/at-0.062*dvpsx_1+0.021*dltt/at-0.044*log(at)+0.102*(indsale/lag(indsale)-1)-0.035*(sale/lag(sale)-1);

	if count=1 then
		ww=.;
run;

"""

mean_orgcap = ccm_jun.rename(columns={'orgcap':'orgcap_mean'}).groupby(['sic2','fyear'])['orgcap_mean'].mean()
std_orgcap = ccm_jun.rename(columns={'orgcap':'orgcap_std'}).groupby(['sic2','fyear'])['orgcap_std'].std()
ccm_jun = pd.merge(ccm_jun, mean_orgcap, on=['sic2','fyear'], how='left')
ccm_jun = pd.merge(ccm_jun, std_orgcap, on=['sic2','fyear'], how='left')
ccm_jun['orgcap_ia'] = (ccm_jun['orgcap']-ccm_jun['orgcap_mean'])/ccm_jun['orgcap_std']

mean_oca = ccm_jun.rename(columns={'oca':'oca_mean'}).groupby(['sic2','fyear'])['oca_mean'].mean()
std_oca = ccm_jun.rename(columns={'oca':'oca_std'}).groupby(['sic2','fyear'])['oca_std'].std()
ccm_jun = pd.merge(ccm_jun, mean_oca, on=['sic2','fyear'], how='left')
ccm_jun = pd.merge(ccm_jun, std_oca, on=['sic2','fyear'], how='left')
ccm_jun['oca_ia'] = (ccm_jun['oca']-ccm_jun['oca_mean'])/ccm_jun['oca_std']

ccm_jun.loc[ccm_jun['dvpsx_f']>0, 'dvpsx_1'] = 1
ccm_jun.loc[ccm_jun['dvpsx_f']<=0, 'dvpsx_1'] = 0
ccm_jun['ww'] = -0.091*(ccm_jun['ib']+ccm_jun['dp'])/ccm_jun['at'] - 0.062*ccm_jun['dvpsx_1'] + 0.021*ccm_jun['dltt']/ccm_jun['at'] -0.044*np.log(ccm_jun['at']) + 0.102*(ccm_jun['indsale']/lag(ccm_jun,'indsale')-1)-0.035*(ccm_jun['sale']/lag(ccm_jun,'sale')-1)
#%%
"""

/* check 2 */
*===========================================================

			Now moving past annual Compustat


	============================================================;
*===========================================================

			Create merge with CRSP


	============================================================;
*======================GET CRSP IDENTIFIER=============================;

proc sort data=cc.ccmxpf_linktable out=lnk;
	where LINKTYPE in ("LU", "LC", "LD", "LF", "LN", "LO", "LS", "LX") and
       					(2018 >=year(LINKDT) or LINKDT=.B) 
		and (1940 <=year(LINKENDDT) or LINKENDDT=.E);
	by GVKEY LINKDT;
run;

proc sql;
	create table temp as select a.lpermno as permno, b.* from lnk a, data b where 
		a.gvkey=b.gvkey and (LINKDT <=b.datadate or LINKDT=.B) 
		and (b.datadate <=LINKENDDT or LINKENDDT=.E) and lpermno ne . and not 
		missing(b.gvkey);
quit;

data temp;
	set temp;
	where not missing(permno);
run;
"""

lnk = conn.raw_sql("""
                     select * from crsp.ccmxpf_linktable
                     """)
lnk = lnk[lnk['linktype'].isin(['LU','LC','LD','LF','LN','LO','LS','LX'])]
lnk = lnk[(2018 >= lnk['linkdt'].astype(str).str[0:4].astype(int)) | (lnk['linkdt'] == '.B') ]
lnk = lnk[(lnk['linkenddt'].isna()) | ("1940" <= lnk['linkenddt'].astype(str).str[0:4])]
lnk = lnk.sort_values(['gvkey','linkdt'])
lnk['linkdt'] = pd.to_datetime(lnk['linkdt'])

ccm_jun2 = pd.merge(lnk[['gvkey','linkdt','linkenddt','lpermno']], ccm_jun, on='gvkey', how='inner')
ccm_jun2 = ccm_jun2[(ccm_jun2['linkdt'] <= ccm_jun2['datadate']) | (ccm_jun2['linkdt'] == '.B') ]
ccm_jun2 = ccm_jun2[(ccm_jun2['datadate'] <= ccm_jun2['linkenddt']) | (ccm_jun2['linkenddt'].isna())]
ccm_jun2 = ccm_jun2[(ccm_jun2['lpermno'] != '.') & ccm_jun2['gvkey'].notna()]

"""
*======================================

						Screen on Stock market information: common stocks and major exchanges


						=======================================;
*----------------------screen for only NYSE, AMEX, NASDAQ, and common stock-------------;

/* check 3 */
*==========================================================================================================


				Finalize first Compustat data set
				This is most of the annual compustat variables plus a couple components that still need additional information


==========================================================================================================;

data temp;
	set temp;
	keep gvkey permno datadate datadate_a fyear sic2 bm BE mve6b cfp ep cashpr dy 
		lev sp roic rd_sale rd_mve chadv agr invest gma chcsho lgr egr chpm chato 
		chinv hire cf acc pctacc absacc age spii spi sgr pchsale_pchinvt 
		pchsale_pchrect pchgm_pchsale pchsale_pchxsga pchcapx ps divi divo obklg 
		chobklg securedind secured convind conv grltnoa chdrc rd rdbias chpmia 
		chatoia chempia bm_ia pchcapx_ia tb cfp_ia mve_ia herf credrat credrat_dwn 
		orgcap m1-m6 grcapx depr pchdepr grGW tang woGW sin currat pchcurrat quick 
		pchquick

		/*HXZ*/
		orgcap_ia adm gad rdm rds ol rca bca etr lfe kz ww cdd roavol_a ala alm ob_a 
		cinvest_a dpia noa dnoa pchcapx3 cdi ivg dcoa dcol dwc dnca dncl dnco dfin ta 
		dsti dfnl poa nef ndf nxf atm cp op nop em ndp ebp rna pm ato cto gpa rmw ole 
		opa ola cop cla os zs bi oca oca_ia ha he pchsale_pchinvt_hxz 
		pchsale_pchrect_hxz pchgm_pchsale_hxz pchsale_pchxsga_hxz realestate_hxz 
		secured_hxz agr_hxz grltnoa_hxz chcsho_hxz pchcapx_ia_hxz acc_hxz egr_hxz 
		pctacc_hxz lev_hxz ep_hxz cfp_hxz tb_hxz salecash salerec saleinv pchsaleinv 
		cashdebt realestate roe operprof mve_f;
	where datadate>='01Jan2015'd;
run;
"""

temp = ccm_jun2[[ 'gvkey', 'permno', 'datadate', 'fyear', 'sic2', 'cfp', 'ep', 'cashpr', 'dy',
		'lev', 'sp', 'roic', 'rd_sale', 'chadv', 'agr', 'invest', 'gma', 'chcsho', 'lgr', 'egr', 'chpm', 'chato',
		'chinv', 'hire', 'cf', 'acc', 'pctacc', 'absacc', 'spii', 'spi', 'sgr', 'pchsale_pchinvt',
		'pchsale_pchrect', 'pchgm_pchsale', 'pchsale_pchxsga', 'pchcapx', 'ps', 'divi', 'divo', 'obklg',
		'chobklg', 'securedind', 'secured', 'convind', 'conv', 'grltnoa', 'chdrc', 'rd', 'rdbias', 'chpmia',
		'chatoia', 'chempia', 'pchcapx_ia', 'tb', 'cfp_ia', 'mve_ia', 'herf', 'credrat', 'credrat_dwn',
		'orgcap', 'grcapx', 'depr', 'pchdepr', 'grGW', 'tang', 'woGW', 'sin', 'currat', 'pchcurrat', 'quick',
		'pchquick', 'orgcap_ia', 'adm', 'gad', 'rdm', 'rds', 'ol', 'ww', 'cdd', 'roavol_a', 'ala', 'alm', 'ob_a',
		'cinvest_a', 'noa', 'dnoa', 'pchcapx3', 'cdi', 'ivg', 'dcoa', 'dcol', 'dwc', 'dnca', 'dncl', 'dnco', 'dfin', 'ta',
		'dsti', 'dfnl', 'poa', 'nef', 'ndf', 'atm', 'cp', 'op', 'nop', 'ndp', 'ebp', 'rna', 'pm', 'ato', 'cto', 'gpa', 'rmw', 'ole',
		'opa', 'ola', 'cop', 'cla', 'os', 'zs', 'bi', 'oca', 'oca_ia', 'ha', 'he', 'pchsale_pchinvt_hxz',
		'pchsale_pchrect_hxz', 'pchgm_pchsale_hxz', 'pchsale_pchxsga_hxz', 'realestate_hxz',
		'secured_hxz', 'agr_hxz', 'grltnoa_hxz', 'chcsho_hxz', 'pchcapx_ia_hxz', 'acc_hxz', 'egr_hxz',
		'pctacc_hxz', 'lev_hxz', 'ep_hxz', 'cfp_hxz', 'tb_hxz', 'salecash', 'salerec', 'pchsaleinv',
		'cashdebt', 'realestate', 'roe', 'operprof', 'mve_f', 'm1','m2','m3','m4','m5','m6']]

#currently missing some variables kz nxf bca
"""
temp = ccm_jun2[[ 'gvkey', 'permno', 'datadate', 'datadate_a', 'fyear', 'sic2', 'bm', 'BE', 'mve6b', 'cfp', 'ep', 'cashpr', 'dy',
		'lev', 'sp', 'roic', 'rd_sale', 'rd_mve', 'chadv', 'agr', 'invest', 'gma', 'chcsho', 'lgr', 'egr', 'chpm', 'chato',
		'chinv', 'hire', 'cf', 'acc', 'pctacc', 'absacc', 'age', 'spii', 'spi', 'sgr', 'pchsale_pchinvt',
		'pchsale_pchrect', 'pchgm_pchsale', 'pchsale_pchxsga', 'pchcapx', 'ps', 'divi', 'divo', 'obklg',
		'chobklg', 'securedind', 'secured', 'convind', 'conv', 'grltnoa', 'chdrc', 'rd', 'rdbias', 'chpmia',
		'chatoia', 'chempia', 'bm_ia', 'pchcapx_ia', 'tb', 'cfp_ia', 'mve_ia', 'herf', 'credrat', 'credrat_dwn',
		'orgcap', 'm1-m6', 'grcapx', 'depr', 'pchdepr', 'grGW', 'tang', 'woGW', 'sin', 'currat', 'pchcurrat', 'quick',
		'pchquick', 'orgcap_ia', 'adm', 'gad', 'rdm', 'rds', 'ol', 'rca', 'bca', 'etr', 'lfe', 'kz', 'ww', 'cdd', 'roavol_a', 'ala', 'alm', 'ob_a',
		'cinvest_a', 'dpia', 'noa', 'dnoa', 'pchcapx3', 'cdi', 'ivg', 'dcoa', 'dcol', 'dwc', 'dnca', 'dncl', 'dnco', 'dfin', 'ta',
		'dsti', 'dfnl', 'poa', 'nef', 'ndf', 'nxf', 'atm', 'cp', 'op', 'nop', 'em', 'ndp', 'ebp', 'rna', 'pm', 'ato', 'cto', 'gpa', 'rmw', 'ole',
		'opa', 'ola', 'cop', 'cla', 'os', 'zs', 'bi', 'oca', 'oca_ia', 'ha', 'he', 'pchsale_pchinvt_hxz',
		'pchsale_pchrect_hxz', 'pchgm_pchsale_hxz', 'pchsale_pchxsga_hxz', 'realestate_hxz',
		'secured_hxz', 'agr_hxz', 'grltnoa_hxz', 'chcsho_hxz', 'pchcapx_ia_hxz', 'acc_hxz', 'egr_hxz',
		'pctacc_hxz', 'lev_hxz', 'ep_hxz', 'cfp_hxz', 'tb_hxz', 'salecash', 'salerec', 'saleiRnv', 'pchsaleinv',
		'cashdebt', 'realestate', 'roe', 'operprof', 'mve_f']]
"""

#temp about equal to TEMP (line 1350) (21186 rows vs 21177 rows)

"""     

*========================================================================================================

		Now align the annual Compustat variables in calendar month with the assumption that
		annual information is available with a lag of 6 months (if we had point-in-time we would use that)

=========================================================================================================;
*---------------------------add returns and monthly CRSP data we need later-----------------------------;

proc sql;
	create table temp2 as select a.ret, /*HXZ*/
	a.retx, abs(prc) as prc, shrout, vol, a.date, a.permno as newpermno, b.* from 
		crsp.msf a left join temp b on a.permno=b.permno and intnx('MONTH', 
		b.datadate, 7)<=a.date<intnx('MONTH', b.datadate, 20);
quit;

data temp2(drop=permno);
	set temp2;
run;

data temp2(rename=(newpermno=permno));
	set temp2;
run;
"""
crsp_msf = conn.raw_sql("""
                      select ret, retx, prc, shrout, vol, date, permno from crsp.msf
                      where date >= '01/01/2015'
                      """) 
crsp_msf = crsp_msf[crsp_msf['permno'].isin(temp['permno'])]
crsp_msf['date'] = pd.to_datetime(crsp_msf['date'])
crsp_msf = crsp_msf.sort_values('date')

z = temp[['datadate','permno']]
z['date_l'] = temp['datadate'] + pd.TimedeltaIndex([7]*len(z), 'M') + pd.TimedeltaIndex([-5]*len(z), 'd')
z['date_u'] = temp['datadate'] + pd.TimedeltaIndex([20]*len(z), 'M')
z = pd.merge(z, crsp_msf, on='permno', how='left')
z['date'] = pd.to_datetime(z['date'])
z = z[(z['date'] >= z['date_l']) & (z['date'] < z['date_u'])]

temp2 = pd.merge(z, temp, on=['permno','datadate'], how='left')
#142,805 rows in Python vs 154,220 rows in SAS
#Difference due to differences in date cut-off (7 months evaluated differently between SAS and Python)

"""
proc sort data=crsp.mseall(keep=date permno exchcd shrcd siccd) out=mseall 
		nodupkey;
	where exchcd in (1, 2, 3) and shrcd in (10, 11);
	by permno exchcd date;
run;

proc sql;
	create table mseall as select *, min(date) as exchstdt, max(date) as exchedt 
		from mseall group by permno, exchcd;
quit;

proc sort data=mseall(rename=(date=time_1)) nodupkey;
	by permno exchcd;
run;

proc sql;
	create table temp2 as select * from temp2 as a left join mseall as b on 
		a.permno=b.permno and exchstdt<=a.date<=exchedt;
quit;

data temp2;
	set temp2;
	where exchcd in (1, 2, 3) and shrcd in (10, 11);
	drop exchstdt exchedt;
run;
"""
crsp_mseall = conn.raw_sql("""
                      select date, permno, exchcd, shrcd, siccd from crsp.mseall
                      where exchcd in (1, 2, 3) and shrcd in (10, 11)
                      """) 
crsp_mseall = crsp_mseall.sort_values(['permno','exchcd','date'])
mseall_min = crsp_mseall.groupby(['permno','exchcd'])['date'].min().reset_index().rename(columns={'date':'exchstdt'})
mseall_max = crsp_mseall.groupby(['permno','exchcd'])['date'].max().reset_index().rename(columns={'date':'exchedt'})

crsp_mseall = pd.merge(crsp_mseall, mseall_min, on=['permno','exchcd'])
crsp_mseall = pd.merge(crsp_mseall, mseall_max, on=['permno','exchcd'])
crsp_mseall = crsp_mseall.rename(columns={'date':'time_1'})
crsp_mseall = crsp_mseall.sort_values(['permno','exchcd'])
crsp_mseall = crsp_mseall.drop_duplicates(['permno','exchcd'])

temp2 = pd.merge(temp2, crsp_mseall, on='permno', how='left')
temp2 = temp2[((temp2['date']>=temp2['exchstdt']) & (temp2['date']<=temp2['exchedt']))] 

#106982 vs 114761

"""
*-----------Included delisted returns in the monthly returns--------------------;

proc sql;
	create table temp2 as select a.*, b.dlret, b.dlstcd, b.exchcd, b.shrcd, 
		b.siccd from temp2 a left join crsp.mseall b on a.permno=b.permno and 
		a.date=b.date;
quit;

data temp2;
	set temp2;

	if missing(dlret) and (dlstcd=500 or (dlstcd>=520 and dlstcd<=584)) and exchcd 
		in (1, 2) then
			dlret=-.35;

	if missing(dlret) and (dlstcd=500 or (dlstcd>=520 and dlstcd<=584)) and exchcd 
		in (3) then
			dlret=-.55;
	*see Johnson and Zhao (2007), Shumway and Warther (1999) etc.;

	if not missing(dlret) and dlret<-1 then
		dlret=-1;

	if missing(dlret) then
		dlret=0;
	ret=ret+dlret;

	if missing(ret) and dlret ne 0 then
		ret=dlret;
run;

proc sort data=temp2;
	by permno date descending datadate;
run;

proc sort data=temp2 nodupkey;
	by permno date;
run;
"""
crsp_mseall_dl = conn.raw_sql("""
                      select dlret, dlstcd, exchcd, shrcd, siccd, date, permno from crsp.mseall
                      where date >= '01/01/2015'
                      """) 
crsp_mseall_dl['date'] = pd.to_datetime(crsp_mseall_dl['date'])
temp2 = pd.merge(temp2, crsp_mseall_dl, on=['date','permno'])

temp2['exchcd'] = temp2['exchcd_x']

temp2.loc[ (temp2['dlret'].isna()) & ((temp2['dlstcd']==500) | ((temp2['dlstcd']>=520) & (temp2['dlstcd']<=584))) & (temp2['exchcd'].isin([1,2])), 'dlret'] = -0.35
temp2.loc[ (temp2['dlret'].isna()) & ((temp2['dlstcd']==500) | ((temp2['dlstcd']>=520) & (temp2['dlstcd']<=584))) & (temp2['exchcd'].isin([3])), 'dlret'] = -0.55
temp2.loc[ (temp2['dlret'].notna()) & (temp2['dlret']<-1), 'dlret'] = -1
temp2.loc[ (temp2['dlret'].isna()), 'dlret'] = 0 #TODO: wtf? this should not be 0... i think this should be not missing...
temp2['ret'] = temp2['ret'] + temp2['dlret']

temp2 = temp2.sort_values(['permno','date','datadate'], ascending=[True, True, False])
temp2 = temp2.drop_duplicates(['permno','date'])

#temp2 about equal to TEMP2 (+ filter for datadate ne .) (99,628 rows in Python, 107,695 rows in SAS)

"""

*can use monthly market cap and price, but need to lag because it is currently

							contemporaneous with the returns we want to predict;

data temp2(rename=(datadate=time_2));
	set temp2;
	by permno date;

	/*market cap measure*/
	*mve=log(mve_f);
	*we use mve_now to replace mve_m

								*we use mve_past to replace mve
								*mve=log(mve_m);
	mve0=abs(prc)*shrout;
	*mvel1=abs(lag(prc))*lag(shrout);
	mvel1=lag(mve0);
	pps=lag(prc);

	if first.permno then
		mvel1=.;

	if first.permno then
		pps=.;
run;
"""
temp2 = temp2.rename(columns={'datadate':'time_2'})
temp2['mve0'] = np.abs(temp2['prc'])*temp2['shrout']
temp2['mvel1'] = lag(temp2, 'mve0')
temp2['pps'] = lag(temp2, 'prc')

#%%

"""
*==============================================================================================================


				Now add in COMPUSTAT QUARTERLY and then add to the monthly aligned dataset


==============================================================================================================;
proc sql;
	create table data as select substr(compress(f.cusip), 1, 6) as cnum, c.gvkey, 
		fyearq, fqtr, datadate, datadate as datadate_q, rdq, cat(sic, 1, 2) as sic2, 
		/*income items*/
		ibq, saleq, txtq, revtq, cogsq, xsgaq, /*balance sheet items*/
		atq, actq, cheq, lctq, dlcq, ppentq, /*HXZ*/
		xrdq, rectq, invtq, ppegtq, txdbq, dlttq, dvpsxq, gdwlq, intanq, txditcq, 
		dpq, oibdpq, cshprq, ajexq, oiadpq, ivaoq, mibq, xintq, drcq, drltq, apq, 
		/*other*/
		abs(prccq) as prccq, abs(prccq)*cshoq as mveq, ceqq, seqq, pstkq, atq, ltq, 
		pstkrq from comp.names as c, comp.fundq as f where f.gvkey=c.gvkey and 
		f.indfmt='INDL' and f.datafmt='STD' and f.popsrc='D' and f.consol='C' and 
		datadate>='01Jan2015'd;
quit;
"""

comp_qtr = conn.raw_sql("""
                        select f.cusip as cnum, c.gvkey, fyearq, fqtr, datadate as datadate_q, rdq, sic as sic2,
                		ibq, saleq, txtq, revtq, cogsq, xsgaq,
                		atq, actq, cheq, lctq, dlcq, ppentq, 
                		xrdq, rectq, invtq, ppegtq, txdbq, dlttq, dvpsxq, gdwlq, intanq, txditcq, 
                		dpq, oibdpq, cshprq, ajexq, oiadpq, ivaoq, mibq, xintq, drcq, drltq, apq, 
                		abs(prccq) as prccq, abs(prccq)*cshoq as mveq, ceqq, seqq, pstkq, atq, ltq, pstkrq
                        from comp.names as c, comp.fundq as f
                        where
                        f.gvkey = c.gvkey
                        and f.indfmt='INDL'
                        and f.datafmt='STD'
                        and f.popsrc='D'
                        and f.consol='C'
                        and datadate >= '01/01/2015'
                        """)
                        
comp_qtr = comp_qtr.loc[:,~comp_qtr.columns.duplicated()]
comp_qtr['cshoq'] = comp_qtr['mveq'] / abs(comp_qtr['prccq'])
"""
proc sort data=data nodupkey;
	by gvkey datadate;
run;

proc sort data=data;
	by gvkey datadate;
run;
"""

comp_qtr = comp_qtr.sort_values(['gvkey','datadate_q'])
comp_qtr = comp_qtr.drop_duplicates(['gvkey', 'datadate_q'])

"""
data data3;
	set data;
	by gvkey datadate;
	retain count;

	if not missing(pstkrq) then
		pstk=pstkrq;
	else
		pstk=pstkq;
	scal=seqq;

	if missing(seqq) then
		scal=ceqq+pstk;

	if missing(seqq) and (missing(ceqq) or missing(pstk)) then
		scal=atq-ltq;
	chtx=(txtq-lag4(txtq))/lag4(atq);
	roaq=ibq/lag(atq);
	roeq=(ibq)/lag(scal);
	rsup=(saleq-lag4(saleq))/mveq;
	sacc=((actq-lag(actq) - (cheq-lag(cheq))) - ((lctq-lag(lctq))-(dlcq-lag(dlcq)) 
		) ) /saleq;
	;
    if saleq<=0 then
	sacc=((actq-lag(actq) - (cheq-lag(cheq))) 
		- ((lctq-lag(lctq))-(dlcq-lag(dlcq)) ) ) /.01;
"""

def lag(df, col, n=1, on='gvkey'):
    z = df.groupby(on)[col].shift(n)
    z = z.reset_index()
    z = z.sort_values('index')
    z = z.set_index('index')
    return z[col]

compq3 = comp_qtr
compq3.loc[compq3['pstkrq'].notna(), 'pstk'] = compq3['pstkrq']
compq3.loc[compq3['pstkrq'].isna(), 'ptsk'] = compq3['pstkq']
compq3['scal'] = compq3['seqq']
compq3.loc[compq3['seqq'].isna(), 'scal'] = compq3['ceqq'] + compq3['pstk']
compq3.loc[(compq3['seqq'].isna()) & ((compq3['ceqq'].isna()) | (compq3['pstk'].isna())),'scal'] = compq3['atq'] - compq3['ltq']

compq3['chtx'] = (compq3['txtq']-lag(compq3, 'txtq', 4))/lag(compq3, 'atq', 4)
compq3['roaq'] = compq3['ibq']/lag(compq3, 'atq')
compq3['roeq'] = compq3['ibq']/lag(compq3, 'scal')
compq3['rsup'] = (compq3['saleq']-lag(compq3, 'saleq', 4))/compq3['mveq']
compq3['sacc'] = ( ((compq3['actq']-lag(compq3, 'actq')) - (compq3['cheq'] - lag(compq3, 'cheq'))) 
                  -((compq3['lctq']-lag(compq3, 'lctq')) - (compq3['dlcq'] - lag(compq3, 'dlcq')))) / compq3['saleq']
compq3.loc[compq3['saleq'] <= 0, 'sacc'] = ( ((compq3['actq']-lag(compq3, 'actq')) - (compq3['cheq'] - lag(compq3, 'cheq'))) 
                                            -((compq3['lctq']-lag(compq3, 'lctq')) - (compq3['dlcq'] - lag(compq3, 'dlcq')))) / 0.01
"""
stdacc=std(sacc, lag(sacc), lag2(sacc), lag3(sacc), lag4(sacc), lag5(sacc), 
	lag6(sacc), lag7(sacc), lag8(sacc), lag9(sacc), lag10(sacc), lag11(sacc), 
	lag12(sacc), lag13(sacc), lag14(sacc), lag15(sacc));
sgrvol=std(rsup, lag(rsup), lag2(rsup), lag3(rsup), lag4(rsup), lag5(rsup), 
	lag6(rsup), lag7(rsup), lag8(rsup), lag9(rsup), lag10(rsup), lag11(rsup), 
	lag12(rsup), lag13(rsup), lag14(rsup));
roavol=std(roaq, lag(roaq), lag2(roaq), lag3(roaq), lag4(roaq), lag5(roaq), 
	lag6(roaq), lag7(roaq), lag8(roaq), lag9(roaq), lag10(roaq), lag11(roaq), 
	lag12(roaq), lag13(roaq), lag14(roaq), lag15(roaq));

scf=(ibq/saleq)-sacc;
if saleq<=0 then
	scf=(ibq/.01)-sacc;
stdcf=std(scf, lag(scf), lag2(scf), lag3(scf), lag4(scf), lag5(scf), 
	lag6(scf), lag7(scf), lag8(scf), lag9(scf), lag10(scf), lag11(scf), 
	lag12(scf), lag13(scf), lag14(scf), lag15(scf));
cash=cheq/atq;
cinvest=((ppentq-lag(ppentq))/saleq)-mean(((lag(ppentq)-lag2(ppentq))/lag(saleq)), 
	((lag2(ppentq)-lag3(ppentq))/lag2(saleq)), 
	((lag3(ppentq)-lag4(ppentq))/lag3(saleq)));
if saleq<=0 then
	cinvest=((ppentq-lag(ppentq))/.01)-mean(((lag(ppentq)-lag2(ppentq))/(.01)), 
		((lag2(ppentq)-lag3(ppentq))/(.01)), ((lag3(ppentq)-lag4(ppentq))/(.01)));
*for sue later and for nincr;
che=ibq-lag4(ibq);
nincr=((ibq>lag(ibq)) + (ibq>lag(ibq))*(lag(ibq)>lag2(ibq)) 
	+ (ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq)) + (ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq))*(lag3(ibq)>lag4(ibq)) + (ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq))*(lag3(ibq)>lag4(ibq))*(lag4(ibq)>lag5(ibq)) + (ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq))*(lag3(ibq)>lag4(ibq))*(lag4(ibq)>lag5(ibq))*(lag5(ibq)>lag6(ibq)) + (ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq))*(lag3(ibq)>lag4(ibq))*(lag4(ibq)>lag5(ibq))*(lag5(ibq)>lag6(ibq))*(lag6(ibq)>lag7(ibq)) + (ibq>lag(ibq))*(lag(ibq)>lag2(ibq))*(lag2(ibq)>lag3(ibq))*(lag3(ibq)>lag4(ibq))*(lag4(ibq)>lag5(ibq))*(lag5(ibq)>lag6(ibq))*(lag6(ibq)>lag7(ibq))*(lag7(ibq)>lag8(ibq)) 
	);
"""

def trailing_std(df, col, n=15, on='gvkey'):
    z = df.groupby(on)[col].rolling(n).std()
    z = z.reset_index()
    z = z.sort_values('level_1')
    z = z.set_index('level_1')
    return z[col]
compq3['stdacc'] = trailing_std(compq3, 'sacc', 16)
compq3['sgrvol'] = trailing_std(compq3, 'rsup', 15)
compq3['roavol'] = trailing_std(compq3, 'roaq', 15)
compq3['scf'] = compq3['ibq']/compq3['saleq'] - compq3['sacc']
compq3.loc[compq3['saleq']<=0, 'scf'] = compq3['ibq']/0.01 - compq3['sacc']
compq3['stdcf'] = trailing_std(compq3, 'scf', 16)
compq3['cash'] = compq3['cheq']/compq3['atq']
compq3['cinvest'] = (compq3['ppentq'] - lag(compq3, 'ppentq'))/compq3['saleq'] - (1/3)*((lag(compq3,'ppentq') - lag(compq3, 'ppentq',2))/lag(compq3, 'saleq')) - (1/3)*((lag(compq3, 'ppentq',2) - lag(compq3, 'ppentq',3))/lag(compq3, 'saleq', 2)) - (1/3)*((lag(compq3,'ppentq',3) - lag(compq3, 'ppentq',4))/lag(compq3, 'saleq',3))
compq3.loc[compq3['saleq'] <= 0, 'cinvest'] = (compq3['ppentq'] - lag(compq3, 'ppentq'))/0.01 - (1/3)*((lag(compq3,'ppentq') - lag(compq3, 'ppentq',2))/0.01) - (1/3)*((lag(compq3, 'ppentq',2) - lag(compq3, 'ppentq',3))/0.01) - (1/3)*((lag(compq3,'ppentq',3) - lag(compq3, 'ppentq',4))/0.01)
compq3['che'] = compq3['ibq'] = lag(compq3, 'ibq',4)
#compq3['nincr']
#TODO: nincr

"""
rdmq=xrdq/mveq;
rdsq=xrdq/saleq;
olq=(cogsq+xsgaq)/atq;
tanq=(cheq + .715*rectq + .54*invtq + .535*ppegtq) / atq;
kzq=-1.002*((ibq+lag(ibq)+lag2(ibq)+lag3(ibq)+dpq) / lag(ppentq) ) 
	+ .283*(atq+mveq-ceqq-txdbq) / atq
		- 3.139*(dlcq+dlttq) / (dlcq+dlttq+seqq) 
	+39.368*(dvpsxq*cshoq + lag(dvpsxq*cshoq) + lag2(dvpsxq*cshoq) 
	+ lag3(dvpsxq*cshoq)) / lag(ppentq) -1.315*cheq / lag(ppentq);

if missing(gdwlq) then
	gdwlq=0;

if missing(intanq) then
	intanq=0;
alaq=cheq + .75*(actq - cheq) + .5*(atq - actq - gdwlq - intanq);
almq=alaq/(atq + mveq - ceqq);
laq=atq/lag(atq) - 1;
"""

compq3['rdmq'] = compq3['xrdq']/compq3['mveq']
compq3['rdsq'] = compq3['xrdq']/compq3['saleq']
compq3['olq'] = (compq3['cogsq'] + compq3['xsgaq'])/compq3['atq']
compq3['tanq'] = (compq3['cheq'] + 0.715*compq3['rectq'] + 0.54*compq3['invtq'] + 0.535*compq3['ppegtq'])/compq3['atq']
compq3['kzq'] = -1.002*((compq3['ibq'] + lag(compq3, 'ibq', 1) + lag(compq3, 'ibq', 2) + lag(compq3, 'ibq', 3) + compq3['dpq']) / lag(compq3,'ppentq')) \
    + 0.283*(compq3['atq'] + compq3['mveq'] - compq3['ceqq'] - compq3['txdbq'])/compq3['atq'] \
    - 3.139 * (compq3['dlcq'] + compq3['dlttq'])/(compq3['dlcq'] + compq3['dlttq'] + compq3['seqq']) \
    + 39.368*(compq3['dvpsxq'] * compq3['cshoq'] + lag(compq3, 'dvpsxq')*lag(compq3, 'cshoq') + lag(compq3, 'dvpsxq',2)*lag(compq3, 'cshoq',2) + lag(compq3, 'dvpsxq',3)*lag(compq3, 'cshoq', 3))/lag(compq3, 'ppentq')\
    - 1.315*compq3['cheq']/lag(compq3, 'ppentq')
compq3.loc[compq3['gdwlq'].isna(), 'gdwlq'] = 0
compq3.loc[compq3['intanq'].isna(), 'intanq'] = 0
compq3['alaq'] = compq3['cheq'] + 0.75*(compq3['actq'] - compq3['cheq']) + 0.5*(compq3['atq'] - compq3['actq'] - compq3['gdwlq'] - compq3['intanq'])
compq3['almq'] = compq3['alaq']/(compq3['atq']+compq3['mveq'] - compq3['ceqq'])
compq3['laq'] = compq3['atq']/lag(compq3, 'atq') - 1

"""
if missing(seqq) then
	seqq=ceqq + pstkq -ltq;

if not missing(txditcq) then
	bmq=(seqq + txditcq - pstkq)/mveq;
dmq=(dlcq + dlttq)/mveq;
amq=atq/mveq;
epq=ibq/mveq;
cpq=(ibq + dpq)/mveq;
emq=(mveq + dlcq + dlttq + pstkq - cheq)/oibdpq;
spq=saleq/mveq;
ndpq=dlttq + dlcq + pstkq -cheq;
ebpq=(ndpq + ceqq)/(ndpq + mveq);
x_1=saleq/(cshprq*ajexq);
rs=(x_1-lag4(x_1))/std(x_1, lag(x_1), lag2(x_1), lag3(x_1), lag4(x_1), 
	lag5(x_1));
droeq=roeq-lag4(roeq);
droaq=roaq - lag4(roaq);
"""

compq3.loc[compq3['seqq'].isna(), 'seqq'] = compq3['ceqq'] + compq3['pstkq'] - compq3['ltq']
compq3.loc[compq3['txditcq'].notna(), 'bmq'] = (compq3['seqq'] + compq3['txditcq'] - compq3['pstkq'])/compq3['mveq']
compq3['dmq'] = (compq3['dlcq'] + compq3['dlttq'])/compq3['mveq']
compq3['amq'] = compq3['atq']/compq3['mveq']
compq3['epq'] = compq3['ibq']/compq3['mveq']
compq3['cpq'] = (compq3['ibq'] + compq3['dpq'])/compq3['mveq']
compq3['emq'] = (compq3['mveq'] + compq3['dlcq'] + compq3['dlttq'] + compq3['pstkq'] - compq3['cheq'])/compq3['oibdpq']
compq3['spq'] = compq3['saleq']/compq3['mveq']
compq3['ndpq'] = compq3['dlttq'] + compq3['dlcq'] + compq3['pstkq'] - compq3['cheq']
compq3['ebpq'] = (compq3['ndpq'] + compq3['ceqq']) / (compq3['ndpq'] + compq3['mveq'])
compq3['x_1'] = compq3['saleq']/(compq3['cshprq']*compq3['ajexq'])
compq3['rs'] = (compq3['x_1'] - lag(compq3, 'x_1', 4)) / trailing_std(compq3, 'x_1', 6)
compq3['droeq'] = compq3['roeq'] - lag(compq3, 'roeq', 4)
compq3['droaq'] = compq3['roaq'] - lag(compq3, 'roaq', 4)

"""
if missing(dlcq) then
	dlcq=0;

if missing(ivaoq) then
	ivaoq=0;

if missing(mibq) then
	mibq=0;

if missing(pstkq) then
	pstkq=0;
noaq=(atq-cheq-ivaoq) - (dlcq-dlttq-mibq-pstkq-ceqq)/lag(atq);
rnaq=oiadpq/lag(noaq);
pmq=oiadpq/saleq;
atoq=saleq/lag(noaq);
ctoq=saleq/lag(atq);
glaq=(revtq-cogsq)/lag(atq);
oleq=(revtq-cogsq-xsgaq-xintq)/lag(bmq);
olaq=(revtq-cogsq-xsgaq+xrdq)/lag(atq);
claq=(revtq-cogsq-xsgaq+xrdq-(rectq-lag(rectq))-(invtq-lag(invtq)) 
	+drcq-lag(drcq)+drltq-lag(drltq)+apq-lag(apq)) / lag(atq);
blq=atq/bmq;
sgq=saleq/lag4(saleq);
"""

compq3.loc[compq3['dlcq'].isna(), 'dlcq'] = 0
compq3.loc[compq3['ivaoq'].isna(), 'ivaoq'] = 0
compq3.loc[compq3['mibq'].isna(), 'mibq'] = 0
compq3.loc[compq3['pstkq'].isna(), 'pstkq'] = 0
compq3['noaq'] = (compq3['atq']-compq3['cheq']-compq3['ivaoq']) - (compq3['dlcq']-compq3['dlttq']-compq3['mibq']-compq3['pstkq']-compq3['ceqq'])/lag(compq3, 'atq')
compq3['rnaq'] = compq3['oiadpq']/lag(compq3, 'noaq')
compq3['pmq'] = compq3['oiadpq']/compq3['saleq']
compq3['atoq'] = compq3['saleq']/lag(compq3,'noaq')
compq3['ctoq'] = compq3['saleq']/lag(compq3,'atq')
compq3['glaq'] = (compq3['revtq'] - compq3['cogsq'])/lag(compq3, 'atq')
compq3['oleq'] = (compq3['revtq'] - compq3['cogsq'] - compq3['xsgaq'] - compq3['xintq'])/lag(compq3, 'bmq')
compq3['olaq'] = (compq3['revtq'] - compq3['cogsq'] - compq3['xsgaq'] + compq3['xrdq'])/lag(compq3, 'atq')
compq3['claq'] = ((compq3['revtq'] - compq3['cogsq'] - compq3['xsgaq'] + compq3['xrdq'] - (compq3['rectq']-lag(compq3, 'rectq')) - (compq3['invtq']-lag(compq3, 'invtq')) \
      + compq3['drcq'] - lag(compq3, 'drcq') + compq3['drltq'] - lag(compq3, 'drltq') + compq3['apq'] - lag(compq3, 'apq'))) / lag(compq3, 'atq')
compq3['blq'] = compq3['atq']/compq3['bmq']
compq3['sgq'] = compq3['saleq']/lag(compq3, 'saleq', 4)

"""
proc sql;
	create table data3 as select *, sum(saleq) as indsaleq from data3 group by 
		sic2, fyearq;
quit;

data data3;
	set data3;

	if dvpsxq>0 then
		dvpsxq_1=1;
	else
		dvpsxq_1=0;
	wwq=-0.091*(ibq+dpq)/atq - 0.062*dvpsxq_1 + 0.021*dlttq/atq - 0.044*log(atq) +
												0.102*(indsaleq/lag(indsaleq)-1) 
		- 0.035*(saleq/lag(saleq)-1);

	if count=1 then
		wwq=.;
run;
"""

compq3.loc[compq3['dvpsxq'] > 0, 'dvpsxq_1'] = 1
compq3.loc[compq3['dvpsxq'] <= 0, 'dvpsxq_1'] = 0
temp_indsaleq = compq3.groupby(['sic2','fyearq'])['saleq'].sum().reset_index()
temp_indsaleq = temp_indsaleq.rename(columns={'saleq':'indsaleq'})
compq3 = pd.merge(compq3, temp_indsaleq, on=['fyearq','sic2'])
compq3['wwq'] = -0.091*(compq3['ibq'] + compq3['dpq'])/compq3['atq'] - 0.062*compq3['dvpsxq_1'] + 0.021*compq3['dlttq']/compq3['atq'] \
    + 0.044 * np.log(compq3['atq']) + 0.102*(compq3['indsaleq']/lag(compq3, 'indsaleq') - 1) - 0.035*(compq3['saleq']/lag(compq3, 'saleq') - 1)

"""
proc sort data=data3;
	by fyearq fqtr sic2;
run;

proc univariate data=data3 noprint;
	by fyearq fqtr sic2;
	var roavol sgrvol;
	output out=indmd median=md_roavol md_sgrvol;
run;

proc sql;
	create table data3 as select * from data3 a left join indmd b on 
		a.fyearq=b.fyearq and a.fqtr=b.fqtr and a.sic2=b.sic2;
quit;

proc sort data=data3 nodupkey;
	by gvkey fyearq fqtr;
run;

data data3;
	set data3;

	if roavol<md_roavol then
		m7=1;
	else
		m7=0;

	if sgrvol<md_sgrvol then
		m8=1;
	elsessssssssssss
		m8=0;
run;
"""

temp_md_roavol = compq3.groupby(['fyearq','fqtr','sic2'])['roavol'].median().reset_index()
temp_md_roavol = temp_md_roavol.rename(columns={'roavol' : 'md_roavol'})
compq3 = pd.merge(compq3, temp_md_roavol, on=['fyearq','fqtr','sic2'])
temp_md_sgrvol = compq3.groupby(['fyearq','fqtr','sic2'])['sgrvol'].median().reset_index()
temp_md_sgrvol = temp_md_sgrvol.rename(columns={'sgrvol' : 'md_sgrvol'})
compq3 = pd.merge(compq3, temp_md_sgrvol, on=['fyearq','fqtr','sic2'])
compq3.loc[compq3['roavol'] < compq3['md_roavol'], 'm7'] = 1
compq3.loc[compq3['roavol'] >= compq3['md_roavol'], 'm7'] = 0
compq3.loc[compq3['sgrvol'] < compq3['md_sgrvol'], 'm8'] = 1
compq3.loc[compq3['sgrvol'] >= compq3['md_sgrvol'], 'm8'] = 0

"""
proc sql;
	create table ibessum as select ticker, cusip, fpedats, statpers, ANNDATS_ACT, 
		numest, ANNTIMS_ACT, medest, actual, stdev from ibes.statsum_epsus where 
		fpi='6'  /*1 is for annual forecasts, 6 is for quarterly*/
		and statpers<ANNDATS_ACT /*only keep summarized forecasts prior to earnings annoucement*/
		and measure='EPS' and not missing(medest) and not missing(fpedats) 
		and (fpedats-statpers)>=0;
quit;

proc sort data=ibessum;
	by cusip fpedats descending statpers;
run;

proc sort data=ibessum nodupkey;
	by cusip fpedats;
run;
"""

ibessum = conn.raw_sql("""
                        select ticker, cusip, fpedats, statpers, ANNDATS_ACT, 
                        numest, ANNTIMS_ACT, medest, actual, stdev 
                        from ibes.statsum_epsus
                        where
                        fpi='6'
                        and statpers<ANNDATS_ACT
                        and measure='EPS' 
                        and (fpedats-statpers)>=0
                        """)
ibessum = ibessum[(ibessum['medest'].notna()) & (ibessum['fpedats'].notna())]
ibessum = ibessum.sort_values(by=['cusip','fpedats','statpers'], ascending=[True,True,False])
ibessum = ibessum.drop_duplicates(['cusip', 'fpedats'])       

"""
** Prepare Compustat-IBES translation file;

proc sort data=crsp.msenames(where=(ncusip ne '')) out=names nodupkey;
	by permno ncusip;
run;

* Add current cusip to IBES (IBES cusip is historical);

proc sql;
	create table ibessum2 as select a.*, substr(compress(b.cusip), 1, 6) as cusip6 
		from ibessum a left join names b on
						(a.cusip=b.ncusip);
quit;

* Merge IBES, CRSP/Compustat;

proc sql;
	create table data4 as select a.*, b.medest, b.actual from data3 a left join 
		ibessum2 b on
					(a.cnum=b.cusip6) and a.datadate=b.fpedats;
quit;

proc sort data=data4 nodupkey;
	by gvkey datadate;
run;

data data4;
	set data4;
	* finish SUE variable;

	if missing(medest) or missing(actual) then
		sue=che/mveq;

	if not missing(medest) and not missing(actual) then
		sue=(actual-medest)/abs(prccq);
run;
"""

crsp_msenames = conn.raw_sql("""select * from crsp.msenames""")
crsp_msenames = crsp_msenames[crsp_msenames['ncusip'].notna()]
crsp_msenames = crsp_msenames.sort_values(['permno','ncusip'])
crsp_msenames = crsp_msenames.drop_duplicates(['permno','ncusip'])
names = crsp_msenames.rename(columns={'cusip':'cusip6'})

ibessum2 = pd.merge(ibessum, names[['ncusip','cusip6']], left_on='cusip', right_on=['ncusip'], how='left')
ibessum2['cusip6'] = ibessum2['cusip6'].astype(str).str[0:6]

compq3['cnum'] = compq3['cnum'].astype(str).str[0:6]
compq4 = pd.merge(compq3, ibessum2[['medest','actual','cusip6','fpedats']], left_on=['cnum','datadate_q'], right_on=['cusip6','fpedats'], how='left')
compq4 = compq4.sort_values(['gvkey','datadate_q'])
compq4 = compq4.drop_duplicates(['gvkey','datadate_q'])

compq4.loc[(compq4['medest'].isna()) | (compq4['actual']).isna(), 'sue'] = compq4['che']/compq4['mveq']
compq4.loc[(compq4['medest'].notna()) & (compq4['actual']).notna(), 'sue'] = (compq4['actual'] - compq4['medest'])/abs(compq4['prccq'])

"""
proc sort data=cc.ccmxpf_linktable out=lnk;
	where LINKTYPE in ("LU", "LC", "LD", "LF", "LN", "LO", "LS", "LX") and
       					(2018 >=year(LINKDT) or LINKDT=.B) 
		and (1940 <=year(LINKENDDT) or LINKENDDT=.E);
	by GVKEY LINKDT;
run;

...

proc sql;
	create table data5 as select a.lpermno as permno, b.* from lnk a, data4 b 
		where a.gvkey=b.gvkey and (LINKDT <=b.datadate or LINKDT=.B) 
		and (b.datadate <=LINKENDDT or LINKENDDT=.E) and lpermno ne . and not 
		missing(b.gvkey);
quit;

data data5;
	set data5;
	where not missing(permno);
run;

data data5;
	set data5;
	where not missing(rdq);
	*seems like a reasonable screen at this point to make sure have at least some of this information;
run;
"""

lnk = conn.raw_sql("""
                     select * from crsp.ccmxpf_linktable
                     """)
lnk = lnk[lnk['linktype'].isin(['LU','LC','LD','LF','LN','LO','LS','LX'])]
lnk = lnk[(2018 >= lnk['linkdt'].astype(str).str[0:4].astype(int)) | (lnk['linkdt'] == '.B') ]
lnk = lnk[(lnk['linkenddt'].isna()) | ("1940" <= lnk['linkenddt'].astype(str).str[0:4])]
lnk = lnk.sort_values(['gvkey','linkdt'])

compq5 = pd.merge(lnk[['gvkey','linkdt','linkenddt','lpermno']], compq4, on='gvkey', how='inner')
compq5 = compq5[(compq5['linkdt'] <= compq5['datadate_q']) | (compq5['linkdt'] == '.B') ]
compq5 = compq5[(compq5['datadate_q'] <= compq5['linkenddt']) | (compq5['linkenddt'].isna())]
compq5 = compq5[(compq5['lpermno'] != '.') & compq5['gvkey'].notna()]
compq5 = compq5[(compq5['lpermno'].notna()) & (compq5['rdq'].notna())]

"""
proc sql;
	create table data5 as select a.*, b.vol from data5 a left join crsp.dsf b on 
		a.permno=b.permno and intnx('WEEKDAY', rdq, -30)<=b.date<=intnx('WEEKDAY', 
		rdq, -10);
quit;

proc sql;
	create table data5 as select *, mean(vol) as avgvol from data5 group by 
		permno, datadate, rdq;
quit;

proc sort data=data5(drop=vol) nodupkey;
	where not missing(rdq);
	by permno datadate rdq;
run;

proc sql;
	create table data6 as select a.*, b.vol, b.ret from data5 a left join crsp.dsf 
		b on a.permno=b.permno and intnx('WEEKDAY', rdq, 
		-1)<=b.date<=intnx('WEEKDAY', rdq, 1);
quit;

proc sql;
	create table data6 as select *, (mean(vol)-avgvol)/avgvol as aeavol, sum(ret) 
		as ear from data6 group by permno, datadate, rdq;
quit;

proc sort data=data6(drop=vol avgvol ret) nodupkey;
	by permno datadate rdq;
run;
"""

crsp_dsf = conn.raw_sql("""
                      select vol, ret, permno, date from crsp.dsf
                      where date >= '01/01/2015'
                      """) 
crsp_dsf = crsp_dsf[crsp_dsf['permno'].isin(compq5['lpermno'])]
crsp_dsf['date'] = pd.to_datetime(crsp_dsf['date'])
crsp_dsf = crsp_dsf.sort_values('date')

compq5['temp_rdq'] = np.busday_offset(compq5['rdq'].values.astype('datetime64[D]'), -10, roll='forward')
crsp_dsf['avgvol'] = trailing_std(crsp_dsf, 'vol', n=21, on='permno')
compq5 = pd.merge(compq5, crsp_dsf[['date', 'permno','avgvol']], how='left', left_on=['temp_rdq','lpermno'], right_on=['date', 'permno'])

compq5['temp_rdq'] = np.busday_offset(compq5['rdq'].values.astype('datetime64[D]'), 1, roll='forward')
crsp_dsf['aeavol'] = trailing_std(crsp_dsf, 'vol', n=3, on='permno')
crsp_dsf['ear'] = lag(crsp_dsf, 'ret', 0, on='permno') + lag(crsp_dsf, 'ret', 1, on='permno') + lag(crsp_dsf, 'ret', 2, on='permno')
compq6 = pd.merge(compq5, crsp_dsf[['date', 'permno', 'aeavol', 'avgvol', 'ear']], how='left', left_on=['temp_rdq','lpermno'], right_on=['date', 'permno'])
compq6['avgvol'] = compq6['avgvol_x'] #TODO
compq6['aeavol'] = (compq6['aeavol'] - compq6['avgvol'])/compq6['avgvol']

"""
data data6;
	set data6;
	keep gvkey permno datadate datadate_q rdq chtx roaq rsup stdacc stdcf sgrvol 
		roavol cash cinvest nincr
		rdmq rdsq olq tanq kzq alaq almq laq bmq dmq amq epq cpq emq spq ndpq ebpq 
		wwq rs droeq droaq noaq rnaq pmq atoq ctoq glaq oleq olaq claq blq sgq sue 
		aeavol ear m7 m8 prccq roeq;
	where datadate>='01Jan2015'd;
run;
"""
compq6 = compq6[['gvkey', 'lpermno', 'datadate_q', 'rdq', 'chtx', 'roaq', 'rsup', 'stdacc', 'stdcf', 'sgrvol', 
		'rdmq', 'rdsq', 'olq', 'tanq', 'kzq', 'alaq', 'almq', 'laq', 'bmq', 'dmq', 'amq', 'epq', 'cpq', 'emq', 'spq', 'ndpq', 'ebpq', 
		'wwq', 'rs', 'droeq', 'droaq', 'noaq', 'rnaq', 'pmq', 'atoq', 'ctoq', 'glaq', 'oleq', 'olaq', 'claq', 'blq', 'sgq', 'sue', 
		'roavol', 'cash', 'cinvest', 'm7', 'm8', 'prccq', 'roeq', 'aeavol', 'ear']]
compq6 = compq6.drop_duplicates()
#Add nincr
"""
compq6 = compq6[['gvkey', 'lpermno', 'datadate_q', 'rdq', 'chtx', 'roaq', 'rsup', 'stdacc', 'stdcf', 'sgrvol', 
		'rdmq', 'rdsq', 'olq', 'tanq', 'kzq', 'alaq', 'almq', 'laq', 'bmq', 'dmq', 'amq', 'epq', 'cpq', 'emq', 'spq', 'ndpq', 'ebpq', 
		'wwq', 'rs', 'droeq', 'droaq', 'noaq', 'rnaq', 'pmq', 'atoq', 'ctoq', 'glaq', 'oleq', 'olaq', 'claq', 'blq', 'sgq', 'sue', 
		'roavol', 'cash', 'cinvest', 'nincr', 'm7', 'm8', 'prccq', 'roeq', 'aeavol', 'ear']]
"""

#compq6 same as data6 in SAS (93,184 rows vs 93,224 rows)

#%%
"""
proc sql;
	create table temp3 as select * from temp2 a left join data6 b on 
		a.permno=b.permno and intnx('MONTH', a.date, -10)<=b.datadate<=intnx('MONTH', 
		a.date, -5) and intnx('WEEKDAY', b.rdq, 1)<intnx('MONTH', a.date, -5, 'beg');
	*allow at least four months for quarterly info to become available;
	*date is the end of the return month;
quit;

proc sort data=temp3;
	by permno date descending datadate;
run;

proc sort data=temp3 nodupkey;
	by permno date;
run;
"""

z = pd.merge(temp2, compq6.rename(columns={'lpermno':'permno'}), how='left', on='permno')

z['date_l'] = z['date'] + pd.TimedeltaIndex([-10]*len(z), 'M')
z['date_u'] = z['date'] + pd.TimedeltaIndex([-5]*len(z), 'M')
z['datadate_q'] = pd.to_datetime(z['datadate_q'])
z = z[(z['datadate_q'] >= z['date_l']) & (z['datadate_q'] <= z['date_u'])]

z['rdq_date_u'] = z['date'] + pd.TimedeltaIndex([-5]*len(z), 'M')
z['rdq'] = pd.to_datetime(z['rdq'])
temp3 = z[np.busday_offset(z['rdq'].values.astype('datetime64[D]'), 1, roll='forward') < z['rdq_date_u']]

temp3 = temp3.sort_values(['permno','date','datadate_q'], ascending=[True,True,False])
temp3 = temp3.drop_duplicates(['permno','date'])

# temp3 98,613 rows vs TEMP3 118,631 rows (filtered for datadate ne .)
# rdq is different between python and SAS

#%%
"""
*since some of the information is also using rdq, need to make sure that the monthly info is after that;

*----------------add eamonth--------------------------;

proc sort data=data6 out=lst(keep=permno rdq) nodupkey;
	where not missing(permno) and not missing(rdq);
	by permno rdq;
run;

proc sql;
	alter table lst add eamonth integer;
	update lst set eamonth=1;
	create table temp3 as select a.*, b.eamonth from temp3 a left join lst b on 
		a.permno=b.permno and year(a.date)=year(b.rdq) and month(a.date)=month(b.rdq);
	update temp3 set eamonth=0 where eamonth=.;
quit;

*finally finish Mohanram score;

data temp3;
	set temp3;
	ms=m1+m2+m3+m4+m5+m6+m7+m8;
	drop m1-m8;
run;

data temp3;
	set temp3;
	where siccd in (7000:9999);
run;
"""

lst = compq6[(compq6['lpermno'].notna()) & (compq6['rdq'].notna())]
lst = lst.sort_values(['lpermno','rdq'])
lst = lst.drop_duplicates(['lpermno','rdq'])
lst = lst[['lpermno','rdq']]

lst['eamonth'] = 1
temp3['date'] = pd.to_datetime(temp3['date'])
lst['rdq'] = pd.to_datetime(lst['rdq'])

temp3['month'] = temp3['date'].dt.month
temp3['year'] = temp3['date'].dt.year
lst['month'] = lst['rdq'].dt.month
lst['year'] = lst['rdq'].dt.year
temp3 = pd.merge(temp3, lst, how='left', left_on=['permno', 'month', 'year'], right_on=['lpermno', 'month', 'year'])

temp3['ms'] = temp3['m1'] + temp3['m2'] + temp3['m3'] + temp3['m4'] + temp3['m5'] + temp3['m6'] + temp3['m7'] + temp3['m8']

temp3 = temp3[(temp3['siccd'] >= 7000) & (temp3['siccd'] <= 9999)]

#temp3 python 29,653, sas 32,841

#%%
"""
/* check 6 */
*=================================================================================================================


				now add RPS that come straight from IBES data:
				set these up in monthly intervals where the IBES variables have the monthly statistical summary


=================================================================================================================;

proc sql;
	create table ibessum as select ticker, cusip, fpedats, statpers, ANNDATS_ACT, 
		numest, ANNTIMS_ACT, medest, meanest, actual, stdev from ibes.statsum_epsus 
		where fpi='1'  /*1 is for annual forecasts, 6 is for quarterly*/
		and statpers<ANNDATS_ACT /*only keep summarized forecasts prior to earnings annoucement*/
		and measure='EPS' and not missing(medest) and not missing(fpedats) 
		and (fpedats-statpers)>=0;
quit;

proc sort data=ibessum;
	by ticker cusip statpers descending fpedats;
	*doing this places all of the missing fpedats at the beginning of the file if not there....;
run;

proc sort data=ibessum nodupkey;
	by ticker cusip statpers;
run;

data ibessum;
	set ibessum;
	by ticker cusip statpers;
	disp=stdev/abs(meanest);

	if meanest=0 then
		disp=stdev/.01;
	*chfeps=meanest-lag(meanest);
	*if first.cusip then chfeps=.;
	chfeps=.;

	/* we do not use this characteristic in the current version */
run;
"""

ibessum = conn.raw_sql("""
                      select ticker, cusip, fpedats, statpers, ANNDATS_ACT, 
              		numest, ANNTIMS_ACT, medest, meanest, actual, stdev from ibes.statsum_epsus 
              		where fpi='1' 
              		and statpers<ANNDATS_ACT 
              		and measure='EPS'
              		and (fpedats-statpers)>=0;
                      """) 
ibessum = ibessum[(ibessum['medest'].notna()) & (ibessum['fpedats'].notna())]
ibessum = ibessum.sort_values(['ticker','cusip','statpers','fpedats'], ascending=[True,True,True,False])
ibessum = ibessum.sort_values(['ticker','cusip','statpers'])
ibessum = ibessum.drop_duplicates(['ticker','cusip','statpers'])

ibessum['disp'] = ibessum['stdev']/np.abs(ibessum['meanest'])
ibessum.loc[ibessum['meanest']==0, 'disp'] = ibessum['stdev']/0.01
ibessum['chfeps'] = np.nan
#exactly the same number of rows as ibessum in sas (1445522 rows)
"""
*add long term forecasts;

proc sql;
	create table ibessum2 as select ticker, cusip, fpedats, statpers, ANNDATS_ACT, 
		numest, ANNTIMS_ACT, medest, meanest, actual, stdev from ibes.statsum_epsus 
		where fpi='0'  /*1 is for annual forecasts, 6 is for quarterly,0 LTG*/
		/*and statpers<ANNDATS_ACT
		and measure='EPS'
		and not missing(medest)
		and (fpedats-statpers)>=0*/
		and not missing(meanest);
quit;

proc sort data=ibessum2 nodupkey;
	by ticker cusip statpers;
run;

proc sql;
	create table ibessum2b as select a.*, b.meanest as fgr5yr from ibessum a left 
		join ibessum2 b on a.ticker=b.ticker and a.cusip=b.cusip and 
		a.statpers=b.statpers;
quit;

"""
ibessum2 = conn.raw_sql("""
                      select ticker, cusip, fpedats, statpers, ANNDATS_ACT, 
        		numest, ANNTIMS_ACT, medest, meanest, actual, stdev from ibes.statsum_epsus 
        		where fpi='0'
                      """) 
ibessum2 = ibessum2[(ibessum2['medest'].notna()) & (ibessum2['meanest'].notna())]
ibessum2 = ibessum2.sort_values(['cusip','statpers'])
ibessum2 = ibessum2.drop_duplicates(['cusip','statpers'])
ibessum2 = ibessum2.rename(columns={'meanest':'fgr5yr'})

ibessum2b = pd.merge(ibessum, ibessum2[['fgr5yr', 'ticker','cusip','statpers']], how='left', on=['ticker','statpers','cusip'])
#exactly the same number as in sas (1445522 rows)
"""

data rec;
	set ibes.recdsum;
	where not missing(statpers) and not missing(meanrec);
run;

proc sql;
	create table ibessum2b as select a.*, b.meanrec from ibessum2b a left join rec 
		b on a.ticker=b.ticker and a.cusip=b.cusip and a.statpers=b.statpers;
quit;

proc sort data=ibessum2b;
	by ticker statpers;
run;

"""
rec = conn.raw_sql("""
                      select * from ibes.recdsum
                      """) 
rec = rec[(rec['statpers'].notna()) & (rec['meanrec'].notna())]

ibessum2b = pd.merge(ibessum2b, rec[['meanrec','ticker','cusip','statpers']], how='left', on=['ticker','cusip','statpers'])
ibessum2b = ibessum2b.sort_values(['ticker','statpers'])

"""
data ibessum2c;
	set ibessum2b;
	by ticker statpers;
	retain count;
	chrec=meanrec-mean(lag(meanrec), lag2(meanrec))-mean(lag3(meanrec), 
		lag4(meanrec), lag5(meanrec));

	if first.ticker then
		count=1;
	else
		count+1;

	if count<5 /*count<6*/
	then
		chrec=.;
run;

*prepare for merge;

proc sort data=crsp.msenames(where=(ncusip ne '')) out=names nodupkey;
	by permno ncusip;
run;

proc sql;
	create table ibessum2b as select a.*, b.permno from ibessum2c a left join 
		names b on
				(a.cusip=b.ncusip);
quit;
"""
ibessum2c = ibessum2b.copy()
ibessum2c['chrec'] = ibessum2b['meanrec'] - (1/2)*lag(ibessum2b, 'meanrec', 1, 'ticker') - (1/2)*lag(ibessum2b, 'meanrec', 2, 'ticker') \
    - (1/3)*lag(ibessum2b, 'meanrec', 3, 'ticker') - (1/3)*lag(ibessum2b, 'meanrec', 4, 'ticker') - (1/3)*lag(ibessum2b, 'meanrec', 5, 'ticker')

names = conn.raw_sql("""
                      select * from crsp.msenames
                      where ncusip != ''
                      """) 
names = names.sort_values(['permno','ncusip'])
names = names.drop_duplicates(['permno','ncusip'])

ibessum2b = pd.merge(ibessum2c, names[['permno','ncusip']], how='left', left_on=['cusip'], right_on=['ncusip'])

"""

*================================================================

Add IBES data to the rest of the data: note that these are not available at the beginning of sample, and the recommendations come even later
disp, chfeps, meanest, numest, sfe, fgr5yr--ltg --1989
meanrec, chrec --1994

	=================================================================;

proc sql;
	create table temp4 as select a.*, b.disp, b.chfeps, b.fgr5yr, b.statpers, 
		b.meanrec, b.chrec, b.numest as nanalyst, b.meanest/abs(a.prccq) as sfe, 
		b.meanest from temp3 a left join ibessum2b b on a.permno=b.permno and 
		intnx('MONTH', a.date, -4, 'beg')<=b.statpers<=intnx('MONTH', a.date, -1, 
		'end');
quit;[]

proc sort data=temp4;
	by permno date descending statpers;
run;

proc sort data=temp4 nodupkey;
	by permno date;
run;
"""
temp4 = pd.merge(temp3, ibessum2b, how='left', on='permno')
temp4['sfe'] = temp4['meanest']/np.abs(temp4['prccq'])
temp4['statpers'] = pd.to_datetime(temp4['statpers'])
temp4 = temp4[temp4['statpers'].isna() | (temp4['statpers'] >= (temp4['date'] + pd.TimedeltaIndex([-4]*len(temp4), 'M')) + MonthBegin(-1))]
temp4 = temp4[temp4['statpers'].isna() | (temp4['statpers'] <= (temp4['date'] + pd.TimedeltaIndex([-1]*len(temp4), 'M')) + MonthEnd(0))]

temp4 = temp4.sort_values(['permno', 'date', 'statpers'], ascending=[True, True, False])
temp4 = temp4.drop_duplicates(['permno','date'])
temp4 = temp4.rename(columns={'numest':'nanalyst'})
#temp4 has 26,391 & TEMP4 has 32,834
#something is wrong with statpers...

"""

*--------a little clean up for IBES variables------------;

data temp4;
	set temp4;

	if year(date)>=1989 and missing(nanalyst) then
		nanalyst=0;

	if year(date)>=1989 and missing(fgr5yr) then
		ltg=0;

	if year(date)>=1989 and not missing(fgr5yr) then
		ltg=1;
	array f{*} disp chfeps meanest nanalyst sfe ltg fgr5yr;
	array s{*} meanrec chrec;

	do i=1 to dim(f);

		if year(date)<1989 then
			f(i)=.;
	end;

	do j=1 to dim(s);

		if year(date)<1994 then
			s(j)=.;
	end;
run;

"""
temp4['year'] = temp4['date'].dt.year
temp4.loc[(temp4['year'] >= 1989) & (temp4['nanalyst'].isna()), 'nanalyst'] = 0
temp4.loc[(temp4['year'] >= 1989) & (temp4['fgr5yr'].isna()), 'ltg'] = 0
temp4.loc[(temp4['year'] >= 1989) & (temp4['fgr5yr'].notna()), 'ltg'] = 1

temp4.loc[(temp4['year'] < 1989), ['disp','chfeps','meanest','nanalyst','sfe','ltg','fgr5yr']] = np.nan
temp4.loc[(temp4['year'] < 1994), ['meanrec','chrec']] = np.nan

"""
/* check 7 */
*======================================================================================================

				There are some other variables that are based on monthly-CRSP information (already in the dataset from monthly CRSP)
				create those variables plus a couple of others

======================================================================================================;

data temp4;
	set temp4;
	*count to make sure we have enough time series for each firm to create variables;
	where not missing(ret);
	by permno date;
	retain count;

	if first.permno then
		count=1;
	else
		count+1;
run;

proc sql;
	create table temp4 as select *, mean(ret) as ewret  /*we have used this before, doesn't seem to make a big difference in the variables*/
	from temp4 group by date;
quit;

proc sort data=temp4;
	by permno date;
run;

data temp5;
	set temp4;
	where not missing(ret);
	by permno date;
	retain count;

	if first.permno then
		count=1;
	else
		count+1;
run;
"""

ewret = temp4.groupby('date')['ret'].mean().reset_index()
temp4 = pd.merge(temp4, ewret, how='left', on=['date'])
temp4 = temp4.sort_values(['permno','date'])
temp4['count']=temp4.groupby(['permno']).cumcount()

"""
data temp6;
	set temp5;
	chnanalyst=nanalyst-lag3(nanalyst);
	mom6m=((1+lag2(ret))*(1+lag3(ret))*(1+lag4(ret))*(1+lag5(ret))*(1+lag6(ret)) ) 
		- 1;
	mom12m=((1+lag2(ret))*(1+lag3(ret))*(1+lag4(ret))*(1+lag5(ret))*(1+lag6(ret))* (1+lag7(ret))*(1+lag8(ret))*(1+lag9(ret))*(1+lag10(ret))*(1+lag11(ret))*(1+lag12(ret)) 
		) - 1;
	mom36m=((1+lag13(ret))*(1+lag14(ret))*(1+lag15(ret))*(1+lag16(ret))*(1+lag17(ret))*(1+lag18(ret)) 
		* (1+lag19(ret))*(1+lag20(ret))*(1+lag21(ret))*(1+lag22(ret))*(1+lag23(ret))*(1+lag24(ret))* (1+lag25(ret))*(1+lag26(ret))*(1+lag27(ret))*(1+lag28(ret))*(1+lag29(ret))*(1+lag30(ret)) 
		* (1+lag31(ret))*(1+lag32(ret))*(1+lag33(ret))*(1+lag34(ret))*(1+lag35(ret))*(1+lag36(ret)) 
		) - 1;
	mom1m=lag(ret);
	dolvol=log(lag2(vol)*lag2(prc));
	chmom=((1+lag(ret))*(1+lag2(ret))*(1+lag3(ret))*(1+lag4(ret))*(1+lag5(ret))*(1+lag6(ret)) 
		) - 1
			- (((1+lag7(ret))*(1+lag8(ret))*(1+lag9(ret))*(1+lag10(ret))*(1+lag11(ret))*(1+lag12(ret)) 
		) - 1);
	turn=mean(lag(vol), lag2(vol), lag3(vol))/shrout;

	if lag(ret)>0 and lag2(ret)>0 and lag3(ret)>0 and lag4(ret)>0 and lag5(ret)>0 
		and lag6(ret)>0 then
			retcons_pos=1;
	else
		retcons_pos=0;

	if lag(ret)<0 and lag2(ret)<0 and lag3(ret)<0 and lag4(ret)<0 and lag5(ret)<0 
		and lag6(ret)<0 then
			retcons_neg=1;
	else
		retcons_neg=0;
"""
def lag(df, col, n=1, on='permno'):
    z = df.groupby(on)[col].shift(n)
    z = z.reset_index()
    z = z.sort_values('index')
    z = z.set_index('index')
    return z[col]

temp6 = temp4.copy()
temp6['chnanalyst'] = temp6['nanalyst'] - lag(temp6, 'nanalyst', 3)
temp6['mom6m'] = ((1+lag(temp6,'ret',2)) * (1+lag(temp6, 'ret', 3)) * (1+lag(temp6,'ret',4)) * (1+lag(temp6, 'ret', 5)) * (1+lag(temp6,'ret',6))) -1
temp6['mom12m'] = ((1+lag(temp6,'ret',2)) * (1+lag(temp6, 'ret', 3)) * (1+lag(temp6,'ret',4)) * (1+lag(temp6, 'ret', 5)) * (1+lag(temp6,'ret',6)) * \
    (1+lag(temp6,'ret',7))*(1+lag(temp6,'ret',8))*(1+lag(temp6,'ret',9))*(1+lag(temp6,'ret',10))*(1+lag(temp6,'ret',11))*(1+lag(temp6,'ret',12)))-1
temp6['mom36m'] = ((1+lag(temp6,'ret',2)) * (1+lag(temp6, 'ret', 3)) * (1+lag(temp6,'ret',4)) * (1+lag(temp6, 'ret', 5)) * (1+lag(temp6,'ret',6)) * \
    (1+lag(temp6,'ret',7))*(1+lag(temp6,'ret',8))*(1+lag(temp6,'ret',9))*(1+lag(temp6,'ret',10))*(1+lag(temp6,'ret',11))*(1+lag(temp6,'ret',12)) * \
    (1+lag(temp6,'ret',13))*(1+lag(temp6,'ret',14))*(1+lag(temp6,'ret',15))*(1+lag(temp6,'ret',16))*(1+lag(temp6,'ret',17))*(1+lag(temp6,'ret',18)) * \
    (1+lag(temp6,'ret',19))*(1+lag(temp6,'ret',20))*(1+lag(temp6,'ret',21))*(1+lag(temp6,'ret',22))*(1+lag(temp6,'ret',23))*(1+lag(temp6,'ret',24)) * \
    (1+lag(temp6,'ret',25))*(1+lag(temp6,'ret',26))*(1+lag(temp6,'ret',27))*(1+lag(temp6,'ret',28))*(1+lag(temp6,'ret',29))*(1+lag(temp6,'ret',30)) * \
    (1+lag(temp6,'ret',31))*(1+lag(temp6,'ret',32))*(1+lag(temp6,'ret',33))*(1+lag(temp6,'ret',34))*(1+lag(temp6,'ret',35))*(1+lag(temp6,'ret',36)))-1
temp6['mom1m'] = lag(temp6, 'ret', 1)
temp6['dolvol'] = np.log(lag(temp6, 'ret', 2) * lag(temp6,'prc',2))
temp6['chmom'] = ((1+lag(temp6,'ret',1)) * (1+lag(temp6,'ret',2)) * (1+lag(temp6,'ret',3)) * (1+lag(temp6,'ret',4)) * (1+lag(temp6,'ret',5)) * (1+lag(temp6,'ret',6)))-1\
    - ((1+lag(temp6,'ret',7))*(1+lag(temp6,'ret',8))*(1+lag(temp6,'ret',9))*(1+lag(temp6,'ret',10))*(1+lag(temp6,'ret',11))*(1+lag(temp6,'ret',12))-1)
temp6['turn'] = (lag(temp6,'vol',1) + lag(temp6,'vol',2) + lag(temp6,'vol',3))/(3*temp6['shrout'])
temp6['retcons_pos'] = 0
temp6.loc[(lag(temp6,'ret',1)>0) & (lag(temp6,'ret',2)>0) & (lag(temp6,'ret',3)>0) & (lag(temp6,'ret',4)>0) & (lag(temp6,'ret',5)>0) & (lag(temp6,'ret',6)>0), 'retcons_pos'] = 1
temp6['retcons_pos'] = 0
temp6.loc[(lag(temp6,'ret',1)<0) & (lag(temp6,'ret',2)<0) & (lag(temp6,'ret',3)<0) & (lag(temp6,'ret',4)<0) & (lag(temp6,'ret',5)<0) & (lag(temp6,'ret',6)<0), 'retcons_pos'] = 1

"""
	/*HXZ*/
	moms12m=(lag(ret)+lag2(ret)+lag3(ret)+lag4(ret)+lag5(ret)+lag6(ret)+lag7(ret)+lag8(ret)+lag9(ret)+lag10(ret)+lag11(ret))/11.0;
	moms60m=(lag13(ret)+lag14(ret)+lag15(ret)+lag16(ret)+lag17(ret)+lag18(ret)+lag19(ret)+lag20(ret)+lag21(ret)+lag22(ret)+lag23(ret)+
lag25(ret)+lag26(ret)+lag27(ret)+lag28(ret)+lag29(ret)+lag30(ret)+lag31(ret)+lag32(ret)+lag33(ret)+lag34(ret)+lag35(ret)+lag37(ret)+lag38(ret)+
lag39(ret)+lag40(ret)+lag41(ret)+lag42(ret)+lag43(ret)+lag44(ret)+lag45(ret)+lag46(ret)+lag47(ret)+lag49(ret)+lag50(ret)+lag51(ret)+lag52(ret)+
lag53(ret)+lag54(ret)+lag55(ret)+lag56(ret)+lag57(ret)+lag58(ret)+lag59(ret))/47.0;
	cei=log(mve6b/lag5(mve6b))-(log(lag(ret)+1)+log(lag1(ret)+1)+log(lag2(ret)+1)+log(lag3(ret)+1)+log(lag4(ret)+1)+log(lag5(ret)+1)+log(lag6(ret)+1) +log(lag7(ret)+1)+log(lag8(ret)+1)+log(lag9(ret)+1)+log(lag10(ret)+1)+log(lag11(ret)+1)+log(lag12(ret)+1)+log(lag13(ret)+1)+log(lag14(ret)+1)+log(lag15(ret)+1) +log(lag16(ret)+1)+log(lag17(ret)+1)+log(lag18(ret)+1)+log(lag19(ret)+1)+log(lag20(ret)+1)+log(lag21(ret)+1)+log(lag22(ret)+1)+log(lag23(ret)+1)+log(lag24(ret)+1) +log(lag25(ret)+1)+log(lag26(ret)+1)+log(lag27(ret)+1)+log(lag28(ret)+1)+log(lag29(ret)+1)+log(lag30(ret)+1)+log(lag31(ret)+1)+log(lag32(ret)+1)+log(lag33(ret)+1) +log(lag34(ret)+1)+log(lag35(ret)+1)+log(lag36(ret)+1)+log(lag37(ret)+1)+log(lag38(ret)+1)+log(lag39(ret)+1)+log(lag40(ret)+1)+log(lag41(ret)+1)+log(lag42(ret)+1) +log(lag43(ret)+1)+log(lag44(ret)+1)+log(lag45(ret)+1)+log(lag46(ret)+1)+log(lag47(ret)+1)+log(lag48(ret)+1)+log(lag49(ret)+1)+log(lag50(ret)+1)+log(lag51(ret)+1) +log(lag52(ret)+1)+log(lag53(ret)+1)+log(lag54(ret)+1)+log(lag55(ret)+1)+log(lag56(ret)+1)+log(lag57(ret)+1)+log(lag58(ret)+1)+log(lag59(ret)+1));

	/*HXZ quarterly*/
	x_1=ret-retx;
	dpqq=(x_1*shrout+lag(x_1)*lag(shrout)+lag2(x_1)*lag2(shrout)) / mvel1;
"""

temp6['moms12m'] = (lag(temp6,'ret',1) + lag(temp6,'ret',2) + lag(temp6,'ret',3) + lag(temp6,'ret',4) + lag(temp6,'ret',5) + lag(temp6,'ret',6) + lag(temp6,'ret',7) + lag(temp6,'ret',8) + lag(temp6,'ret',9) + lag(temp6,'ret',10) + lag(temp6,'ret',11))/11.0
temp6['x_1'] = temp6['ret'] - temp6['retx']
temp6['dpqq'] = (temp6['x_1']*temp6['shrout']+lag(temp6,'x_1')*lag(temp6,'shrout')+lag(temp6,'x_1',2)*lag(temp6,'shrout'))/temp6['mvel1']
"""
	if count=1 then
		mom1m=.;

	if count<3 then
		do;
			dpqq=.;
		end;

	if count<5 /*count<13*/
	then
		do;
			mom12m=.;
			chmom=.;
		end;

	if count<5 /*count<7*/
	then
		mom6m=.;

	if count<5 /*count<37*/
	then
		mom36m=.;

	if count<3 then
		dolvol=.;

	if count<4 then
		turn=.;

	if count<4 then
		chnanalyst=.;

	if count<5 /*count<7*/
	then
		retcons_pos=.;

	if count<5 /*count<7*/
	then
		retcons_neg=.;

	if count<5 /*count<=12*/
	then
		IPO=1;
	else
		IPO=0;

	if count<5 /*count<60*/
	then
		cei=.;
run;

proc sql;
	create table temp5 as select *, mean(mom12m) as indmom from temp6 group by 
		sic2, date;
quit;
"""




"""

/* check 8 */
*=====================================================================================================================


			finally, a few more directly from daily CRSP to create monthly variables


======================================================================================================================;

proc sql;
	create table dcrsp as select permno, year(date) as yr, month(date) as month, 
		max(ret) as maxret, std(ret) as retvol, /*skewness(of ret) as ts,*/
		mean((askhi-bidlo)/((askhi+bidlo)/2)) as baspread, std(log(abs(prc*vol))) as 
		std_dolvol, std(vol/shrout) as std_turn, mean(abs(ret)/(abs(prc)*vol)) as 
		ill, sum(vol=0) as countzero, n(permno) as ndays, sum(vol/shrout) as turn 
		from crsp.dsf group by permno, year(date), month(date);
quit;

proc sort data=dcrsp nodupkey;
	by permno yr month;
run;

data dcrsp;
	set dcrsp;
	zerotrade=(countzero+((1/turn)/480000))*21/ndays;
run;

*match to prior month to use lagged variables to predict returns;

proc sql;
	create table temp6 as select a.*, b.month, b.yr, b.maxret, b.retvol, /*HXZ*/
	/*b.ts,*/
	baspread, std_dolvol, std_turn, ill, zerotrade from temp5 a left join dcrsp b 
		on a.permno=b.permno and year(intnx('MONTH', date, -1))=b.yr and 
		month(intnx('MONTH', date, -1))=b.month;
quit;

proc sort data=temp6 nodupkey;
	by permno date;
run;

*---create beta from weekly returns---;

proc sql;
	create table dcrsp as select permno, intnx('WEEK', date, 0, 'end') as wkdt, 
		exp(sum(log(1+(ret))))-1 as wkret from crsp.dsf group by permno, calculated 
		wkdt;
quit;

proc sort data=dcrsp nodupkey;
	by permno wkdt;
run;

proc sql;
	create table dcrsp as select *, mean(wkret) as ewret from dcrsp group by wkdt;
quit;

data dcrsp;
	set dcrsp;
	where not missing(wkret) and not missing(ewret);
run;

proc sort data=temp6 out=lst(keep=permno date) nodupkey;
	where not missing(permno) and not missing(date);
	by permno date;
run;

proc sql;
	create table betaest as select a.*, b.wkret, b.ewret as ewmkt, b.wkdt from lst 
		a left join dcrsp b on a.permno=b.permno and intnx('MONTH', date, 
		-36)<=wkdt<=intnx('MONTH', date, -1);
quit;

*3 years of weekly returns;

proc sql;
	create table betaest as select * from betaest group by permno, date having 
		n(wkret)>=52;
quit;

*require at least 1 year of weekly returns;

proc sort data=betaest;
	by permno date wkdt;
run;

data betaest;
	set betaest;
	by permno date wkdt;
	retain count;
	ewmkt_l1=lag(ewmkt);
	ewmkt_l2=lag2(ewmkt);
	ewmkt_l3=lag3(ewmkt);
	ewmkt_l4=lag4(ewmkt);

	if first.date then
		count=1;
	else
		count+1;

	if count<5 then
		do;
			ewmkt_l1=.;
			ewmkt_l2=.;
			ewmkt_l3=.;
			ewmkt_l4=.;
		end;
run;

proc reg data=betaest outest=est noprint;
	by permno date;
	model wkret=ewmkt/adjrsq;
	output out=idiovolest residual=idioret;
	run;
	*two different approaches, one typical, the other including lagged market values to use as price delay measure;

proc reg data=betaest outest=est2 noprint;
	by permno date;
	model wkret=ewmkt ewmkt_l1 ewmkt_l2 ewmkt_l3 ewmkt_l4/adjrsq;
	output out=idiovolest residual=idioret;
	run;

proc sql;
	create table idiovolest as select permno, date, std(idioret) as idiovol from 
		idiovolest group by permno, date;
quit;

proc sort data=idiovolest nodupkey;
	where not missing(idiovol);
	by permno date;
run;

data est;
	set est;
	where not missing(permno) and not missing(date);
	beta=ewmkt;
run;

proc sql;
	create table temp7 as select a.*, b.beta, b.beta*b.beta as betasq, _adjrsq_ as 
		rsq1 from temp6 a left join est b on a.permno=b.permno and a.date=b.date;
quit;

proc sql;
	create table temp7 as select a.*, 1-(rsq1 / _adjrsq_) as pricedelay from temp7 
		a left join est2 b on a.permno=b.permno and a.date=b.date;
quit;

proc sql;
	create table temp7 as select a.*, b.idiovol from temp7 a left join idiovolest 
		b on a.permno=b.permno and a.date=b.date;
quit;

proc sort data=temp7 nodupkey;
	by permno date;
run;

*=============================================================================================


	So there we go, save this monster,
	after this we do a little clean up, but it is all in here now

=============================================================================================;

data temp;
	set temp7;

	if missing(eamonth) then
		eamonth=0;

	if missing(IPO) then
		IPO=0;
run;

/*===============================================================================================


This is to clean up the data a little,
primarily we want to limit the influence of extreme outliers which are in a lot of these RPS variables

=================================================================================================*/
/* delete the part for removing outliers */
proc sql;
	create table final3 as select permno, date, mvel1, ret, prc, shrout, count, 
		beta, betasq, chmom, chnanalyst, dolvol, idiovol, indmom, ipo, mom1m, mom6m, 
		mom12m, mom36m, mve0, pricedelay, turn, time_1, shrcd, exchcd, siccd, time_2, 
		sic2, gvkey, fyear, absacc, acc, age, agr, bm, BE, mve6b, bm_ia, cashdebt, 
		cashpr, cfp, cfp_ia, chatoia, chcsho, chempia, chinv, chpmia, convind, 
		currat, depr, divi, divo, dy, egr, ep, gma, grcapx, grltnoa, herf, hire, 
		invest, lev, lgr, mve_ia, operprof, mve_f, orgcap, pchcapx_ia, pchcurrat, 
		pchdepr, pchgm_pchsale, pchquick, pchsale_pchinvt, pchsale_pchrect, 
		pchsale_pchxsga, pchsaleinv, pctacc, ps, quick, rd, rd_mve, rd_sale, 
		realestate, roic, salecash, saleinv, salerec, secured, securedind, sgr, sin, 
		sp, tang, tb, datadate, datadate_a, datadate_q, aeavol, cash, chtx, cinvest, 
		ear, nincr, roaq, roavol, roeq, rsup, stdacc, stdcf, sue, rdq, ms, statpers, 
		chfeps, disp, fgr5yr, nanalyst, sfe, yr, month, baspread, ill, maxret, 
		retvol, std_dolvol, std_turn, zerotrade, /*HXZ*/
		orgcap_ia, adm, gad, rdm, rds, ol, rca, bca, etr, lfe, kz, ww, cdd, roavol_a, 
		ala, alm, moms12m, moms60m, ob_a, cinvest_a, dpia, noa, dnoa, pchcapx, 
		pchcapx3, cei, cdi, ivg, dcoa, dcol, dwc, dnca, dncl, dnco, dfin, ta, dsti, 
		dfnl, poa, nef, ndf, nxf, atm, cp, op, nop, BE/mve0 as bmj, em, ndp, ebp, 
		rna, pm, ato, cto, gpa, rmw, ole, opa, ola, cop, cla, os, zs, bi, oca, 
		oca_ia, ha, he, pchsale_pchinvt_hxz, pchsale_pchrect_hxz, pchgm_pchsale_hxz, 
		pchsale_pchxsga_hxz, realestate_hxz, secured_hxz, agr_hxz, grltnoa_hxz, 
		chcsho_hxz, pchcapx_ia_hxz, acc_hxz, egr_hxz, pctacc_hxz, lev_hxz, ep_hxz, 
		cfp_hxz, tb_hxz, /*HXZ quaterly*/
		rdmq, rdsq, olq, tanq, kzq, alaq, almq, laq, bmq, dmq, amq, epq, cpq, emq, 
		spq, ndpq, ebpq, wwq, rs, droeq, droaq, noaq, rnaq, pmq, atoq, ctoq, glaq, 
		oleq, olaq, claq, blq, sgq, dpqq, /*ts,*/
		/*HXZ monthly*/
		pps from temp where datadate>='01Jan2015'd;
quit;

data stemp(keep=gvkey sic sic2);
	set comp.names;
	sic2=substr(sic, 1, 2);
run;

proc sql;
	create table stemp1 as select a.*, b.lpermno as permno, b.linkprim, b.linkdt, 
		b.linkenddt from stemp as a, cc.CCMXPF_LNKHIST as b where a.gvkey=b.gvkey and 
		substr(b.linktype, 1, 1)='L' and
(2018>=year(LINKDT) or LINKDT=.B) and (1940<=year(LINKENDDT) or LINKENDDT=.E);
quit;

proc sort data=stemp1 out=stemp2;
	by permno descending linkprim gvkey;
run;

proc sql;
	create table ndata1 as select a.*, b.sic2 as sic2_new, b.linkprim from final3 
		a left join stemp2 b on a.permno=b.permno and (a.date>=b.linkdt or 
		b.linkdt=.B) and (b.linkenddt>=a.date or b.linkenddt=.E);
quit;

proc sort data=ndata1 out=ndata2;
	by date permno descending linkprim;
run;

data ndata3;
	set ndata2;
	by date permno descending linkprim;

	if first.permno;
run;

data ndata4;
	set ndata3;

	if missing(sic2_new) and not missing(sic2) then
		sic2_new=sic2;
run;

data ndata5(drop=linkprim sic2);
	set ndata4;
run;

data final3(rename=(sic2_new=sic2));
	set ndata5;
run;

*==============================================================================

I finally download and save the data here,

	if you are using this program, you need to save to a different location


	==============================================================================;

/*    Save data   */
proc export data=final3 outfile=outputpath dbms=csv replace;
run;"""


