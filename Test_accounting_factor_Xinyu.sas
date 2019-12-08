********************************

Major adjustments to get this code run:

1. getting ride of signin and signout part as we use sas studio directly
2. adding libname library reference for compr, which contains rating information
3. adding ibes 
4. updating comp.company to comp.names to acess naics and sic
5. adding work space cleaning process in the beginning 
6. remember to set you own path before running this file
7. delete one redundent 'dp' input from sql
********************************;

proc datasets library=work kill nolist;
quit;
run;

libname comp '/wrds/comp/sasdata/d_na';
libname cc '/wrds/crsp/sasdata/a_ccm';
libname compr '/wrds/comp/sasdata/d_na/rating';
libname ibes '/wrds/ibes/sasdata';
%Let outputpath="CHANGE THIS TO YOU OWN PATH PLEASE";
*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


Create signals (RPS) that are aligned in calendar month from firm characteristics pulled from Compustat, CRSP, and I/B/E/S
used to predict monthly cross-sections of stock returns

Creating this data includes every RPS in the paper plus a few more that we didn't use, mainly because there were a lot of missing observations for those RPS
	the ones not in the paper are from the RPS in the earlier paper Green, Hand, and Zhang (2013)


Some sort of liability disclaimer here...
If you use the program, make sure you understand it,
If you spot some errors, let is know,
and if you use it, please cite us :)

This is the data creation program for Green, Hand, and Zhang (2016)

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


*==============================================================================================================


				start with COMPUSTAT ANNUAL information


==============================================================================================================;

proc sql;
	create table data_0 as select 		/*header info*/
	substr(compress(f.cusip), 1, 6) as cnum, c.gvkey, datadate, datadate as 
		datadate_a, fyear, c.cik, cat(sic, 1, 2) as sic2, sic, naics, 
		/*firm variables*/
		/*income statement*/
		sale, revt, cogs, xsga, xrd, xad, ib, ebitda, ebit, nopi, spi, pi, txp, 
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
		datadate>='01Jan2016'd and f.indfmt='INDL' and f.datafmt='STD' and 
		f.popsrc='D' and f.consol='C';
quit;

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

/*sort and clean up*/
proc sort data=data_1 nodupkey;
	by gvkey datadate;
run;

%let dsevars=ticker ncusip shrcd exchcd siccd dlprc ;
%let dsfvars =prc vol ret retx shrout;

%macro crspmerge (s=, start=, end=, sfvars=, sevars=, filters=, final_ds=);
	/* Check Series: Daily or Monthly and define datasets - Default is Monthly  */
%if &s=M %then
		%let s=m;
	%else %if &s ne m %then
		%let s=d;

	%if (%sysfunc(libref(crsp))) %then
		%do;
			libname crsp"/wrds/crsp/sasdata/a_stock";
		%end;
	%let sf       = crsp.&s.sf ;
	%let se       = crsp.&s.seall ;
	%let senames  = crsp.&s.senames ;
	%put ;
	%put ;
	%put ;
	%put ;
	%put ;
	%put #### ## # Merging CRSP Stock File (&s.sf) and Event File (&s.se) # ## #### ;
	options nonotes;
	%let sdate = %sysfunc(putn("&start"d, 5.)) ;
	%let edate = %sysfunc(putn("&end"d, 5.)) ;

	%if %length(&sevars) > 2 %then
		%let sevars = %sysfunc(compbl(&sevars));
	%else
		%let sevars = ticker ncusip exchcd shrcd siccd dlprc;
	%let sevars  = %sysfunc(lowcase(&sevars));
	%let nsevars = %eval(%sysfunc(length(&sevars))-%sysfunc(length(%sysfunc(compress(&sevars))))+1);
	%* create lag event variable names to be used in the RETAIN statement;
	%let sevars_l = lag_%sysfunc(tranwrd(&sevars, %str( ), %str( lag_)));

	%if %length(&filters) > 2 %then
		%let filters = and &filters;
	%else
		%let filters = %str( );

	%if &final_ds = %str() %then
		%let final_ds = work.crsp_&s.;
	%put #### ### ## # ;

	/* Get stock data */
	proc sql;
		create table sfdata_1 as select *, abs(prc)*shrout as meq 
			from &sf (keep=permco permno date &sfvars) where &sdate<=date<=&edate and 
			missing(prc)=0 and permno in
(select distinct permno from
&senames(WHERE=(&edate>=NAMEDT and &sdate<=NAMEENDT) keep=permno namedt 
			nameendt) ) order by date, permco;
	quit;

	proc sort data=sfdata_1 out=sfdata_2;
		by date permco permno descending meq;
	run;

	proc sort data=sfdata_2 out=sfdata_3 nodupkey;
		by date permco permno;
	run;

	proc sort data=sfdata_3;
		by date permco meq;
	run;

	data sfdata_4 (drop=meq);
		set sfdata_3;
		by date permco meq;
		retain ME;

		if first.permco and last.permco then
			do;
				ME=meq;
				output;
				* most common case where a firm has a unique permno;
			end;
		else
			do;

				if first.permco then
					ME=meq;
				else
					ME=sum(meq, ME);

				if last.permco then
					output;
			end;
	run;

	proc sql;
		create table sfdata_5 as select a.*, b.ME from sfdata_1 as a left join 
			sfdata_4 as b on a. permco=b. permco and intnx('day', a.date, 0, 
			'E')=intnx('day', b.date, 0, 'E') order by permno, date;
	quit;

	proc sort data=sfdata_5 out=sfdata_6 nodupkey;
		by permno date;
	run;

	data sfdata_7;
		set sfdata_6;
		retain prca lprc LME;
		prca=abs(prc);
		lprc=lag(prca);
		LME=lag(me);
	run;

	proc sql;
		create table sfdata_8 as select a.*, b.prca, b.lprc, b.LME from sfdata_5 as a 
			left join sfdata_7 as b on a.permno=b.permno and intnx('day', a.date, 0, 
			'E')=intnx('day', b.date, 0, 'E') order by permno, date;
	quit;

	proc sql;
		create table sfdata as select * from sfdata_8 where missing(ret)=0 or ret=.B 
			order by permno, date;
	quit;

	%put #### ### ## # ;

	/* Get event data */
	proc sql;
		create table sedata as select a.* from &se (keep=permco permno 
			date &sevars exchcd) as a, (select distinct permno, min(namedt) as minnamedt 
			from
&senames(WHERE=(&edate>=NAMEDT and &sdate<=NAMEENDT) keep=permno namedt 
			nameendt) group by permno) as b where a.date >=b.minnamedt and 
			a.date <=&edate and a.permno=b.permno order by a.permno, a.date;
	quit;

	%put #### ### ## # ;
	** get the distribution info;

	data distr_0 (keep=date permno distcd divamt facpra facshr);
		set &se;
		where &sdate <=date <=&edate;
		retain facpra;
		facpra=facpr*1;
	run;

	proc sort data=distr_0;
		by date permno;
	run;

	data distr;
		set distr_0;

		if missing(DISTCD)=1 then
			delete;
	run;

	proc sql;
		update distr set date='17sep2001'd where '11sep2001'd <=date <='17sep2001'd;
		update distr set date='31oct2012'd where '29oct2012'd <=date <='31oct2012'd;
		update distr set date=intnx('day', date, 2) where weekday(date)=7;
		update distr set date=intnx('day', date, 1) where weekday(date)=1;
	quit;

	proc sort data=distr;
		by date permno;
	run;

	%put #### ### ## # ;

	/* Merge stock and event data */
%let eventvars = ticker comnam ncusip shrout siccd exchcd shrcls shrcd shrflg trtscd nmsind mmcnt nsdinx ;
	* variables whose values need to be retain to fill the blanks;

	data &final_ds (keep=permco permno date &sfvars &sevars exchcd ME prca lprc 
			LME);
		merge sedata (in=eventdata) sfdata (in=stockdata);
		by permno date;
		retain &sevars_l;

		%do i=1 %to &nsevars;
			%let var   = %scan(&sevars, &i, %str( ));
			%let var_l = %scan(&sevars_l, &i, %str( ));

			%if %sysfunc(index(&eventvars, &var))>0 %then
				%do;

					if eventdata or first.permno then
						&var_l=&var.;
					else if not eventdata then
						&var=&var_l.;
				%end;
		%end;

		if eventdata and not stockdata then
			delete;
		drop  &sevars_l;
	run;

	%put #### ### ## # ;

	/* ------------------------------------------------------------------------------ */
	/* The following sort is included to handle duplicate observations when a company */
	/* has more than one distribution on a given date. For example, a stock and cash  */
	/* distribution on the same date will generate two records, identical except for  */
	/* different DISTCD and DISTAMT (and possibly other) values. The NODUPLICATES     */
	/* option only deletes a record if all values for all variables are the same as   */
	/* those in another record. So, in the above example, if DISTCD is included in    */
	/* &sevars a record will not be deleted, but a redundant record will be deleted   */
	/* if DISTCD and DISTAMT are not included in &sevars.                             */
	/* ------------------------------------------------------------------------------ */
	proc sort data=&final_ds noduplicates;
		/* the "exchcd" condition below removes rows with empty stock price data created  */
		/* because CRSP event file track some event information even before the stock     */
		/* is trading in major stock exchange                                             */
		where 1 &filters;
		by date permno;
	run;

	%put #### ### ## # ;
	options notes;
	%put #### ## # Done: Dataset &final_ds Created! # ## ####;
	%put ;
%mend crspmerge;

%crspmerge(s=m, start=01July2017, end=30jun2018, sfvars=&dsfvars, 
	sevars=&dsevars, filters=exchcd in (1, 2, 3) and shrcd in (10, 11));

/*new*/
/* Add CRSP delisting returns */
proc sql;
	create table crsp_m2 as select a.*, b.dlret, sum(1, ret)*sum(1, dlret)-1 as 
		retadj, abs(a.prc)*a.shrout as meq from CRSP_M as a left join 
		crsp.dsedelist (where=(missing(dlret)=0) ) as b on a.permno=b.permno and 
		intnx('day', a.date, 0, 'E')=intnx('day', b.DLSTDT, 0, 'E') order by permno, 
		date;
quit;

/* The next step does 2 things:                                        */
/* - Create weights for later calculation of VW returns.               */
/*   Each firm's daily return RET t will be weighted by             */
/*   ME t-1 = ME t-2 * (1 + RETX t-1)                                  */
/*     where RETX is the without-dividend return.                      */
/* - Create a File with December t-1 Market Equity (ME)                */
data crsp_m3 (drop=ncusip LME Lpermno me_base retx year) decme (keep=permno 
		year date ME rename=(me=mve6b) );
	set crsp_m2;
	by permno date;
	retain weight_port cumretx me_base year;
	Lpermno=lag(permno);
	LME=lag(me);
	year=year(date);

	if month(date)=7 then
		do;
			weight_port=LME;
			me_base=LME;

			/* lag ME also at the end of June */
			cumretx=sum(1, retx);
		end;
	else
		do;

			if LME>0 then
				weight_port=cumretx*me_base;
			else
				weight_port=.;
			cumretx=cumretx*sum(1, retx);
		end;
	output crsp_m3;

	if month(date)=12 and ME>0 then
		output decme;
run;

/* Cleaning DECME Data    */
proc sort data=decme out=decme1;
	by permno year descending date;
run;

/* Keep the last day ME of December*/
proc sort data=decme1 out=decme2 nodupkey;
	by permno year;
run;

/* End cleaning DECME       */
/* Create a file with data for each June with ME from previous December */
proc sql;
	create table crsp_june as select a.*, b.mve6b /*, year(a.date) as year */
	from crsp_m3(where=(month(date)=6)) as a, decme2 as b where a.permno=b.permno 
		and intck('month', b.date, a.date)=6;
quit;

run;

/* Add Permno to Compustat sample */
proc sql;
	create table ccm1_beme as select a.*, b.lpermno as permno, b.linkprim from 
		data_1 as a, cc.CCMXPF_LNKHIST as b where a.gvkey=b.gvkey and 
		substr(b.linktype, 1, 1)='L' and linkprim in ('P', 'C') and (intnx('month', 
		intnx('year', a.datadate, 0, 'E'), 6, 'E') >=b.linkdt) 
		and (b.linkenddt >=intnx('month', intnx('year', a.datadate, 0, 'E'), 6, 'E') 
		or missing(b.linkenddt)) order by a.datadate, permno, b.linkprim desc;
quit;

/*  Cleaning Compustat Data for no relevant duplicates                      */
/*  Eliminating overlapping matching : few cases where different gvkeys     */
/*  for same permno-date --- some of them are not 'primary' matches in CCM  */
/*  Use linkprim='P' for selecting just one gvkey-permno-date combination   */
proc sort data=ccm1_beme out=ccm2_beme;
	by datadate permno descending linkprim gvkey;
run;

data ccm1a_beme;
	set ccm2_beme;
	by datadate permno descending linkprim gvkey;

	if first.permno;
run;

/* Sanity Check -- No Duplicates */
proc sort data=ccm1a_beme nodupkey;
	by datadate permno;
run;

/* 2. However, there other type of duplicates within the year                */
/* Some companiess change fiscal year end in the middle of the calendar year */
/* In these cases, there are more than one annual record for accounting data */
/* We will be selecting the last annual record in a given calendar year      */
proc sort data=ccm1a_beme;
	by permno year datadate;
run;

data ccm2a_beme;
	set ccm1a_beme;
	by permno year datadate;

	if last.year=1;
run;

/* Sanity Check -- No Duplicates */
proc sort data=ccm2a_beme nodupkey;
	by permno datadate;
run;

/* Finalize Compustat Sample                              */
/* Merge CRSP with Compustat data, at June of every year  */
/* Match fiscal year ending calendar year t-1 with June t */
proc sql;
	create table data as select a.*, b.mve6b, b.date, (1000*a.BE)/b.mve6b as bm, 
		intck('month', a.datadate, b.date) as dist from ccm2a_beme a left join 
		crsp_june b on a.permno=b.permno and intnx('month', b.date, 0, 
		'E')=intnx('month', intnx('year', a.datadate, 0, 'E'), 6, 'E');
quit;

proc sort data=data;
	by permno datadate descending date;
run;

proc sort data=data nodupkey;
	by permno datadate;
run;

/* check 1 */
/*sort and clean up*/
proc sort data=data nodupkey;
	by gvkey datadate;
run;

/*prep for clean-up and using time series of variables*/
data data;
	set data;
	retain count;
	by gvkey datadate;

	if first.gvkey then
		count=1;
	else
		count+1;
run;

data data;
	set data;
	*do some clean up, several of these variables have lots of missing values;

	if not missing(drc) and not missing(drlt) then
		dr=drc+drlt;

	if not missing(drc) and missing(drlt) then
		dr=drc;

	if not missing(drlt) and missing(drc) then
		dr=drlt;

	if missing(dcvt) and not missing(dcpstk) and not missing(pstk) and dcpstk>pstk 
		then
			dc=dcpstk-pstk;

	if missing(dcvt) and not missing(dcpstk) and missing(pstk) then
		dc=dcpstk;

	if missing(dc) then
		dc=dcvt;

	if missing(xint) then
		xint=0;

	if missing(xsga) then
		xsga=0;
run;


/*--------------------------------------------------------

more clean-up and create first pass of variables

----------------------------------------------------------*/
data data2;
	set data;

	/*create simple-just annual Compustat variables*/
	ep=ib/mve_f;
	cashpr=((mve_f+dltt-at)/che);
	dy=dvt/mve_f;
	lev=lt/mve_f;
	sp=sale/mve_f;
	roic=(ebit-nopi)/(ceq+lt-che);
	rd_sale=xrd/sale;
	rd_mve=xrd/mve_f;
	*agr=1-(at/lag(at));
	lat=lag(at);
	agr=(lat-at)/lat;

	if missing(at) or missing(lat) then
		agr=.;
	gma=(revt-cogs)/lag(at);
	chcsho=(csho/lag(csho))-1;
	lgr=(lt/lag(lt))-1;
	acc=(ib-oancf) / ((at+lag(at))/2);

	if missing(oancf) then
		acc=((act-lag(act) - (che-lag(che))) 
			- ((lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp) )/ ((at+lag(at))/2);
	pctacc=(ib-oancf)/abs(ib);

	if ib=0 then
		pctacc=(ib-oancf)/.01;

	if missing(oancf) then
		pctacc=((act-lag(act) - (che-lag(che))) 
			- ((lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp) )/abs(ib);

	if missing(oancf) and ib=0 then
		pctacc=((act-lag(act) - (che-lag(che))) 
			- ((lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp) )/.01;
	cfp=(ib-((act-lag(act) - (che-lag(che))) 
		- ((lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp) ))/mve_f;







	if not missing(oancf) then
		cfp=oancf/mve_f;
	absacc=abs(acc);
	age=count;
	chinv=(invt-lag(invt))/((at+lag(at))/2);

	if spi ne 0 and not missing(spi) then
		spii=1;
	else
		spii=0;
	spi=spi/ ((at+lag(at))/2);
	cf=oancf/((at+lag(at))/2);

	if missing(oancf) then
		cf=(ib-((act-lag(act) - (che-lag(che))) 
			- ((lct-lag(lct))-(dlc-lag(dlc))-(txp-lag(txp))-dp) ))/((at+lag(at))/2);
	hire=(emp-lag(emp))/lag(emp);

	if missing(emp) or missing(lag(emp)) then
		hire=0;
		
		
	sgr=(sale/lag(sale))-1;
	chpm=(ib/sale)-(lag(ib)/lag(sale));
	chato=(sale/((at+lag(at))/2)) - (lag(sale)/((lag(at)+lag2(at))/2));
	pchsale_pchinvt=((sale-lag(sale))/lag(sale))-((invt-lag(invt))/lag(invt));
	pchsale_pchrect=((sale-lag(sale))/lag(sale))-((rect-lag(rect))/lag(rect));
	pchgm_pchsale=(((sale-cogs)-(lag(sale)-lag(cogs)))/(lag(sale)-lag(cogs)))-((sale-lag(sale))/lag(sale));
	pchsale_pchxsga=((sale-lag(sale))/lag(sale) )-((xsga-lag(xsga)) /lag(xsga) );
	depr=dp/ppent;
	pchdepr=((dp/ppent)-(lag(dp)/lag(ppent)))/(lag(dp)/lag(ppent));
	chadv=log(1+xad)-log((1+lag(xad)));
	*had error here before, might work better now...;
	invest=((ppegt-lag(ppegt)) +  (invt-lag(invt)) ) / lag(at);

	if missing(ppegt) then
		invest=((ppent-lag(ppent)) +  (invt-lag(invt)) ) / lag(at);
	egr=((ceq-lag(ceq))/lag(ceq) );

	if missing(capx) and count>=2 then
		capx=ppent-lag(ppent);
	pchcapx=(capx-lag(capx))/lag(capx);
	grcapx=(capx-lag2(capx))/lag2(capx);
	grGW=(gdwl-lag(gdwl))/lag(gdwl);

	if missing(gdwl) or gdwl=0 then
		grGW=0;

	if gdwl ne 0 and not missing(gdwl) and missing(grGW) then
		grGW=1;

	if (not missing(gdwlia) and gdwlia ne 0) or (not missing(gdwlip) and gdwlip ne 
		0) or (not missing(gwo) and gwo ne 0) then
			woGW=1;
	else
		woGW=0;
	tang=(che+rect*0.715+invt*0.547+ppent*0.535)/at;

	if (2100<=sic<=2199) or (2080<=sic<=2085) or (naics in ('7132', '71312', 
		'713210', '71329', '713290', '72112', '721120')) then
			sin=1;
	else
		sin=0;

	if missing(act) then
		act=che+rect+invt;

	if missing(lct) then
		lct=ap;
	currat=act/lct;
	pchcurrat=((act/lct)-(lag(act)/lag(lct)))/(lag(act)/lag(lct));
	quick=(act-invt)/lct;
	pchquick=((act-invt)/lct - (lag(act)-lag(invt))/lag(lct) )/ 
		((lag(act)-lag(invt) )/lag(lct) );
	salecash=sale/che;
	salerec=sale/rect;
	saleinv=sale/invt;
	pchsaleinv=((sale/invt)-(lag(sale)/lag(invt)) ) / (lag(sale)/lag(invt));
	cashdebt=(ib+dp)/((lt+lag(lt))/2);
	realestate=(fatb+fatl)/ppegt;

	if missing(ppegt) then
		realestate=(fatb+fatl)/ppent;

	if (not missing(dvt) and dvt>0) and (lag(dvt)=0 or missing(lag(dvt))) then
		divi=1;
	else
		divi=0;

	if (missing(dvt) or dvt=0) and (lag(dvt)>0 and not missing(lag(dvt))) then
		divo=1;
	else
		divo=0;
	obklg=ob/((at+lag(at))/2);
	chobklg=(ob-lag(ob))/((at+lag(at))/2);

	if not missing(dm) and dm ne 0 then
		securedind=1;
	else
		securedind=0;
	secured=dm/dltt;

	if not missing(dc) and dc ne 0 or (not missing(cshrc) and CSHRC ne 0) then
		convind=1;
	else
		convind=0;
	conv=dc/dltt;
	grltnoa=((rect+invt+ppent+aco+intan+ao-ap-lco-lo)-(lag(rect)+lag(invt)+lag(ppent)+lag(aco)+lag(intan)+lag(ao)-lag(ap)-lag(lco)-lag(lo)) 
		-(rect-lag(rect)+invt-lag(invt)+aco-lag(aco)-(ap-lag(ap)+lco-lag(lco)) 
		-dp))/((at+lag(at))/2);
	chdrc=(dr-lag(dr))/((at+lag(at))/2);

	if ((xrd/at)-(lag(xrd/lag(at))))/(lag(xrd/lag(at))) >.05 then
		rd=1;
	else
		rd=0;
	rdbias=(xrd/lag(xrd))-1-ib/lag(ceq);
	roe=ib/lag(ceq);

	
	ps_beme=coalesce(pstkrv, pstkl, pstk, 0);

	if missing(txditc) then
		txditc=0;
	BE=ceq + txditc - ps_beme;

	if BE<0 then
		BE=.;
	operprof=(revt-cogs-xsga-xint)/BE;

	if missing(revt) then
		operprof=.;

	if missing(cogs)=1 and missing(xsga)=1 and missing(xint)=1 then
		operprof=.;

	if missing(BE) then
		operprof=.;
	ps=(ni>0)+(oancf>0)+(ni/at > lag(ni)/lag(at))+(oancf>ni)+(dltt/at < lag(dltt)/lag(at))+(act/lct > lag(act)/lag(lct)) +((sale-cogs)/sale > (lag(sale)-lag(cogs))/lag(sale))+ (sale/at > lag(sale)/lag(at))+ (scstkc=0);
	*-----Lev and Nissim (2004);

	if fyear<=1978 then
		tr=.48;

	if 1979<=fyear<=1986 then
		tr=.46;

	if fyear=1987 then
		tr=.4;

	if 1988<=fyear<=1992 then
		tr=.34;

	if 1993<=fyear then
		tr=.35;
	tb_1=((txfo+txfed)/tr)/ib;

	if missing(txfo) or missing(txfed) then
		tb_1=((txt-txdi)/tr)/ib;
	*they rank within industries;

	if (txfo+txfed>0 or txt>txdi) and ib<=0 then
		tb_1=1;






	*variables that will be used in subsequent steps to get to final RPS;
	*--prep for for Mohanram (2005) score;
	roa=ni/((at+lag(at))/2);
	cfroa=oancf/((at+lag(at))/2);

	if missing(oancf) then
		cfroa=(ib+dp)/((at+lag(at))/2);
	xrdint=xrd/((at+lag(at))/2);
	capxint=capx/((at+lag(at))/2);
	xadint=xad/((at+lag(at))/2);

	/*HXZ*/
	adm=xad/mve6b;
	gad=(xad-lag(xad))/lag(xad);
	rdm=xrd/mve6b;
	rds=xrd/sale;
	ol=(cogs+xsga)/at;
	rc_1=xrd+0.8*lag(xrd)+0.6*lag2(xrd)+0.4*lag3(xrd)+0.2*lag4(xrd);
	rca=rc_1/at;
	x_1=txt/(pi+am);
	eps_1=ajex/prcc_f;
	etr=(x_1-(lag(x_1)+lag2(x_1)+lag3(x_1))/3)*(eps_1-lag(eps_1));
	x_2=sale/emp;
	lfe=(x_2-lag(x_2))/lag(x_2);
	kz=-1.002*(ib+dp)/lag(ppent)+0.283*(at+mve6b-ceq-txdb)/at+3.139*(dlc+dltt)/(dlc+dltt+seq)-39.368*(dvc+dvp)/lag(ppent)-1.315*che/lag(ppent);
	cdd=dcvt/(dlc+dltt);
	roaq_a=ib/lag(at);
	roavol_1=std(roaq_a, lag(roaq_a), lag2(roaq_a), lag3(roaq_a), lag4(roaq_a), 
		lag5(roaq_a), lag6(roaq_a), lag7(roaq_a), lag8(roaq_a), lag9(roaq_a));
	cs_1=(ib-(act-lag(act)-(lct-lag(lct))-(che-lag(che))+dlc-lag(dlc)))/lag(at);
	roavol_2=std(cs_1, lag(cs_1), lag2(cs_1), lag3(cs_1), lag4(cs_1), lag5(cs_1), 
		lag6(cs_1), lag7(cs_1), lag8(cs_1), lag9(cs_1));
	roavol_a=roavol_1/roavol_2;

	if missing(gdwl) then
		gdwl=0;

	if missing(intan) then
		intan=0;
	ala=che+0.75*(act-che)-0.5*(at-act-gdwl-intan);
	alm=ala/(at+prcc_f*csho-ceq);
	ob_a=ob/(.5*at+.5*lag(at));
	x_3=capx/sale;
	cinvest_a=x_3/((lag(x_3)+lag2(x_3)+lag3(x_3))/3)-1;
	dpia=(ppegt-lag(ppegt)+invt-lag(invt))/lag(at);

	if missing(dlc) then
		dlc=0;

	if missing(dltt) then
		dltt=0;

	if missing(mib) then
		mib=0;

	if missing(pstk) then
		pstk=0;

	if missing(ceq) then
		ceq=0;
	noa=((at-che)-(at-dlc-dltt-mib-pstk-ceq))/lag(at);
	dnoa=noa-lag(noa);
	pchcapx3=capx/lag3(capx)-1;
	x_4=dlc+dltt;
	cdi=log(x_4/lag5(x_4));
	ivg=invt/lag(invt)-1;
	dcoa=(act-lag(act)-(che-lag(che)))/lag(at);
	dcol=(lct-lag(lct)-(dlc-lag(dlc)))/lag(at);
	dwc=(dcoa-dcol)/lag(at);
	dnca=(at-act-ivao-(lag(at)-lag(act)-lag(ivao)))/lag(at);
	dncl=(lt-lct-dltt-(lag(lt)-lag(lct)-lag(dltt)))/lag(at);
	dnco=(dnca-dncl)/lag(at);
	dfin=(ivst+ivao-dltt-dlc-pstk-(lag(ivst)+lag(ivao)-lag(dltt)-lag(dlc)-lag(pstk)))/lag(at);
	ta=(dwc+dnco+dfin)/lag(at);
	dsti=(ivst-lag(ivst))/lag(at);
	dfnl=(dltt+dlc+pstk-(lag(dltt)+lag(dlc)+lag(pstk)))/lag(at);
	egr_hxz=(ceq-lag(ceq))/lag(at);

	if missing(txp) then
		txp=0;
	poa=(act-lag(act)-(che-lag(che))-(lct-lag(lct)-(dlc-lag(dlc))-(txp-lag(txp)))-dp)/abs(ni);
	nef=(sstk-prstkc-dv)/((at+lag(at))/2);

	if missing(dlcch) then
		dlcch=0;
	ndf=(dltis-dltr+dlcch)/((at+lag(at))/2);
	nxf=nef+ndf;
	atm=at/mve6b;
	cp=(ib+dp)/mve6b;
	op=(dvc+prstkc-(pstkrv-lag(pstkrv)))/mve6b;
	nop=(dvc+prstkc-(pstkrv-lag(pstkrv))-sstk+pstkrv-lag(pstkrv))/mve6b;
	em=(mve6b+dlc+dltt+pstkrv-che)/oibdp;

	if missing(dvpa) then
		dvpa=0;

	if missing(tstkp) then
		tstkp=0;
	ndp=dltt+dlc+pstk+dvpa-tstkp-che;
	ebp=(ndp+ceq+tstkp-dvpa)/(ndp+mve6b);
	rna=oiadp/lag(noa);
	pm=rna/sale;
	ato=sale/lag(noa);
	cto=sale/lag(at);
	gpa=(revt-cogs)/at;
	rmw=(revt-cogs-xsga-xint)/(ceq+pstk);
	ole=(revt-cogs-xsga-xint)/(lag(ceq)+lag(pstk));
	opa=(revt-cogs-xsga+xrd)/at;
	ola=(revt-cogs-xsga+xrd)/lag(at);
	cop=(revt-cogs-xsga+xrd-(rect-lag(rect))-(invt-lag(invt))-(xpp-lag(xpp))+drc-lag(drc)+drlt-lag(drlt)+ap-lag(ap)+xacc-lag(xacc))/at;
	cla=(revt-cogs-xsga+xrd-(rect-lag(rect))-(invt-lag(invt))-(xpp-lag(xpp))+drc-lag(drc)+drlt-lag(drlt)+ap-lag(ap)+xacc-lag(xacc))/lag(at);

	if lt>at then
		i_1=1;
	else
		i_1=0;

	if ni<0 and lag(ni)<0 then
		i_2=1;
	else
		i_2=0;
	os=-1.32-0.407*log(at)+6.03*(dlc+dltt)/at-1.43*(act-lct)/at+0.076*(lct/act)-1.72*i_1-2.37*ni/at-1.83*(pi+dp)/lt+0.285*i_2-0.521*(ni+lag(ni))/(abs(ni)+abs(lag(ni)));
	zs=1.2*(act-lct)/at+1.4*re/at+3.3*oiadp/at+0.6*mve6b/lt+sale/at;

	if BE>0 then
		bi=at/BE;

	if missing(bi) then
		bi=.;
	pchsale_pchinvt_hxz=(sale-lag(sale))/(0.5*sale+0.5*lag(sale))-(invt-lag(invt))/(0.5*invt+0.5*lag(invt));
	pchsale_pchrect_hxz=(sale-lag(sale))/(0.5*sale+0.5*lag(sale))-(rect-lag(rect))/(0.5*rect+0.5*lag(rect));
	gm_1=sale-cogs;
	pchgm_pchsale_hxz=(gm_1-lag(gm_1))/(0.5*(gm_1+lag(gm_1)))-(sale-lag(sale))/(0.5*sale+0.5*lag(sale));
	pchsale_pchxsga_hxz=(sale-lag(sale))/(0.5*sale+0.5*lag(sale))-(xsga-lag(xsga))/(0.5*xsga+0.5*lag(xsga));
	realestate_hxz=(fatb+fatl)/ppegt;

	if missing(fatb) then
		realestate_hxz=(ppenb+ppenls)/ppent;
	secured_hxz=dm/(dltt+dlc);
	agr_hxz=at/lag(at)-1;
	x_1=ppent+intan+ao-lo+dp;
	grltnoa_hxz=(x_1-lag(x_1))/((at+lag(at))/2);
	chcsho_hxz=log(csho*ajex)-log(lag(csho)*lag(ajex));
	pchcapx_hxz=(capxv-0.5*lag(capxv)-0.5*lag2(capxv))/(0.5*lag(capxv)+0.5*lag2(capxv));

	if missing(txp) then
		txp=0;
	acc_hxz=(act-lag(act)-(che-lag(che))-(lct-lag(lct)-(dlc-lag(dlc))-(txp-lag(txp)))-dp)/lag(at);
	pctacc_hxz=(dwc+dnco+dfin)/abs(ni);
	lev_hxz=(dlc+dltt)/mve6b;
	ep_hxz=ib/mve6b;
	cfp_hxz=(fopt-(wcap-lag(wcap)))/mve6b;

	if not missing(oancf) then
		cfp_hxz=(fopt-oancf)/mve6b;
	tb_hxz=pi/ni;

	/*clean up for observations that do not have lagged observations to create variables*/
	array req{*} chadv agr invest gma chcsho lgr egr chpm chinv hire cf acc pctacc 
		absacc spi sgr pchsale_pchinvt pchsale_pchrect pchgm_pchsale pchsale_pchxsga 
		pchcapx ps roa cfroa xrdint capxint xadint divi divo

		/*HXZ*/
		gad lfe kz ob_a dpia noa ivg dcoa dcol dwc dnca dncl dnco dfin ta dsti dfnl 
		poa nef ndf nxf op nop cto ole ola cop cla os pchsale_pchinvt_hxz 
		pchsale_pchrect_hxz pchgm_pchsale_hxz pchsale_pchxsga_hxz agr_hxz grltnoa_hxz 
		chcsho_hxz acc_hxz egr_hxz pctacc_hxz cfp_hxz obklg chobklg grltnoa chdrc rd 
		pchdepr grGW pchcurrat pchquick pchsaleinv roe operprof;

	if count=1 then
		do;

			do b=1 to dim(req);
				req(b)=.;
			end;
		end;

	if count<3 then
		do;
			chato=.;
			grcapx=.;
			dnoa=.;
			rna=.;
			pm=.;
			ato=.;
			pchcapx_hxz=.;
		end;

	if count<4 then
		do;
			etr=.;
			cinvest_a=.;
			pchcapx3=.;
		end;

	if count<5 then
		do;
			rca=.;
		end;

	if count<5 /*count<6*/
	then
		do;
			cdi=.;
		end;

	if count<5 /*count<10*/
	then
		do;
			roavol_a=.;
		end;
run;
