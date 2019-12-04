********************************

Major adjustments to get this code run:

1. getting ride of signin and signout part as we use sas studio directly
2. adding libname library reference for compr, which contains rating information
3. adding ibes 
4. updating comp.company to comp.names to acess naics and sic
5. adding work space cleaning process in the beginning 
6. remember to set you own path before running this file

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
		datadate>='01Jan2018'd and f.indfmt='INDL' and f.datafmt='STD' and 
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

%crspmerge(s=m, start=01Jan2016, end=30jun2018, sfvars=&dsfvars, 
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