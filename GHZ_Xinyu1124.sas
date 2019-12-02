********************************

Major adjustments to get this code run:

1. getting ride of signin and signout part as we use sas studio directly
2. adding libname library reference for compr, which contains rating information
3. adding ibes 
4. updating comp.company to comp.names to acess naics and sic
5. adding work space cleaning process in the beginning 
6. correcting some missing semicolons and typos

Last note: do remember to set you own path before running this file


********************************;

proc datasets library=work kill nolist;
quit;
run;

libname comp '/wrds/comp/sasdata/d_na';
libname cc '/wrds/crsp/sasdata/a_ccm';
libname compr '/wrds/comp/sasdata/d_na/rating';
libname ibes '/wrds/ibes/sasdata';
outputpath="CHANGE THIS TO YOU OWN PATH PLEASE";
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
		datadate>='01Jan2015'd and f.indfmt='INDL' and f.datafmt='STD' and 
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

%crspmerge(s=m, start=01Jan2015, end=30jun2018, sfvars=&dsfvars, 
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

/*
look at how many missing;
data test;
set data;
where datadate>='01Jan2015'd;
array lst{*} xrd nopi xad dvt ob dm dc aco ap intan ao lco lo rect invt drc drlt dr spi gdwl emp dm dcvt fatb fatl che dp lct act xsga at;
array dms{*} dxrd dnopi dxad ddvt dob ddm ddc daco dap dintan dao dlco dlo drect dinvt ddrc ddrlt ddr dspi dgdwl demp ddm ddcvt dfatb dfatl dche ddp dlct dact dxsga dat;
do i=1 to dim(lst);
if lst(i)=. then dms(i)=1; else dms(i)=0;
end;
run;
proc means data=test mean sum ;
var dxrd dnopi dxad ddvt dob ddm ddc daco dap dintan dao dlco dlo drect dinvt ddrc ddrlt ddr dspi dgdwl demp ddm ddcvt dfatb dfatl dche ddp dlct dact dxsga dat;
run;
endrsubmit;	*/
*xsga also has a fair amount missing...;

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

/*other preparation steps for annual variables: industry adjustments*/
proc sql;
	create table data2 as select *, chpm-mean(chpm) as chpmia, chato-mean(chato) 
		as chatoia, sum(sale) as indsale, hire-mean(hire) as chempia, bm-mean(bm) as 
		bm_ia, pchcapx-mean(pchcapx) as pchcapx_ia, tb_1-mean(tb_1) as tb, 
		cfp-mean(cfp) as cfp_ia, mve_f-mean(mve_f) as mve_ia, /*HXZ*/
		sum(at) as indat, sum(BE) as indbe, pchcapx_hxz-mean(pchcapx_hxz) as 
		pchcapx_ia_hxz from data2 group by sic2, fyear;
quit;

proc sql;
	create table data2 as select *, sum((sale/indsale)*(sale/indsale) ) as herf, 
		/*HXZ*/
		sum((at/indat)*(at/indat)) as ha, sum((BE/indbe)*(BE/indbe)) as he from data2 
		group by sic2, fyear;
quit;

*---industry measures for ms----;

proc sort data=data2;
	by fyear sic2;
run;

proc univariate data=data2 noprint;
	by fyear sic2;
	var roa cfroa xrdint capxint xadint;
	output out=indmd median=md_roa md_cfroa md_xrdint md_capxint md_xadint;
run;

proc sql;
	create table data2 as select * from data2 a left join indmd b on 
		a.fyear=b.fyear and a.sic2=b.sic2;
quit;

proc sort data=data2 nodupkey;
	by gvkey datadate;
run;

data data2;
	set data2;
	*more for Mohanram score;

	if roa>md_roa then
		m1=1;
	else
		m1=0;

	if cfroa>md_cfroa then
		m2=1;
	else
		m2=0;

	if oancf>ni then
		m3=1;
	else
		m3=0;

	if xrdint>md_xrdint then
		m4=1;
	else
		m4=0;

	if capxint>md_capxint then
		m5=1;
	else
		m5=0;

	if xadint>md_xadint then
		m6=1;
	else
		m6=0;
	*still need to add another thing for Mohanram (2005) score;
run;

*----add credit rating--------;

proc sql;
	create table data2 as select a.*, b.splticrm from data2 a left join 
		compr.adsprate b on a.gvkey=b.gvkey and year(a.datadate)=year(b.datadate);
quit;

proc sort data=data2 nodupkey;
	by gvkey datadate;
run;

*consumer price index to create orgcap measure

						from Bureau of Labor Statistics website;

data cpi;
	infile datalines;
	input yr 4.0 cpi 10.3;
	datalines;
2017	246.19
2016	242.23
2015	236.53
2014	229.91
2013	229.17
2012	229.594
2011	224.939
2010	218.056
2009	214.537
2008	215.303
2007	207.342
2006	201.6
2005	195.3
2004	188.9
2003	183.96
2002	179.88
2001	177.1
2000	172.2
1999	166.6
1998	163.00
1997	160.5
1996	156.9
1995	152.4
1994	148.2
1993	144.5
1992	140.3
1991	136.2
1990	130.7
1989	124.00
1988	118.3
1987	113.6
1986	109.6
1985	107.6
1984	103.9
1983	99.6
1982	96.5
1981	90.9
1980	82.4
1979	72.6
1978	65.2
1977	60.6
1976	56.9
1975	53.8
1974  49.3
1973	44.2
1972	41.8
1971	40.6
1970	38.9
1969	36.7
1968	34.8
1967	33.4
1966	32.5
1965	31.6
1964	31.0
1963	30.7
1962	30.2
1961	29.9
1960	29.6
1959	29.2
1958	28.9
1957	28.2
1956	27.3
1955	26.8
1954	26.9
1953	26.8
1952	26.7
1951	25.9
1950	24.0
1949	23.7
1948	24.4
1947	22.2
1946	19.7
1945	18.1
1944	17.6
1943	17.3
1942	16.3
1941	14.7
1940	14.0
1939	13.8
1938	14.1
1937	14.4
1936	13.9
1935	13.6
1934	13.3
1933	13.1
1932	13.6
1931	15.1
1930	16.6
1929	17.2
1928	17.1
1927	17.2
1926	17.5
1925	17.7
1924	17.0
;
run;

proc sql;
	create table data2 as select a.*, b.cpi from data2 a left join cpi b on 
		a.fyear=b.yr;
quit;

proc sort data=data2 nodupkey;
	by gvkey datadate;
run;

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

/* check 4 */
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

proc sort data=data nodupkey;
	by gvkey datadate;
run;

proc sort data=data;
	by gvkey datadate;
run;

*create first set of quarterly compustat variables;

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

	/*HXZ*/
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
	*clean up;

	if first.gvkey then
		count=1;
	else
		count+1;

	if first.gvkey then
		do;
			roaq=.;
			roeq=.;
			kzq=.;
			laq=.;
			noaq=.;
			atoq=.;
			ctoq=.;
			glaq=.;
			oleq=.;
			olaq=.;
			claq=.;
		end;

	if count<3 then
		do;
			rnaq=.;
		end;

	if count<5 then
		do;
			chtx=.;
			che=.;
			cinvest=.;
			roeq=.;
			roaq=.;
			sgq=.;
		end;

	if count<5 /*count<6*/
	then
		do;
			rs=.;
		end;

	if count<5 /*count<17*/
	then
		do;
			stdacc=.;
			stdcf=.;
			sgrvol=.;
			roavol=.;
		end;
run;

/*HXZ industry adjustments*/
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

*finally finish Mohanram score components;

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
	else
		m8=0;
run;

/*  to get earnings forecasts from I/B/E/S that matches to quarterly Compustat variables  */
proc sql;
	create table ibessum as select ticker, cusip, fpedats, statpers, ANNDATS_ACT, 
		numest, ANNTIMS_ACT, medest, actual, stdev from ibes.statsum_epsus where 
		fpi='6'  /*1 is for annual forecasts, 6 is for quarterly*/
		and statpers<ANNDATS_ACT /*only keep summarized forecasts prior to earnings annoucement*/
		and measure='EPS' and not missing(medest) and not missing(fpedats) 
		and (fpedats-statpers)>=0;
quit;

*forecasts closest prior to fiscal period end date;

proc sort data=ibessum;
	by cusip fpedats descending statpers;
run;

proc sort data=ibessum nodupkey;
	by cusip fpedats;
run;

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

*get permno for CRSP data;

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

*=============================================

		Some of the RPS require daily CRSP data in conjunction with Compustat quarterly,
		so add daily CRSP info to create these RPS


	==============================================;
*this is for abnormal trading volume and returns around earings announcements;

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

*================================================================================


		First Compustat quarterly data set

================================================================================;

data data6;
	set data6;
	keep gvkey permno datadate datadate_q rdq chtx roaq rsup stdacc stdcf sgrvol 
		roavol cash cinvest nincr

		/*HXZ*/
		rdmq rdsq olq tanq kzq alaq almq laq bmq dmq amq epq cpq emq spq ndpq ebpq 
		wwq rs droeq droaq noaq rnaq pmq atoq ctoq glaq oleq olaq claq blq sgq sue 
		aeavol ear m7 m8 prccq roeq;
	where datadate>='01Jan2015'd;
run;

*==============================================================================

	add quarterly compustat data to monthly returns and annual compustat data


	===============================================================================;

proc sql;
	create table temp3 as select * from temp2 a left join data6 b on 
		a.permno=b.permno and intnx('MONTH', a.date, -10)<=b.datadate<=intnx('MONTH', 
		a.date, -5) and intnx('WEEKDAY', b.rdq, 1)<intnx('MONTH', a.date, -5, 'beg');
	*allow at least four months for quarterly info to become available;
	*date is the end of the return month;
quit;

*since some of the information is also using rdq, need to make sure that the monthly info is after that;

proc sort data=temp3;
	by permno date descending datadate;
run;

proc sort data=temp3 nodupkey;
	by permno date;
run;

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
quit;

proc sort data=temp4;
	by permno date descending statpers;
run;

proc sort data=temp4 nodupkey;
	by permno date;
run;

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
run;